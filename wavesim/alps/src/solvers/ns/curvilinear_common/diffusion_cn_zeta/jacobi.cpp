#include "jacobi.h"

#include "reduced_common.h"
#include <common/base/logging.h>
#include <common/container/view_utils.h>
#include <decomp/ghost_cell_exchange.h>
#include <solvers/ns/curvilinear_common/diffusion_cn_zeta_coeff.h>
#include <spectral/spectral.h>

#include <Kokkos_Core.hpp>
#include <mpipp/collectives.h>

#include <cmath>
#include <limits>
#include <stdexcept>

namespace alps::solver {

template<typename T>
using host_buf_t = MDView<T, alps::PoolSpace<Kokkos::SharedHostPinnedSpace>>;

namespace {
/** Performs one Jacobi sweep and returns the local squared update norm. */
double jacobi_iteration(MDView<Real***> const& first_line,
                        MDView<Real**> const&  dd0,
                        MDView<Real**> const&  x_prev,
                        MDView<Real**> const&  x_next,
                        bool                   has_prev,
                        bool                   has_next,
                        LoopPolicy<2> const&   policy)
{
  double norm = 0;
  Kokkos::parallel_reduce(
    "jacobi_iter",
    policy,
    KOKKOS_LAMBDA(int i, int j, double& val) {
      Real x_new = dd0(i, j);
      if (has_prev) {
        auto aa = first_line(i, j, 0);
        x_new -= aa * x_prev(i, j);
      }
      if (has_next) {
        auto cc = first_line(i, j, 1);
        x_new -= cc * x_next(i, j);
      }
      auto diff           = double(x_new - first_line(i, j, 2));
      first_line(i, j, 2) = x_new;
      val += diff * diff;
    },
    norm);
  return norm;
}

void exchange_sol(MDView<Real**> const&                sol,
                  MDView<Real**> const&                recv_prev,
                  MDView<Real**> const&                recv_next,
                  host_buf_t<Real***> const&           host_buf,
                  bool                                 has_prev,
                  bool                                 has_next,
                  Kokkos::DefaultExecutionSpace const& stream,
                  mpipp::communicator const&           comm)
{
  int constexpr tag{DistTridiagMpiTag::Reduced};

  int const  rank      = comm.rank();
  auto const prev_rank = has_prev ? rank - 1 : MPI_PROC_NULL;
  auto const next_rank = has_next ? rank + 1 : MPI_PROC_NULL;

  mpipp::irequest_pool requests;
  auto*                recv_prev_ptr = recv_prev.data();
  auto*                recv_next_ptr = recv_next.data();
  auto*                send_ptr      = sol.data();
  auto                 msg_size      = recv_prev.span();
  if constexpr (!alps::mpi_can_access_v<default_memory_pool>) {
    recv_prev_ptr = &host_buf(0, 0, 0);
    recv_next_ptr = &host_buf(0, 0, 1);
    if (has_next || has_prev) {
      Kokkos::deep_copy(stream, trailing_subview(host_buf, 2), sol);
      stream.fence();
    }
    send_ptr = &host_buf(0, 0, 2);
  }

  requests.push(
    mpipp::irecv(nonstd::span(recv_prev_ptr, msg_size), prev_rank, tag, comm));
  requests.push(
    mpipp::isend(nonstd::span(send_ptr, msg_size), prev_rank, tag, comm));
  requests.push(
    mpipp::irecv(nonstd::span(recv_next_ptr, msg_size), next_rank, tag, comm));
  requests.push(
    mpipp::isend(nonstd::span(send_ptr, msg_size), next_rank, tag, comm));
  requests.waitall();

  if constexpr (!alps::mpi_can_access_v<default_memory_pool>) {
    if (has_prev) {
      Kokkos::deep_copy(stream, recv_prev, trailing_subview(host_buf, 0));
    }
    if (has_next) {
      Kokkos::deep_copy(stream, recv_next, trailing_subview(host_buf, 1));
    }
    stream.fence();
  }
}

} // namespace

template<typename CoeffFunctor>
void solve_diffusion_cn_zeta_jacobi(MDView<Real***> const&              dst,
                                    MDView<Real const***> const&        rhs,
                                    CoeffFunctor const&                 coeff,
                                    int                                 nx,
                                    int                                 ny,
                                    int                                 n_eqns,
                                    DiffusionCNZetaSolverOptions const& options,
                                    mpipp::communicator const&          comm,
                                    Kokkos::DefaultExecutionSpace const& stream)
{
  auto const rank         = comm.rank();
  auto const nproc        = comm.size();
  auto const n_eqns_local = n_eqns;

  auto const create_view = [=](std::string name, int dim2) {
    return MDView<Real***, default_memory_pool>(
      Kokkos::view_alloc(name, Kokkos::WithoutInitializing), nx, ny, dim2);
  };
  auto aa         = create_view("cn_zeta_aa", n_eqns_local);
  auto cc         = create_view("cn_zeta_cc", n_eqns_local);
  auto first_line = create_view("cn_zeta_first", 3);
  auto last_line  = create_view("cn_zeta_last", 3);
  auto mpi_buf    = create_view("cn_zeta_jacobi_recv_buf", 3);

  using local_elim_t = DistTridiagLocalElimination<false, CoeffFunctor>;
  auto const local_elim =
    local_elim_t(dst, rhs, coeff, aa, cc, first_line, last_line, comm);
  Kokkos::parallel_for(
    "cn_zeta_forward_reduce",
    GridPolicy<typename local_elim_t::Forward>(stream, ny, Kokkos::AUTO()),
    local_elim);

  // Send each rank's last reduced equation to rank + 1. For rank r > 0, combine
  // rank r's first equation
  //   a_r x_{r-1,last} + x_{r,first} + c_r x_{r+1,first} = d_r
  // with rank r - 1's last equation
  //   A_{r-1} x_{r-1,first} + x_{r-1,last} + C_{r-1} x_{r,first} = D_{r-1}
  // to eliminate x_{r-1,last}. The retained couplings are x_{r-1,first}, x_{r,
  // first} and x_{r+1,first}, except on rank 0.
  stream.fence();
  alps::exchange(last_line,
                 rank + 1 < nproc ? rank + 1 : MPI_PROC_NULL,
                 mpi_buf,
                 rank > 0 ? rank - 1 : MPI_PROC_NULL,
                 comm,
                 DistTridiagMpiTag::Initial);
  auto policy = LoopPolicy<2>(
    stream, {0, 0}, {first_line.extent_int(0), first_line.extent_int(1)});
  if (rank > 0) {
    Kokkos::parallel_for(
      "cn_zeta_jacobi_initial_collapse", policy, KOKKOS_LAMBDA(int i, int j) {
        Real am1 = mpi_buf(i, j, 0);
        Real cm1 = mpi_buf(i, j, 1);
        Real dm1 = mpi_buf(i, j, 2);

        auto aa0 = first_line(i, j, 0);
        auto cc0 = first_line(i, j, 1);
        auto dd0 = first_line(i, j, 2);

        Real const bbi      = Real(1) / (Real(1) - aa0 * cm1);
        first_line(i, j, 2) = bbi * (dd0 - aa0 * dm1);
        first_line(i, j, 0) = bbi * -aa0 * am1;
        first_line(i, j, 1) = cc0 * bbi;
      });
  }

  bool   converged  = false;
  bool   diverged   = false;
  double norm0      = -1;
  double norm       = -1;
  int    iter_count = 0;
  auto   logger     = get_logger("cn_eqn");

  // reuse slot 0 of last_line to cache the pre-Jacobi first_sol value
  auto dd0 = trailing_subview(last_line, 0);
  Kokkos::deep_copy(stream, dd0, local_elim.first_sol);

  host_buf_t<Real***> host_buf;
  if constexpr (!alps::mpi_can_access_v<default_memory_pool>) {
    host_buf =
      create_mirror_view(Kokkos::WithoutInitializing,
                         alps::PoolSpace<Kokkos::SharedHostPinnedSpace>(),
                         mpi_buf);
  }
  auto recv_prev = trailing_subview(mpi_buf, 0);
  auto recv_next = trailing_subview(mpi_buf, 1);
  stream.fence();
  for (int iter = 0; iter < options.jacobi_max_iters; ++iter) {
    iter_count        = iter + 1;
    double local_norm = 0;
    if (rank > 0) {
      stream.fence();

      bool const has_prev = rank > 1;
      bool const has_next = rank + 1 < nproc;

      exchange_sol(local_elim.first_sol,
                   recv_prev,
                   recv_next,
                   host_buf,
                   has_prev,
                   has_next,
                   stream,
                   comm);

      local_norm = jacobi_iteration(
        first_line, dd0, recv_prev, recv_next, has_prev, has_next, policy);
    }

    mpipp::allreduce(local_norm, norm, mpipp::plus<double>(), comm);
    if (!std::isfinite(norm) || norm < 0) {
      diverged = true;
      logger->debug(
        "Jacobi reduced system iter {}: diverged with non-finite norm", iter);
      break;
    }
    norm = std::sqrt(norm / double(nproc));
    if (norm0 < 0) {
      norm0 = norm;
    }
    logger->debug(
      "Jacobi reduced system iter {}: norm={}, abs_tol={}, rel_norm={}, "
      "rel_tol={}",
      iter,
      norm,
      options.jacobi_abs_tol,
      norm / norm0,
      options.jacobi_rel_tol);
    if (norm <= options.jacobi_abs_tol
        || norm / norm0 <= options.jacobi_rel_tol) {
      converged = true;
      break;
    }
  }

  if (diverged) {
    throw std::runtime_error(
      "DiffusionCNZetaEqn Jacobi solver diverged (non-finite norm)");
  }

  if (!converged) {
    logger->warn(
      "Jacobi reduced system did not reach tolerance after {} iterations "
      "(norm={}, rel_norm={}, desired abs_tol={}, desired rel_tol={})",
      iter_count,
      norm,
      (norm0 > 0 ? norm / norm0 : 0.0),
      options.jacobi_abs_tol,
      options.jacobi_rel_tol);
  }

  stream.fence(); // ensure all Jacobi GPU work is visible before MPI sends
  alps::exchange(local_elim.first_sol,
                 rank > 0 ? rank - 1 : MPI_PROC_NULL,
                 local_elim.next_sol,
                 rank + 1 < nproc ? rank + 1 : MPI_PROC_NULL,
                 comm,
                 DistTridiagMpiTag::Backward);

  Kokkos::parallel_for(
    "cn_zeta_backward",
    GridPolicy<typename local_elim_t::Backward>(stream, ny, Kokkos::AUTO()),
    local_elim);
  stream.fence();
}

template void solve_diffusion_cn_zeta_jacobi<
  DiffusionCNZetaCoeffProvider<CenterPt, BottomWaveMesh>>(
  MDView<Real***> const&,
  MDView<Real const***> const&,
  DiffusionCNZetaCoeffProvider<CenterPt, BottomWaveMesh> const&,
  int,
  int,
  int,
  DiffusionCNZetaSolverOptions const&,
  mpipp::communicator const&,
  Kokkos::DefaultExecutionSpace const&);

template void solve_diffusion_cn_zeta_jacobi<
  DiffusionCNZetaCoeffProvider<NodePt, BottomWaveMesh>>(
  MDView<Real***> const&,
  MDView<Real const***> const&,
  DiffusionCNZetaCoeffProvider<NodePt, BottomWaveMesh> const&,
  int,
  int,
  int,
  DiffusionCNZetaSolverOptions const&,
  mpipp::communicator const&,
  Kokkos::DefaultExecutionSpace const&);

} // namespace alps::solver
