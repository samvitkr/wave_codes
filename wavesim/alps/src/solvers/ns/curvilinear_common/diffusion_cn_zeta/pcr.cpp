#include "pcr.h"

#include "reduced_common.h"
#include <common/kokkos_abstraction/pool_space.h>
#include <decomp/ghost_cell_exchange.h>
#include <solvers/mesh/curvilinear_mesh.h>
#include <solvers/ns/curvilinear_common/diffusion_cn_zeta_coeff.h>
#include <spectral/spectral.h>

#include <Kokkos_Core.hpp>

#include <stdexcept>

namespace alps::solver {

namespace {

/**
 * Performs one parallel cyclic reduction update on active equations.
 * When has_prev or has_next is false, the corresponding neighbor is treated as
 * absent and its received values are not read.
 */
void pcr_iteration(MDView<Real***> const&       active,
                   MDView<Real const***> const& prev,
                   MDView<Real const***> const& next,
                   bool                         has_prev,
                   bool                         has_next,
                   LoopPolicy<2> const&         policy)
{
  Kokkos::parallel_for(
    "pcr", policy, KOKKOS_LAMBDA(int i, int j) {
      Real cm1        = has_prev ? prev(i, j, 1) : 0;
      Real dm1        = has_prev ? prev(i, j, 2) : 0;
      Real ap1        = has_next ? next(i, j, 0) : 0;
      Real dp1        = has_next ? next(i, j, 2) : 0;
      auto aa0        = active(i, j, 0);
      auto cc0        = active(i, j, 1);
      auto dd0        = active(i, j, 2);
      Real bbi        = Real(1) / (Real(1) - aa0 * cm1 - cc0 * ap1);
      active(i, j, 2) = bbi * (dd0 - aa0 * dm1 - cc0 * dp1);
      if (has_prev) {
        active(i, j, 0) = -bbi * aa0 * prev(i, j, 0);
      }
      if (has_next) {
        active(i, j, 1) = -bbi * cc0 * next(i, j, 1);
      }
    });
}

/** Exchanges one contiguous buffer with both neighboring ranks. */
void bidirection_exchange(MDView<Real***> const&               send,
                          int                                  prev_rank,
                          int                                  next_rank,
                          MDView<Real***> const&               recv_prev,
                          MDView<Real***> const&               recv_next,
                          int                                  tag,
                          mpipp::communicator const&           comm,
                          Kokkos::DefaultExecutionSpace const& stream)
{
  bool const has_prev = prev_rank != MPI_PROC_NULL;
  bool const has_next = next_rank != MPI_PROC_NULL;
  bool constexpr is_gpu_direct =
    alps::mpi_can_access_v<Kokkos::DefaultExecutionSpace::memory_space>;

  auto const create_host_pinned_mirror = [=](MDView<Real***> const& view)
    -> MDView<Real***, Kokkos::SharedHostPinnedSpace> {
    if constexpr (alps::mpi_can_access_v<
                    Kokkos::DefaultExecutionSpace::memory_space>) {
      return {};
    } else {
      return create_mirror_view(
        Kokkos::WithoutInitializing,
        alps::PoolSpace<Kokkos::SharedHostPinnedSpace>(),
        view);
    }
  };
  auto send_host      = create_host_pinned_mirror(send);
  auto recv_prev_host = create_host_pinned_mirror(recv_prev);
  auto recv_next_host = create_host_pinned_mirror(recv_next);

  auto* send_ptr = send.data();
  if constexpr (!is_gpu_direct) {
    Kokkos::deep_copy(stream, send_host, send);
    send_ptr = send_host.data();
  }

  mpipp::irequest_pool requests;
  if (has_prev) {
    auto* ptr  = is_gpu_direct ? recv_prev.data() : recv_prev_host.data();
    auto  span = is_gpu_direct ? recv_prev.span() : recv_prev_host.span();
    requests.push(mpipp::irecv(nonstd::span(ptr, span), prev_rank, tag, comm));
  }
  if (has_next) {
    auto* ptr  = is_gpu_direct ? recv_next.data() : recv_next_host.data();
    auto  span = is_gpu_direct ? recv_next.span() : recv_next_host.span();
    requests.push(mpipp::irecv(nonstd::span(ptr, span), next_rank, tag, comm));
  }

  stream.fence(); // ensure GPU writes to send buffer are visible before sends
  if (has_prev) {
    requests.push(
      mpipp::isend(nonstd::span(send_ptr, send.span()), prev_rank, tag, comm));
  }
  if (has_next) {
    requests.push(
      mpipp::isend(nonstd::span(send_ptr, send.span()), next_rank, tag, comm));
  }
  requests.waitall();

  if constexpr (!is_gpu_direct) {
    if (has_prev) {
      Kokkos::deep_copy(stream, recv_prev, recv_prev_host);
    }
    if (has_next) {
      Kokkos::deep_copy(stream, recv_next, recv_next_host);
    }
    stream.fence();
  }
}

} // namespace

template<typename CoeffFunctor>
void solve_diffusion_cn_zeta_pcr(MDView<Real***> const&               dst,
                                 MDView<Real const***> const&         rhs,
                                 CoeffFunctor const&                  coeff,
                                 int                                  nx,
                                 int                                  ny,
                                 int                                  n_eqns,
                                 mpipp::communicator const&           comm,
                                 Kokkos::DefaultExecutionSpace const& stream)
{
  auto const rank         = comm.rank();
  auto const nproc        = comm.size();
  auto const n_eqns_local = n_eqns;

  if (n_eqns_local < 2) {
    throw std::invalid_argument(
      "DiffusionCNZetaEqn PCR requires at least two local z equations");
  }

  auto const create_view = [=](std::string name, int dim2) {
    return MDView<Real***, default_memory_pool>(
      Kokkos::view_alloc(name, Kokkos::WithoutInitializing), nx, ny, dim2);
  };
  auto        aa         = create_view("cn_zeta_aa", n_eqns_local);
  auto        cc         = create_view("cn_zeta_cc", n_eqns_local);
  auto        first_line = create_view("cn_zeta_first", 3);
  auto        last_line  = create_view("cn_zeta_last", 3);
  auto        recv_prev  = create_view("cn_zeta_recv_prev", 3);
  auto const& recv_next  = last_line; // reuse last_line

  using local_elim_t = DistTridiagLocalElimination<true, CoeffFunctor>;
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
                 recv_prev,
                 rank > 0 ? rank - 1 : MPI_PROC_NULL,
                 comm,
                 DistTridiagMpiTag::Initial);
  auto policy =
    LoopPolicy<2>(stream, {0, 0}, {first_line.extent(0), first_line.extent(1)});
  if (rank > 0) {
    Kokkos::parallel_for(
      "pcr_initial_collapse", policy, KOKKOS_LAMBDA(int i, int j) {
        Real am1 = recv_prev(i, j, 0);
        Real cm1 = recv_prev(i, j, 1);
        Real dm1 = recv_prev(i, j, 2);

        auto aa0 = first_line(i, j, 0);
        auto cc0 = first_line(i, j, 1);
        auto dd0 = first_line(i, j, 2);

        Real const bbi      = Real(1) / (Real(1) - aa0 * cm1);
        first_line(i, j, 2) = bbi * (dd0 - aa0 * dm1);
        first_line(i, j, 0) = bbi * -aa0 * am1;
        first_line(i, j, 1) = cc0 * bbi;
      });
  }

  int const p_iters = std::max(0, int(std::ceil(std::log2(double(nproc)))));
  for (int p = 0, s = 1; p < p_iters; ++p, s <<= 1) {
    int const  rank_ms  = rank - s;
    int const  rank_ps  = rank + s;
    bool const has_prev = rank_ms >= 0;
    bool const has_next = rank_ps < nproc;
    if (!has_prev && !has_next) {
      continue;
    }
    bidirection_exchange(first_line,
                         has_prev ? rank_ms : MPI_PROC_NULL,
                         has_next ? rank_ps : MPI_PROC_NULL,
                         recv_prev,
                         recv_next,
                         DistTridiagMpiTag::Reduced,
                         comm,
                         stream);
    pcr_iteration(first_line, recv_prev, recv_next, has_prev, has_next, policy);
  }

  stream.fence(); // ensure PCR results are visible before MPI sends
  alps::exchange(local_elim.first_sol,
                 rank > 0 ? rank - 1 : MPI_PROC_NULL,
                 local_elim.next_sol, // this is an alias of last_line
                 rank + 1 < nproc ? rank + 1 : MPI_PROC_NULL,
                 comm,
                 DistTridiagMpiTag::Backward);

  Kokkos::parallel_for(
    "cn_zeta_backward",
    GridPolicy<typename local_elim_t::Backward>(stream, ny, Kokkos::AUTO()),
    local_elim);
  stream.fence();
}

template void solve_diffusion_cn_zeta_pcr<
  DiffusionCNZetaCoeffProvider<CenterPt, BottomWaveMesh>>(
  MDView<Real***> const&,
  MDView<Real const***> const&,
  DiffusionCNZetaCoeffProvider<CenterPt, BottomWaveMesh> const&,
  int,
  int,
  int,
  mpipp::communicator const&,
  Kokkos::DefaultExecutionSpace const&);

template void solve_diffusion_cn_zeta_pcr<
  DiffusionCNZetaCoeffProvider<NodePt, BottomWaveMesh>>(
  MDView<Real***> const&,
  MDView<Real const***> const&,
  DiffusionCNZetaCoeffProvider<NodePt, BottomWaveMesh> const&,
  int,
  int,
  int,
  mpipp::communicator const&,
  Kokkos::DefaultExecutionSpace const&);

} // namespace alps::solver
