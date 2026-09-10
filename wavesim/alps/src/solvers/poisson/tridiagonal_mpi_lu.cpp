#include "tridiagonal_mpi_lu.h"

#include <common/base/macros.h>
#include <common/device/device_traits.h>
#include <common/kokkos_abstraction/exec_policy.h>
#include <common/kokkos_abstraction/pool_space.h>
#include <common/runtime/async_utils.h>

#include <Kokkos_Core.hpp>
#include <mpipp/comm.h>
#include <mpipp/point2point.h>

#include <stdexcept>

namespace alps {
namespace solver {

template<class DataType, class... Properties>
using MDView = Kokkos::View<DataType, Kokkos::LayoutLeft, Properties...>;
using Kokkos::ALL;
using Kokkos::subview;

template<typename ValueT>
TridiagonalMpiLU<ValueT>::TridiagonalMpiLU(int batch_count_1,
                                           int batch_count_2,
                                           int n_eqns,
                                           const mpipp::communicator& comm)
  : base_t(batch_count_1, batch_count_2, n_eqns, comm)
{}

template<typename ValueT>
void TridiagonalMpiLU<ValueT>::setup_impl(coeff_t const& d,
                                          coeff_t const& dl,
                                          coeff_t const& du)
{
  auto n1 = this->n1_;
  auto n2 = this->n2_;
  auto n3 = this->nz_;

  if (d.extent_int(0) < n1 || d.extent_int(1) < n2 || d.extent_int(2) != n3) {
    throw std::invalid_argument("Poisson solver coefficient d size mismatch.");
  }
  if (dl.extent_int(0) < n1 || dl.extent_int(1) < n2
      || dl.extent_int(2) != n3) {
    throw std::invalid_argument("Poisson solver coefficient dl size mismatch.");
  }
  if (du.extent_int(0) < n1 || du.extent_int(1) < n2
      || du.extent_int(2) != n3) {
    throw std::invalid_argument("Poisson solver coefficient du size mismatch.");
  }

  auto policy2d = LoopPolicy<2>(get_next_stream(), {0, 0}, {n1, n2});

  MDView<ValueT***, default_memory_pool> const recvbuf(
    "MPI LU recv buffer", d.extent(0), d.extent(1), 2);
  MDView<ValueT***, default_memory_pool> const sendbuf(
    "MPI LU send buffer", d.extent(0), d.extent(1), 2);

  if (this->rank_ != 0) {
    if constexpr (!alps::mpi_can_access_v<default_memory_pool>) {
      auto recvbuf_host =
        create_mirror_view(Kokkos::WithoutInitializing,
                           alps::PoolSpace<Kokkos::SharedHostPinnedSpace>(),
                           recvbuf);
      mpipp::recv(nonstd::span(recvbuf_host.data(), recvbuf.span()),
                  this->rank_ - 1,
                  0,
                  this->comm_,
                  mpipp::status_ignore);
      deep_copy(policy2d.space(), recvbuf, recvbuf_host);
      policy2d.space().fence(); // wait for the copy to complete
    } else {
      mpipp::recv(nonstd::span(recvbuf.data(), recvbuf.span()),
                  this->rank_ - 1,
                  0,
                  this->comm_,
                  mpipp::status_ignore);
    }
    Kokkos::parallel_for(
      policy2d, KOKKOS_LAMBDA(int i, int j) {
        dl(i, j, 0) /= recvbuf(i, j, 0);
        d(i, j, 0) -= dl(i, j, 0) * recvbuf(i, j, 1);
      });
  }

  using policy_t = Kokkos::TeamPolicy<Kokkos::IndexType<int>>;
  using member_t = policy_t::member_type;
  using Kokkos::AUTO;
  using Kokkos::parallel_for;
  using Kokkos::TeamThreadRange;
  auto policy_n2 = policy_t(policy2d.space(), n2, AUTO());
  parallel_for(
    policy_n2, KOKKOS_LAMBDA(member_t team) {
      int j = team.league_rank();
      for (int k = 1; k < n3; ++k) {
        parallel_for(
          TeamThreadRange(team, n1), KOKKOS_TR_LAMBDA(int& i) {
            dl(i, j, k) = dl(i, j, k) / d(i, j, k - 1);
            d(i, j, k)  = d(i, j, k) - dl(i, j, k) * du(i, j, k - 1);
          });
      }
    });

  if (this->rank_ != this->nproc_ - 1) {
    Kokkos::parallel_for(
      policy2d, KOKKOS_LAMBDA(int i, int j) {
        sendbuf(i, j, 0) = d(i, j, n3 - 1);
        sendbuf(i, j, 1) = du(i, j, n3 - 1);
      });
    if constexpr (!alps::mpi_can_access_v<default_memory_pool>) {
      auto sendbuf_host =
        create_mirror_view(Kokkos::WithoutInitializing,
                           alps::PoolSpace<Kokkos::SharedHostPinnedSpace>(),
                           sendbuf);
      deep_copy(policy2d.space(), sendbuf_host, sendbuf);
      policy2d.space().fence();
      mpipp::send(nonstd::span(sendbuf_host.data(), sendbuf.span()),
                  this->rank_ + 1,
                  0,
                  this->comm_);
    } else {
      policy2d.space().fence();
      mpipp::send(nonstd::span(sendbuf.data(), sendbuf.span()),
                  this->rank_ + 1,
                  0,
                  this->comm_);
    }
  }
}

template<typename ValueT>
void TridiagonalMpiLU<ValueT>::solve_impl(solution_t const& x,
                                          coeff_t const&    d,
                                          coeff_t const&    dl,
                                          coeff_t const&    du)
{
  auto n1 = this->n1_;
  auto n2 = this->n2_;
  auto n3 = this->nz_;

  auto policy2d = LoopPolicy<2>(get_next_stream(), {0, 0}, {n1, n2});

  MDView<ValueT**, default_memory_pool> const recvbuf(
    "LU recv buffer", d.extent(0), d.extent(1));

  // forward substitution
  // receive the solution from below
  if (this->rank_ != 0) {
    if constexpr (!alps::mpi_can_access_v<default_memory_pool>) {
      auto recvbuf_host =
        create_mirror_view(Kokkos::WithoutInitializing,
                           alps::PoolSpace<Kokkos::SharedHostPinnedSpace>(),
                           recvbuf);
      mpipp::recv(nonstd::span(recvbuf_host.data(), recvbuf.span()),
                  this->rank_ - 1,
                  0,
                  this->comm_,
                  mpipp::status_ignore);
      deep_copy(policy2d.space(), recvbuf, recvbuf_host);
      policy2d.space().fence(); // wait for the copy to complete
    } else {
      mpipp::recv(nonstd::span(recvbuf.data(), recvbuf.span()),
                  this->rank_ - 1,
                  0,
                  this->comm_,
                  mpipp::status_ignore);
    }
    Kokkos::parallel_for(
      policy2d, KOKKOS_LAMBDA(int i, int j) {
        x(i, j, 0) -= dl(i, j, 0) * recvbuf(i, j);
      });
  }

  using policy_t = Kokkos::TeamPolicy<Kokkos::IndexType<int>>;
  using member_t = policy_t::member_type;
  using Kokkos::AUTO;
  using Kokkos::parallel_for;
  using Kokkos::TeamThreadRange;

  auto policy_n2 = policy_t(policy2d.space(), n2, AUTO());
  // continue local substitution
  parallel_for(
    policy_n2, KOKKOS_LAMBDA(member_t team) {
      int j = team.league_rank();
      for (int k = 1; k < n3; ++k) {
        parallel_for(
          TeamThreadRange(team, n1), KOKKOS_TR_LAMBDA(int& i) {
            x(i, j, k) -= dl(i, j, k) * x(i, j, k - 1);
          });
      }
    });
  if (this->rank_ != this->nproc_ - 1) {
    auto x_top = subview(x, ALL, ALL, n3 - 1);
    if constexpr (!alps::mpi_can_access_v<default_memory_pool>) {
      auto x_top_host =
        create_mirror_view(Kokkos::WithoutInitializing,
                           alps::PoolSpace<Kokkos::SharedHostPinnedSpace>(),
                           x_top);
      deep_copy(policy_n2.space(), x_top_host, x_top);
      policy_n2.space().fence();
      mpipp::send(nonstd::span(x_top_host.data(), x_top_host.span()),
                  this->rank_ + 1,
                  0,
                  this->comm_);
    } else {
      policy_n2.space().fence();
      mpipp::send(nonstd::span(x_top.data(), x_top.span()),
                  this->rank_ + 1,
                  0,
                  this->comm_);
    }
  }

  // backward substitution
  if (this->rank_ == this->nproc_ - 1) {
    Kokkos::parallel_for(
      policy2d,
      KOKKOS_LAMBDA(int i, int j) { x(i, j, n3 - 1) /= d(i, j, n3 - 1); });
  } else {
    if constexpr (!alps::mpi_can_access_v<default_memory_pool>) {
      auto recvbuf_host =
        create_mirror_view(Kokkos::WithoutInitializing,
                           alps::PoolSpace<Kokkos::SharedHostPinnedSpace>(),
                           recvbuf);
      mpipp::recv(nonstd::span(recvbuf_host.data(), recvbuf.span()),
                  this->rank_ + 1,
                  1,
                  this->comm_,
                  mpipp::status_ignore);
      deep_copy(policy2d.space(), recvbuf, recvbuf_host);
      policy2d.space().fence(); // wait for the copy to complete
    } else {
      mpipp::recv(nonstd::span(recvbuf.data(), recvbuf.span()),
                  this->rank_ + 1,
                  1,
                  this->comm_,
                  mpipp::status_ignore);
    }
    Kokkos::parallel_for(
      policy2d, KOKKOS_LAMBDA(int i, int j) {
        x(i, j, n3 - 1) -= du(i, j, n3 - 1) * recvbuf(i, j);
        x(i, j, n3 - 1) /= d(i, j, n3 - 1);
      });
  }

  parallel_for(
    policy_n2, KOKKOS_LAMBDA(member_t team) {
      int j = team.league_rank();
      for (int k = n3 - 2; k >= 0; --k) {
        parallel_for(
          TeamThreadRange(team, n1), KOKKOS_TR_LAMBDA(int& i) {
            x(i, j, k) -= du(i, j, k) * x(i, j, k + 1);
            x(i, j, k) /= d(i, j, k);
          });
      }
    });

  if (this->rank_ != 0) {
    auto x_bot = subview(x, ALL, ALL, 0);
    if constexpr (!alps::mpi_can_access_v<default_memory_pool>) {
      auto x_bot_host =
        create_mirror_view(Kokkos::WithoutInitializing,
                           alps::PoolSpace<Kokkos::SharedHostPinnedSpace>(),
                           x_bot);
      deep_copy(policy_n2.space(), x_bot_host, x_bot);
      policy_n2.space().fence();
      mpipp::send(nonstd::span(x_bot_host.data(), x_bot_host.span()),
                  this->rank_ - 1,
                  1,
                  this->comm_);
    } else {
      policy_n2.space().fence();
      mpipp::send(nonstd::span(x_bot.data(), x_bot.span()),
                  this->rank_ - 1,
                  1,
                  this->comm_);
    }
  }

  policy_n2.space().fence();
}

template class TridiagonalMpiLU<double>;
template class TridiagonalMpiLU<float>;

} // namespace solver
} // namespace alps
