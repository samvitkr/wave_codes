#include "tridiagonal_single.h"

#include <common/base/macros.h>
#include <common/runtime/async_utils.h>

#include <Kokkos_Core.hpp>

#include <stdexcept>

namespace alps {
namespace solver {

template<typename ValueT>
TridiagonalSeqLU<ValueT>::TridiagonalSeqLU(int batch_count_1,
                                           int batch_count_2,
                                           int n_eqns,
                                           const mpipp::communicator& comm)
  : base_t(batch_count_1, batch_count_2, n_eqns, comm)
{
  if (this->nproc_ != 1) {
    throw std::runtime_error("The sequential LU solver supports only single "
                             "process.");
  }
}

template<typename ValueT>
void TridiagonalSeqLU<ValueT>::setup_impl(coeff_t const& d,
                                          coeff_t const& dl,
                                          coeff_t const& du)
{
  auto const n1 = this->n1_;
  auto const n2 = this->n2_;
  auto const n3 = this->nz_;

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

  using policy_t = Kokkos::TeamPolicy<Kokkos::IndexType<int>>;
  using member_t = policy_t::member_type;
  using Kokkos::AUTO;
  using Kokkos::parallel_for;
  using Kokkos::TeamThreadRange;

  auto stream = get_next_stream();
  parallel_for(
    policy_t(stream, n2, AUTO()), KOKKOS_LAMBDA(member_t team) {
      int j = team.league_rank();
      for (int k = 1; k < n3; ++k) {
        parallel_for(
          TeamThreadRange(team, n1), KOKKOS_TR_LAMBDA(int& i) {
            dl(i, j, k) = dl(i, j, k) / d(i, j, k - 1);
            d(i, j, k)  = d(i, j, k) - dl(i, j, k) * du(i, j, k - 1);
          });
        team.team_barrier();
      }
    });

  stream.fence();
}

template<typename ValueT>
void TridiagonalSeqLU<ValueT>::solve_impl(solution_t const& x,
                                          coeff_t const&    d,
                                          coeff_t const&    dl,
                                          coeff_t const&    du)
{
  auto const n1  = this->n1_;
  auto const n2  = this->n2_;
  auto const n3  = this->nz_;
  using policy_t = Kokkos::TeamPolicy<Kokkos::IndexType<int>>;
  using member_t = policy_t::member_type;
  using Kokkos::AUTO;
  using Kokkos::parallel_for;
  using Kokkos::TeamVectorRange;

  // forward substitution
  auto stream = get_next_stream();
  parallel_for(
    policy_t(stream, n2, AUTO()), KOKKOS_LAMBDA(member_t team) {
      const int j = team.league_rank();
      for (int k = 1; k < n3; ++k) {
        parallel_for(
          TeamVectorRange(team, n1), KOKKOS_TR_LAMBDA(int& i) {
            x(i, j, k) -= dl(i, j, k) * x(i, j, k - 1);
          });
        team.team_barrier();
      }
      parallel_for(
        TeamVectorRange(team, n1),
        KOKKOS_TR_LAMBDA(int& i) { x(i, j, n3 - 1) /= d(i, j, n3 - 1); });
    });

  // backward substitution
  parallel_for(
    policy_t(stream, n2, AUTO()), KOKKOS_LAMBDA(member_t team) {
      const int j = team.league_rank();
      for (int k = n3 - 2; k >= 0; --k) {
        parallel_for(
          TeamVectorRange(team, n1), KOKKOS_TR_LAMBDA(int& i) {
            x(i, j, k) =
              (x(i, j, k) - du(i, j, k) * x(i, j, k + 1)) / d(i, j, k);
          });
        team.team_barrier();
      }
    });

  stream.fence();
}

template class TridiagonalSeqLU<double>;
template class TridiagonalSeqLU<float>;

} // namespace solver
} // namespace alps
