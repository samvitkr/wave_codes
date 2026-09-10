//
// Created by xuananqing on 4/1/23.
//

#include "mean_variance.h"

#include <common/base/macros.h>
#include <common/kokkos_abstraction/pool_space.h>
#include <common/real_type.h>

#include <Kokkos_Core.hpp>
#include <mpipp/collectives.h>

namespace alps {

template<typename T, typename MemSpace>
MeanVariance<T, MemSpace>::MeanVariance(
  Kokkos::View<T const***, Kokkos::LayoutLeft, MemSpace> const& f_)
  : f(f_)
  , mean(f.label() + "_mean", f_.extent(0), f_.extent(2))
  , M2(f.label() + "_M2", f_.extent(0), f_.extent(2))
{}

template<typename T, typename MemSpace>
MeanVariance<T, MemSpace>::MeanVariance(
  Kokkos::View<T**, Kokkos::LayoutLeft, Kokkos::HostSpace> const& mean_,
  Kokkos::View<T**, Kokkos::LayoutLeft, Kokkos::HostSpace> const& M2_,
  Kokkos::View<T const***, Kokkos::LayoutLeft, MemSpace> const&   f_)
  : f(f_)
  , mean(mean_)
  , M2(M2_)
{
  if (f.extent(0) != mean.extent(0) || f.extent(2) != mean.extent(1)
      || f.extent(0) != M2.extent(0) || f.extent(2) != M2.extent(1)) {
    throw std::invalid_argument("MeanVariance::MeanVariance(): Dimensions of "
                                "f, mean and M2 do not match.");
  }
}

template<typename T, typename MemSpace>
KOKKOS_FUNCTION void
MeanVariance<T, MemSpace>::operator()(member_t const& team) const
{
  using Kokkos::parallel_for;
  using Kokkos::TeamVectorRange;
  const auto k = team.league_rank();

  /*
   * Compute mean and M2 using one-pass online algorithm
   * see https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
   */
  for (int j = 0; j < f.extent_int(1); ++j) {
    parallel_for(
      TeamVectorRange(team, f.extent_int(0)), KOKKOS_TR_LAMBDA(int i) {
        const auto val       = f(i, j, k);
        const auto delta     = val - mean_d(i, k);
        const auto mean_new  = mean_d(i, k) + delta / (j + 1);
        const auto delta_new = val - mean_new;

        mean_d(i, k) = mean_new;
        M2_d(i, k) += delta * delta_new;
      });
  }
}

template<typename T, typename MemSpace>
void MeanVariance<T, MemSpace>::calculate_local()
{
  mean_d = Kokkos::create_mirror_view(MemSpace(), mean);
  M2_d   = Kokkos::create_mirror_view(MemSpace(), M2);
  Kokkos::deep_copy(mean_d, 0);
  Kokkos::deep_copy(M2_d, 0);

  policy_t const policy(
    policy_t::execution_space(), f.extent_int(2), Kokkos::AUTO());
  Kokkos::parallel_for("mean_variance", policy, *this);

  Kokkos::deep_copy(policy.space(), mean, mean_d);
  Kokkos::deep_copy(policy.space(), M2, M2_d);

  policy.space().fence();

  // cleanup
  mean_d = {};
  M2_d   = {};
}

template<typename T, typename MemSpace>
void MeanVariance<T, MemSpace>::gather(int                        root,
                                       const mpipp::communicator& comm)
{
  using gather_view_t =
    Kokkos::View<Real***, Kokkos::LayoutLeft, Kokkos::HostSpace>;
  const int  ny      = f.extent_int(1);
  const bool is_root = (comm.rank() == root);
  const int  n_procs = comm.size();

  gather_view_t const mean_gather("mean_gather",
                                  is_root ? mean.extent(0) : 1,
                                  is_root ? mean.extent(1) : 1,
                                  n_procs);
  gather_view_t const M2_gather("var_gather", mean_gather.layout());

  mpipp::gather(
    nonstd::span(mean.data(), mean.span()), mean_gather.data(), root, comm);
  mpipp::gather(
    nonstd::span(M2.data(), M2.span()), M2_gather.data(), root, comm);

  if (!is_root) {
    mpipp::barrier(comm);
    return;
  }

  Kokkos::deep_copy(mean, 0);
  Kokkos::deep_copy(M2, 0);

  /** accumulate the mean and variance
   * Numerically Stable Parallel Computation of (Co-)Variance, Schubert & Gertz
   */
  const int ny_total = ny * n_procs;
#pragma omp parallel for
  for (int k = 0; k < mean.extent_int(1); ++k) {
    for (int p = 0; p < n_procs; ++p) {
#pragma omp simd
      for (int i = 0; i < mean.extent_int(0); ++i) {
        auto const delta = mean_gather(i, k, p) - mean(i, k);
        mean(i, k) += delta / (T)(p + 1);
        auto const delta_new = mean_gather(i, k, p) - mean(i, k);
        M2(i, k) += M2_gather(i, k, p) + delta * delta_new * ny;
      }
    }
#pragma omp simd
    for (int i = 0; i < mean.extent_int(0); ++i) {
      M2(i, k) /= ny_total;
    }
  }

  mpipp::barrier(comm);
}

template struct MeanVariance<Real, default_memory_pool>;
} // namespace alps
