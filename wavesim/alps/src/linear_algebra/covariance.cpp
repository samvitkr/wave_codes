//
// Created by xuananqing on 4/4/23.
//

#include "covariance.h"

#include <common/base/macros.h>
#include <common/kokkos_abstraction/pool_space.h>
#include <common/real_type.h>

#include <Kokkos_Core.hpp>
#include <mpipp/collectives.h>

namespace alps {

template<typename T, typename MemSpace>
Covariance<T, MemSpace>::Covariance(
  Kokkos::View<T const***, Kokkos::LayoutLeft, MemSpace> const& f_,
  Kokkos::View<T const***, Kokkos::LayoutLeft, MemSpace> const& g_)
  : f(f_)
  , g(g_)
  , mean_f(f.label() + "_mean", f_.extent(0), f_.extent(2))
  , mean_g(g.label() + "_mean", g_.extent(0), g_.extent(2))
  , M2(f.label() + g.label() + "_cov", f_.extent(0), f_.extent(2))
{
  if (f.layout() != g.layout()) {
    throw std::invalid_argument(
      "Covariance::Covariance(): Dimensions of f and g do not match.");
  }
}

template<typename T, typename MemSpace>
Covariance<T, MemSpace>::Covariance(
  const Kokkos::View<T**, Kokkos::LayoutLeft, Kokkos::HostSpace>& cov_,
  const Kokkos::View<const T***, Kokkos::LayoutLeft, MemSpace>&   f_,
  const Kokkos::View<const T***, Kokkos::LayoutLeft, MemSpace>&   g_)
  : f(f_)
  , g(g_)
  , mean_f("" + f_.label() + "_mean", f_.extent(0), f_.extent(2))
  , mean_g("" + g_.label() + "_mean", g_.extent(0), g_.extent(2))
  , M2(cov_)
{
  if (f.layout() != g.layout()) {
    throw std::invalid_argument(
      "Covariance::Covariance(): Dimensions of f and g do not match.");
  }
  if (mean_f.layout() != M2.layout()) {
    throw std::invalid_argument("Covariance::Covariance(): Dimensions of Cov do"
                                "not match f or g.");
  }
}

template<typename T, typename MemSpace>
KOKKOS_FUNCTION void
Covariance<T, MemSpace>::operator()(member_t const& team) const
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
        const auto val_f       = f(i, j, k);
        const auto delta_f     = val_f - mean_f_d(i, k);
        const auto mean_f_new  = mean_f_d(i, k) + delta_f / (j + 1);
        const auto delta_f_new = val_f - mean_f_new;

        const auto val_g      = g(i, j, k);
        const auto delta_g    = val_g - mean_g_d(i, k);
        const auto mean_g_new = mean_g_d(i, k) + delta_g / (j + 1);

        mean_f_d(i, k) = mean_f_new;
        mean_g_d(i, k) = mean_g_new;
        M2_d(i, k) += delta_f_new * delta_g;
      });
  }
}

template<typename T, typename MemSpace>
void Covariance<T, MemSpace>::calculate_local()
{
  mean_f_d = Kokkos::create_mirror_view(MemSpace(), mean_f);
  mean_g_d = Kokkos::create_mirror_view(MemSpace(), mean_g);
  M2_d     = Kokkos::create_mirror_view(MemSpace(), M2);
  Kokkos::deep_copy(mean_f_d, 0);
  Kokkos::deep_copy(mean_g_d, 0);
  Kokkos::deep_copy(M2_d, 0);

  policy_t const policy(
    policy_t::execution_space(), f.extent_int(2), Kokkos::AUTO());
  Kokkos::parallel_for("covariance", policy, *this);

  Kokkos::deep_copy(policy.space(), mean_f, mean_f_d);
  Kokkos::deep_copy(policy.space(), mean_g, mean_g_d);
  Kokkos::deep_copy(policy.space(), M2, M2_d);

  policy.space().fence();

  // cleanup
  mean_f_d = {};
  mean_g_d = {};
  M2_d     = {};
}

template<typename T, typename MemSpace>
void Covariance<T, MemSpace>::gather(int root, const mpipp::communicator& comm)
{
  using gather_view_t =
    Kokkos::View<Real***, Kokkos::LayoutLeft, Kokkos::HostSpace>;
  const int  ny      = f.extent_int(1);
  const bool is_root = (comm.rank() == root);
  const int  n_procs = comm.size();

  gather_view_t const mean_f_gather("mean_f_gather",
                                    is_root ? mean_f.extent(0) : 1,
                                    is_root ? mean_f.extent(1) : 1,
                                    n_procs);
  gather_view_t const mean_g_gather("mean_g_gather", mean_f_gather.layout());
  gather_view_t const M2_gather("cov_gather", mean_f_gather.layout());

  mpipp::gather(nonstd::span(mean_f.data(), mean_f.span()),
                mean_f_gather.data(),
                root,
                comm);
  mpipp::gather(nonstd::span(mean_g.data(), mean_g.span()),
                mean_g_gather.data(),
                root,
                comm);
  mpipp::gather(
    nonstd::span(M2.data(), M2.span()), M2_gather.data(), root, comm);

  if (!is_root) {
    comm.barrier();
    return;
  }

  Kokkos::deep_copy(mean_f, 0);
  Kokkos::deep_copy(mean_g, 0);
  Kokkos::deep_copy(M2, 0);

  /** accumulate the mean and variance
   * Numerically Stable Parallel Computation of (Co-)Variance, Schubert & Gertz
   */
  auto const ny_total = ny * n_procs;
#pragma omp parallel for
  for (int k = 0; k < mean_f.extent_int(1); ++k) {
    for (int p = 0; p < n_procs; ++p) {
#pragma omp simd
      for (int i = 0; i < mean_f.extent_int(0); ++i) {
        auto const delta_f = mean_f_gather(i, k, p) - mean_f(i, k);
        auto const delta_g = mean_g_gather(i, k, p) - mean_g(i, k);
        mean_f(i, k) += delta_f / (T)(p + 1);
        mean_g(i, k) += delta_g / (T)(p + 1);
        M2(i, k) +=
          M2_gather(i, k, p) + delta_f * delta_g * ny * p / (T)(p + 1);
      }
    }
#pragma omp simd
    for (int i = 0; i < mean_f.extent_int(0); ++i) {
      M2(i, k) /= ny_total;
    }
  }

  mpipp::barrier(comm);
}

template struct Covariance<Real, default_memory_pool>;

} // namespace alps
