#include "series.h"

#include <common/base/macros.h>
#include <common/container/view_utils.h>
#include <common/kokkos_abstraction/exec_policy.h>
#include <common/kokkos_abstraction/pool_space.h>
#include <spectral/spectral.h>

#include <Kokkos_Core.hpp>

namespace alps::solver::hos {

void taylor_series_coeff_async(MDView<Real***>                      zp_hos,
                               MDView<const Real**>                 eta,
                               int                                  npw,
                               Grid const&                          grid,
                               Kokkos::DefaultExecutionSpace const& space)
{
  if (zp_hos.extent_int(2) < npw - 1) {
    throw std::runtime_error("Taylor series coefficients array too small");
  }

  auto ends = local_extents(eta);

  Kokkos::deep_copy(space, subview(zp_hos, ALL, ALL, 0), eta);

  GridPolicy<> policy(space, ends[1], Kokkos::AUTO);
  for (int k = 1; k < npw - 1; ++k) {
    Kokkos::parallel_for(
      "taylor_coeff_" + std::to_string(k),
      policy,
      KOKKOS_LAMBDA(decltype(policy)::member_type team) {
        const int j = team.league_rank();
        Kokkos::parallel_for(
          Kokkos::TeamVectorRange(team, ends[0]), KOKKOS_TR_LAMBDA(int i) {
            zp_hos(i, j, k) = zp_hos(i, j, k - 1) * eta(i, j) / Real(k + 1);
          });
      });
    spectral::dealias(
      subview(zp_hos, ALL, ALL, index_range(k, k + 1)), grid, space);
  }
}

void surface_vp_expansion(MDView<Real***>                      r_hat,
                          MDView<const Real**>                 vps,
                          MDView<const Real***>                zp_hos,
                          MDView<const Real***>                wvn_hos,
                          Grid const&                          grid,
                          Kokkos::DefaultExecutionSpace const& space)
{
  const int npw = r_hat.extent_int(2);

  using Kokkos::view_alloc;
  MDView<Real***, default_memory_pool> dr_hat(
    view_alloc("r_hat", Kokkos::WithoutInitializing),
    wvn_hos.extent(0),
    wvn_hos.extent(1),
    zp_hos.extent(2)); // zp_hos.extent(2) == npw - 1
  MDView<Real***, default_memory_pool> dr(
    view_alloc("dr", Kokkos::WithoutInitializing), zp_hos.layout());

  spectral::fft_r2c_xy(
    subview(r_hat, ALL, ALL, index_range(0, 1)),
    MDView<const Real** [1]>(vps.data(), vps.extent(0), vps.extent(1)),
    grid,
    space);

  using Kokkos::parallel_for;
  using Kokkos::TeamVectorRange;
  auto tile_size = []() -> Kokkos::Array<size_t, 2> {
    using Device = Kokkos::DefaultExecutionSpace;
    if constexpr (is_cuda_execution_space_v<Device>) return {32, 8};
    if constexpr (is_hip_execution_space_v<Device>) return {64, 4};
    return {};
  }();
  LoopPolicy<2> policy_y(
    space, {0, 0}, {r_hat.extent(0), r_hat.extent(1)}, tile_size);
  LoopPolicy<2> policy_x(
    space, {0, 0}, {dr.extent(0), dr.extent(1)}, tile_size);
  for (int k = 1; k < npw; ++k) {
    parallel_for(
      "r_hat*wvn", policy_y, KOKKOS_LAMBDA(int i, int j) {
        for (int k1 = 0; k1 < k; ++k1) {
          int k2           = k - 1 - k1;
          dr_hat(i, j, k1) = r_hat(i, j, k2) * wvn_hos(i, j, k1);
        }
      });

    spectral::fft_c2r_xy(subview(dr, ALL, ALL, index_range(0, k)),
                         subview(dr_hat, ALL, ALL, index_range(0, k)),
                         false,
                         grid,
                         space);

    parallel_for(
      "rk", policy_x, KOKKOS_LAMBDA(int i, int j) {
        dr(i, j, 0) = -zp_hos(i, j, 0) * dr(i, j, 0);
        for (int k1 = 1; k1 < k; ++k1) {
          dr(i, j, 0) =
            Kokkos::fma(-zp_hos(i, j, k1), dr(i, j, k1), dr(i, j, 0));
        }
      });

    spectral::fft_r2c_xy(subview(r_hat, ALL, ALL, index_range(k, k + 1)),
                         subview(dr, ALL, ALL, index_range(0, 1)),
                         grid,
                         space);
  }

  space.fence();
}

void surface_w(MDView<Real**>                       ws,
               MDView<Real***>                      r_hat,
               MDView<const Real***>                zp_hos,
               MDView<const Real***>                wvn_hos,
               Grid const&                          grid,
               Kokkos::DefaultExecutionSpace const& space)
{
  using member_t = GridPolicy<>::member_type;

  MDView<Real***, default_memory_pool> pk(
    Kokkos::view_alloc("pk", Kokkos::WithoutInitializing),
    ws.extent(0),
    ws.extent(1),
    r_hat.extent(2));

  using Kokkos::parallel_for;
  GridPolicy<> policy_y(space, wvn_hos.extent_int(1), Kokkos::AUTO);
  parallel_for(
    "pk*wvn", policy_y, KOKKOS_LAMBDA(member_t const& team) {
      int       j    = team.league_rank();
      const int n_ky = r_hat.extent_int(0);
      const int npw  = r_hat.extent_int(2);
      for (int k = 1; k < npw; ++k) {
        parallel_for(
          Kokkos::TeamVectorRange(team, n_ky), KOKKOS_TR_LAMBDA(int i) {
            auto tmp = r_hat(i, j, k - 1);
            r_hat(i, j, k) += tmp;
            r_hat(i, j, k - 1) *= wvn_hos(i, j, npw - k);
          });
      }
      team.team_barrier();
      parallel_for(
        Kokkos::TeamVectorRange(team, n_ky),
        KOKKOS_TR_LAMBDA(int i) { r_hat(i, j, npw - 1) *= wvn_hos(i, j, 0); });
    });

  spectral::fft_c2r_xy(pk, r_hat, false, grid, space);

  GridPolicy<> policy_x(space, pk.extent_int(1), Kokkos::AUTO);
  parallel_for(
    "ws add", policy_x, KOKKOS_LAMBDA(member_t const& team) {
      int       j   = team.league_rank();
      const int nx  = ws.extent_int(0);
      const int npw = pk.extent_int(2);
      parallel_for(
        Kokkos::TeamVectorRange(team, nx),
        KOKKOS_TR_LAMBDA(int i) { ws(i, j) = pk(i, j, npw - 1); });
      for (int k = 0; k < npw - 1; ++k) {
        parallel_for(
          Kokkos::TeamVectorRange(team, nx), KOKKOS_TR_LAMBDA(int i) {
            ws(i, j) =
              Kokkos::fma(zp_hos(i, j, k), pk(i, j, npw - 2 - k), ws(i, j));
          });
      }
    });

  spectral::dealias(
    MDView<Real** [1]>(ws.data(), ws.extent(0), ws.extent(1)), grid, space);

  space.fence();
}

void calc_wavenumbers(MDView<Real***> wvn_hos,
                      Real            kx0,
                      Real            ky0,
                      Grid const&     grid)
{
  const auto offset_x = grid.offset(1, Pencil::Y);

  GridPolicy<> policy(wvn_hos.extent_int(1), Kokkos::AUTO);
  Kokkos::parallel_for(
    "calc wavenumbers",
    policy,
    KOKKOS_LAMBDA(decltype(policy)::member_type team) {
      const int j  = team.league_rank();
      double    ax = (double)kx0 * int((j + offset_x) / 2);
      for (int n = 0; n < wvn_hos.extent_int(2); ++n) {
        Kokkos::parallel_for(
          Kokkos::TeamVectorRange(team, wvn_hos.extent_int(0)),
          KOKKOS_TR_LAMBDA(int i) {
            double ay = (double)ky0 * int(i / 2);
            double a  = Kokkos::hypot(ax, ay);
            // calculate a^(n+1)
            double an = a;
            for (int m = 0; m < n; ++m) {
              an *= a;
            }
            wvn_hos(i, j, n) = static_cast<Real>(an);
          });
      }
    });

  policy.space().fence();
}

} // namespace alps::solver::hos
