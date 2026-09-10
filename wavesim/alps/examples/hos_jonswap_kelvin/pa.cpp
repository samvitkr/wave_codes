#include "pa.h"

#include <common/device/devices.h>
#include <common/kokkos_abstraction/exec_policy.h>
#include <common/math.h>
#include <solvers/hos/solver.h>
#include <spectral/spectral.h>

#include <Kokkos_Core.hpp>

namespace alps {

void radial_moving_pressure(MDView<Real**>                       pa,
                            double                               time,
                            double                               Pmax,
                            double                               R,
                            std::pair<double, double>            x0_initial,
                            double                               U,
                            Grid const&                          grid,
                            Kokkos::DefaultExecutionSpace const& space)
{
  // x0 = x0_initial + U * time
  auto const x0       = x0_initial.first + U * time;
  auto const y0       = x0_initial.second;
  auto const Lx       = 2 * Kokkos::numbers::pi_v<double> / grid.pex;
  auto const Ly       = 2 * Kokkos::numbers::pi_v<double> / grid.pey;
  auto const i_offset = grid.offset(0);
  auto const j_offset = grid.offset(1);

  Kokkos::parallel_for(
    "set_radial_moving_pressure",
    LoopPolicy<2>(space, {0, 0}, {pa.extent(0), pa.extent(1)}),
    KOKKOS_LAMBDA(int i, int j) {
      auto const x = Lx * (i + i_offset) / pa.extent_int(0);
      auto const y = Ly * (j + j_offset) / pa.extent_int(1);
      auto const r = Kokkos::hypot(x - x0, y - y0);
      if (r <= R) {
        auto s     = r / R;
        auto s2    = s * s;
        auto s3    = s * s * s;
        auto s4    = s2 * s2;
        auto s8    = s4 * s4;
        auto Gamma = 1 - 462 * s3 * s3 + 1980 * s4 * s3 - 3465 * s8
                   + 3080 * s8 * s - 1386 * s8 * s2 + 252 * s8 * s3;
        pa(i, j) = static_cast<Real>(Pmax * Gamma);
      } else {
        pa(i, j) = 0;
      }
    });

  ALPS_CHECK_LAST_DEVICE_ERROR();
  space.fence();
}
} // namespace alps
