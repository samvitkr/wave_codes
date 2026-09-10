//
// Created by xuanx004 on 5/26/23.
//

#include "log_law_wall_model.h"

#include <common/container/view_utils.h>
#include <common/kokkos_abstraction/exec_policy.h>
#include <common/kokkos_abstraction/pool_space.h>
#include <common/runtime/async_utils.h>
#include <solvers/field/flow_field.h>
#include <spectral/spectral.h>

#include <fmt/format.h>

#include <stdexcept>

namespace alps::solver {

LogLawWallModel::LogLawWallModel(LogLawWallModelOptions options)
  : options_{std::move(options)}
{}

void calculate_wall_shear_flux(MDView<Real**> const&           shear_flux_mag,
                               MDView<Real const** [2]> const& u_relative,
                               Real                            delta_z,
                               LogLawWallModel const&          wall_model)
{
  using Kokkos::parallel_for;

  // check if z0 is larger than the 1st grid
  if (wall_model.options_.z0 > (double)delta_z) {
    throw std::invalid_argument(
      fmt::format("Bottom wall model z0={} is larger than the off-wall "
                  "matching height {}",
                  wall_model.options_.z0,
                  delta_z));
  }

  const auto stream = get_next_stream();

  auto const          z0    = static_cast<Real>(wall_model.options_.z0);
  auto const          kappa = static_cast<Real>(wall_model.options_.kappa);
  const LoopPolicy<2> policy(
    stream, begins(shear_flux_mag), ends(shear_flux_mag));
  parallel_for(
    "wall_model_shear_stress", policy, KOKKOS_LAMBDA(int i, int j) {
      auto alpha  = Kokkos::log(delta_z / z0);
      auto u_mag  = Kokkos::hypot(u_relative(i, j, 0), u_relative(i, j, 1));
      auto u_star = u_mag * kappa / alpha;
      shear_flux_mag(i, j) = u_star * u_star;
    });

  stream.fence();
}

void apply_bottom_wall_shear_flux(Tensor33Field<Real***> const& fluxes,
                                  FlowField const&              flow,
                                  LogLawWallModel const&        wall_model)
{
  const auto& mesh      = flow.mesh;
  const auto& grid      = mesh.grid;
  const auto  is_bottom = grid.comm().is_first(2);

  if (!is_bottom) return;

  const auto stream1 = get_next_stream();
  const auto nx      = local_extent(fluxes.xx, 0);
  const auto ny      = local_extent(fluxes.xx, 1);

  // Variables for wall model
  MDView<Real**, default_memory_pool> wall_stress_flux_mag(
    "wall_stress_flux_mag", nx, ny);
  MDView<Real** [2], default_memory_pool> u_relative_wall("u relative", nx, ny);

  // calculate the velocity magnitudes
  auto const&         u = flow.u.x;
  auto const&         v = flow.u.y;
  const LoopPolicy<2> policy_2d(stream1, {0, 0}, {nx, ny});
  parallel_for(
    policy_2d, KOKKOS_LAMBDA(int i, int j) {
      u_relative_wall(i, j, 0) = u(i, j, 1) - u(i, j, 0);
      u_relative_wall(i, j, 1) = v(i, j, 1) - v(i, j, 0);
    });

  /* Filter the velocity supplied to the wall model, see the discussions in:
   *
   * Bou-Zeid, E., Meneveau, C., & Parlange, M. (2005). A scale-dependent
   * Lagrangian dynamic model for large eddy simulation of complex turbulent
   * flows. Physics of Fluids, 17(2), 025105.
   *
   * Yang, X. I. A., Park, G. I., & Moin, P. (2017). Log-layer mismatch and
   * modeling of the fluctuating wall stress in wall-modeled large-eddy
   * simulations. Physical Review Fluids, 2(10), 104601.
   */
  spectral::cutoff_xy(u_relative_wall,
                      grid.global_extent(0) / 4,
                      grid.global_extent(1) / 4,
                      grid,
                      policy_2d.space());
  policy_2d.space().fence();

  calculate_wall_shear_flux(wall_stress_flux_mag,
                            u_relative_wall,
                            mesh.zz_h(1) * mesh.hbar,
                            wall_model);

  // set the boundary flux to the wall stress given by the wall model
  const auto f_13 = create_inner_view(fluxes.xz);
  const auto f_23 = create_inner_view(fluxes.yz);
  parallel_for(
    "nu*S13 wall", policy_2d, KOKKOS_LAMBDA(int i, int j) {
      auto inv_u_mag =
        alps::rhypot(u_relative_wall(i, j, 0), u_relative_wall(i, j, 1));
      f_13(i, j, 0) =
        wall_stress_flux_mag(i, j) * u_relative_wall(i, j, 0) * inv_u_mag;
      f_23(i, j, 0) =
        wall_stress_flux_mag(i, j) * u_relative_wall(i, j, 1) * inv_u_mag;
    });

  policy_2d.space().fence();
}

std::string LogLawWallModel::info() const
{
  return fmt::format(
    "{} (kappa={}, z0={:e})", name, options_.kappa, options_.z0);
}

} // namespace alps::solver
