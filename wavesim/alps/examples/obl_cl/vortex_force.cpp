//
// Created by xuanx004 on 12/27/23.
//

#include "vortex_force.h"

#include <common/container/view_utils.h>
#include <common/device/device_traits.h>
#include <common/kokkos_abstraction/exec_policy.h>
#include <common/kokkos_abstraction/pool_space.h>
#include <common/runtime/async_utils.h>
#include <decomp/pencil_plan.h>
#include <spectral/spectral.h>

#include <Kokkos_Core.hpp>
#include <fmt/format.h>

VortexForce::VortexForce(const alps::solver::ChannelFlowSolverAB2& flow_solver,
                         StokesDriftMonochromaticWave              Us_model)
  : solver{flow_solver}
  , stokes_drift{std::move(Us_model)}
{}

std::string VortexForce::info() const
{
  return fmt::format("Vortex force with Stokes drift of a monochromatic wave: "
                     "k_wave = {}, Us0 = {}",
                     stokes_drift.k_wave,
                     stokes_drift.Us0);
}

void VortexForce::add_forces(const alps::Vector3Field<alps::Real***>& fu,
                             const Kokkos::DefaultExecutionSpace& space) const
{
  using Kokkos::parallel_for;

  const auto& mesh = solver.flow_field.mesh;
  const auto& grid = mesh.grid;
  const auto& hbar = mesh.hbar;

  const auto u         = solver.flow_field.u.x;
  const auto v         = solver.flow_field.u.y;
  const auto w         = solver.flow_field.u.z;
  const auto nx        = alps::local_extent(u, 0);
  const auto ny        = alps::local_extent(u, 1);
  const auto nz        = alps::local_extent(u, 2);
  const auto is_top    = grid.comm().is_last(2);
  const auto is_bottom = grid.comm().is_first(2);

  space.fence();
  auto stream = alps::get_next_stream();
  auto constexpr tile =
    []() -> Kokkos::Array<alps::LoopPolicy<3>::array_index_type, 3> {
    using Device = Kokkos::DefaultExecutionSpace;
    if constexpr (alps::is_cuda_execution_space_v<Device>) return {16, 1, 8};
    if constexpr (alps::is_hip_execution_space_v<Device>) return {64, 1, 4};
    return {0, 0, 0};
  }();

  using tmp_array_t = alps::MDView<alps::Real***, alps::default_memory_pool>;
  tmp_array_t const tmp1("tmp1", alps::create_local_layout(grid.pencil));
  tmp_array_t const tmp2("tmp2", alps::create_local_layout(grid.pencil));

  {
    alps::spectral::ddy(tmp1, alps::create_inner_view(u).view(), grid, stream);
    alps::spectral::ddx(tmp2, alps::create_inner_view(v).view(), grid, stream);

    const auto& zz = mesh.zz;
    const auto& f  = fu.y;
    const auto& Us = stokes_drift;
    // vortex force is not applied at the boundaries
    alps::LoopPolicy<3> const policy(
      stream, {0, 0, is_bottom ? 1 : 0}, {nx, ny, is_top ? nz - 1 : nz}, tile);
    parallel_for(
      "vortex force fy", policy, KOKKOS_LAMBDA(int i, int j, int k) {
        // fy = -Us * omg_z = Us * (uy - vx)
        f(i, j, k) +=
          Us.Us((zz(k) - 1) * hbar) * (tmp1(i, j, k) - tmp2(i, j, k));
      });
  }

  {
    alps::spectral::ddx(tmp1, alps::create_inner_view(w).view(), grid, stream);

    const auto& dz = mesh.dz;
    const auto& zw = mesh.zw;
    const auto& f  = fu.z;
    const auto& Us = stokes_drift;
    // vortex force is not applied at the boundaries
    alps::LoopPolicy<3> const policy(
      stream, {0, 0, is_bottom ? 1 : 0}, {nx, ny, is_top ? nz - 2 : nz}, tile);
    parallel_for(
      "vortex force fz", policy, KOKKOS_LAMBDA(int i, int j, int k) {
        // fz = Us * omg_y = Us * (uz - wx)
        auto alpha = 1 / (dz(k) * hbar);
        auto uz    = (u(i, j, k + 1) - u(i, j, k)) * alpha;
        f(i, j, k) += Us.Us((zw(k) - 1) * hbar) * (uz - tmp1(i, j, k));
      });
  }

  stream.fence();
}

StokesDriftMonochromaticWave
create_stokes_drift_from_config(alps::ConfigTable const& config)
{
  auto const Us0 = config.get_value<double>("Us0");
  auto const kw  = config.get_value<double>("k");

  return {kw, Us0};
}
