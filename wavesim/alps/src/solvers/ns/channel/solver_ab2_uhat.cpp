//
// Created by xuanx004 on 6/30/24.
//

#include "solver_ab2.h"

#include "convection_diffusion.h"
#include "strain_rate.h"
#include <common/base/macros.h>
#include <common/container/matrix_field.h>
#include <common/container/vector_field.h>
#include <common/container/view_utils.h>
#include <common/device/devices.h>
#include <common/kokkos_abstraction/exec_policy.h>
#include <common/kokkos_abstraction/pool_space.h>
#include <common/math.h>
#include <common/runtime/async_utils.h>
#include <decomp/mdcomm.h>
#include <solvers/field/flow_field.h>
#include <solvers/operators/div.h>
#include <solvers/operators/grad.h>
#include <solvers/turbulence_model/log_law_wall_model.h>
#include <solvers/turbulence_model/models.h>
#include <solvers/turbulence_model/sgs_flux.h>

#include <Kokkos_Core.hpp>
#include <Kokkos_Profiling_ScopedRegion.hpp>

namespace alps::solver {

void ChannelFlowSolverAB2::calc_uhat(Vector3Field<Real***> Ru)
{
  // First update the scalar fields
  advance_scalars_ab2();

  // cleanup SGS model if needed
  if (auto* model = std::get_if<DynamicSmagorinsky>(&turbulence_model);
      model != nullptr) {
    model->cleanup_Smag();
  }
  if (auto* model =
        std::get_if<AnisotropicMinimumDissipation>(&turbulence_model);
      model != nullptr) {
    model->cleanup_grad();
  }

  using Kokkos::parallel_for;
  using Kokkos::TeamVectorRange;

  const auto& u = flow_field.u.x;
  const auto& v = flow_field.u.y;
  const auto& w = flow_field.u.z;

  auto ends = local_extents(u);

  Kokkos::Profiling::pushRegion("U stepping");
  // Advance the velocity field
  auto policy =
    GridPolicy<>(get_next_stream(), ends[1] * ends[2], Kokkos::AUTO);
  using member_t = GridPolicy<>::member_type;
  if (steps_initialized) {
    // Second-order Adam Bashforth
    const auto  alpha = Real(-dt / 2);
    const auto  beta  = Real(1.5 * dt);
    const auto& Ru_s  = this->Ru_saved;
    parallel_for(
      "u_stepping", policy, KOKKOS_LAMBDA(member_t team) {
        int k = team.league_rank() / ends[1];
        int j = team.league_rank() % ends[1];
        parallel_for(
          TeamVectorRange(team, ends[0]), KOKKOS_TR_LAMBDA(int& i) {
            u(i, j, k) =
              Kokkos::fma(beta,
                          Ru.x(i, j, k),
                          Kokkos::fma(alpha, Ru_s.x(i, j, k), u(i, j, k)));
          });
        parallel_for(
          TeamVectorRange(team, ends[0]), KOKKOS_TR_LAMBDA(int& i) {
            v(i, j, k) =
              Kokkos::fma(beta,
                          Ru.y(i, j, k),
                          Kokkos::fma(alpha, Ru_s.y(i, j, k), v(i, j, k)));
          });
        parallel_for(
          TeamVectorRange(team, ends[0]), KOKKOS_TR_LAMBDA(int& i) {
            w(i, j, k) =
              Kokkos::fma(beta,
                          Ru.z(i, j, k),
                          Kokkos::fma(alpha, Ru_s.z(i, j, k), w(i, j, k)));
          });
      });
  } else {
    // Forward Euler for the first step
    auto delta_t = (Real)dt;
    parallel_for(
      "u_stepping_0", policy, KOKKOS_LAMBDA(member_t team) {
        int k = team.league_rank() / ends[1];
        int j = team.league_rank() % ends[1];
        parallel_for(
          TeamVectorRange(team, ends[0]), KOKKOS_TR_LAMBDA(int& i) {
            u(i, j, k) = Kokkos::fma(Ru.x(i, j, k), delta_t, u(i, j, k));
            v(i, j, k) = Kokkos::fma(Ru.y(i, j, k), delta_t, v(i, j, k));
            w(i, j, k) = Kokkos::fma(Ru.z(i, j, k), delta_t, w(i, j, k));
          });
      });
  }
  policy.space().fence();

  ALPS_CHECK_LAST_DEVICE_ERROR();

  // Save the right-hand side
  std::swap(Ru_saved, Ru);

  // set boundary condition
  apply_bc(policy.space());

  ALPS_CHECK_LAST_DEVICE_ERROR();

  if (!steps_initialized) {
    steps_initialized = !steps_initialized;
  }

  Kokkos::Profiling::popRegion();
}

Vector3Field<Real***> ChannelFlowSolverAB2::calc_explicit_rhs() const
{
  ++step;

  const Mesh& mesh = flow_field.mesh;

  auto const rhs_region = Kokkos::Profiling::ScopedRegion("U rhs");

  Vector3Field<Real***, default_memory_pool> flux_u(
    Kokkos::view_alloc("U_flux", Kokkos::WithoutInitializing),
    mesh.extents(),
    {0, 0, 1});
  Vector3Field<Real***, default_memory_pool> flux_v(
    Kokkos::view_alloc("V_flux", Kokkos::WithoutInitializing),
    mesh.extents(),
    {0, 0, 1});
  Vector3Field<Real***, default_memory_pool> flux_w(
    Kokkos::view_alloc("W_flux", Kokkos::WithoutInitializing),
    mesh.extents(),
    {0, 0, 1});
  Tensor33Field<Real***> fluxes(flux_u, flux_v, flux_w);

  calculate_convection_fluxes(fluxes, flow_field);

  if (std::holds_alternative<std::monostate>(turbulence_model)) {
    (void)0; // Do nothing
  } else if (std::holds_alternative<AnisotropicMinimumDissipation>(
               turbulence_model)) {
    auto grad_u = grad(flow_field.u.x, mesh, CenterPt());
    auto grad_v = grad(flow_field.u.y, mesh, CenterPt());
    auto grad_w = grad(flow_field.u.z, mesh, NodePt());
    Tensor33Field<Real***> const grad_velocities(grad_u, grad_v, grad_w);

    const auto* model =
      std::get_if<AnisotropicMinimumDissipation>(&turbulence_model);

    // Look for the buoyancy term and pass it to the model
    using BuoyancyType = BoussinesqForce<SolverType>;
    BuoyancyType const* buoyancy{nullptr};
    for (const auto& body_force : body_forces) {
      if (auto const* term =
            dynamic_cast<BuoyancyType const*>(body_force.get());
          term != nullptr) {
        buoyancy = term;
        break;
      }
    }
    calculate_eddy_viscosity(
      flow_field.nu_t, *model, grad_velocities, flow_field, buoyancy);

    add_sgs_momentum_fluxes_from_gradu(
      fluxes, flow_field, flow_field.nu_t, 0, grad_u, grad_v, grad_w);
  } else if (std::holds_alternative<ConstantSmagorinsky>(turbulence_model)) {
    SymmTensor33Field<Real***, default_memory_pool> const Sij(
      Kokkos::view_alloc("Sij", Kokkos::WithoutInitializing),
      mesh.extents(),
      {0, 0, 1});
    calculate_strain_rate(Sij, flow_field.u, mesh);

    const auto* model = std::get_if<ConstantSmagorinsky>(&turbulence_model);
    calculate_eddy_viscosity(flow_field.nu_t, *model, Sij, flow_field);
    add_sgs_momentum_fluxes_from_Sij(fluxes, flow_field, flow_field.nu_t, Sij);
  } else if (std::holds_alternative<DynamicSmagorinsky>(turbulence_model)) {
    const auto* model = std::get_if<DynamicSmagorinsky>(&turbulence_model);
    SymmTensor33Field<Real***, default_memory_pool> const Sij(
      Kokkos::view_alloc("Sij", Kokkos::WithoutInitializing),
      mesh.extents(),
      {0, 0, 1});
    calculate_strain_rate(Sij, flow_field.u, mesh);

    calculate_eddy_viscosity(flow_field.nu_t,
                             *model,
                             Sij,
                             flow_field,
                             (step % options.C0_update_frequency != 0));
    add_sgs_momentum_fluxes_from_Sij(fluxes, flow_field, flow_field.nu_t, Sij);
  }

  // add molecular viscous fluxes
  add_viscous_fluxes_from_u(fluxes, flow_field, 1 / Real(options.Re));

  // Adjust the boundary stress according to the boundary condition
  set_boundary_stress_flux(fluxes);

  ALPS_CHECK_LAST_DEVICE_ERROR();

  fluxes = {};

  Vector3Field<Real***> fu;
  fu.z   = div(flux_w, flow_field.mesh, NodePt(), "div(w_flux)");
  flux_w = {};
  fu.y   = div(flux_v, flow_field.mesh, CenterPt(), "div(v_flux)");
  flux_v = {};
  fu.x   = div(flux_u, flow_field.mesh, CenterPt(), "div(u_flux)");
  flux_u = {};

  for (const auto& body_force : body_forces) {
    body_force->add_forces(fu, Kokkos::DefaultExecutionSpace());
  }
  Kokkos::DefaultExecutionSpace().fence();

  ALPS_CHECK_LAST_DEVICE_ERROR();

  return fu;
}

void ChannelFlowSolverAB2::set_boundary_stress_flux(
  Tensor33Field<Real***> const& fluxes,
  WhichBoundary                 boundary) const
{
  auto const& stream1 = get_next_stream();

  auto const& flow = this->flow_field;
  auto const& mesh = flow.mesh;

  auto const is_top    = mesh.comm().is_last(2);
  auto const is_bottom = mesh.comm().is_first(2);

  if (is_top
      && (boundary == WhichBoundary::TopBC
          || boundary == WhichBoundary::Both)) {
    if (auto const* bc =
          dynamic_cast<TangentialStressWall const*>(flow.top_bc.get());
        bc != nullptr) {
      auto const  ends       = local_ends(fluxes.xz);
      auto const& tau_13_top = subview(fluxes.xz, ALL, ALL, ends[2] - 2).view();
      auto const& tau_23_top = subview(fluxes.yz, ALL, ALL, ends[2] - 2).view();
      auto const  tau_1      = bc->tau_1;
      auto const  tau_2      = bc->tau_2;

      Kokkos::parallel_for(
        LoopPolicy<2>(stream1, {0, 0}, {ends[0], ends[1]}),
        KOKKOS_LAMBDA(int i, int j) {
          tau_13_top(i, j) = tau_1;
          tau_23_top(i, j) = tau_2;
        });
    }
  }

  // apply the stress from the wall model
  if (is_bottom
      && (boundary == WhichBoundary::BottomBC
          || boundary == WhichBoundary::Both)) {
    auto const* wall_model =
      dynamic_cast<LogLawWallModel const*>(wall_model_bottom.get());

    if (wall_model != nullptr) {
      apply_bottom_wall_shear_flux(fluxes, flow, *wall_model);
    }
  }

  stream1.fence();
}
} // namespace alps::solver
