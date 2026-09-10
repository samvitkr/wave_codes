//
// Created by xuanx004 on 10/4/22.
//

#include "solver.h"

#include <common/device/devices.h>
#include <common/kokkos_abstraction/pool_space.h>
#include <solvers/ns/convection_invoke.h>
#include <solvers/ns/curvilinear_common/convection_functor.h>
#include <solvers/ns/curvilinear_common/strain_rate.h>
#include <solvers/ns/curvilinear_common/viscous.h>
#include <solvers/operators/div.h>
#include <solvers/operators/grad_curvilinear.h>
#include <solvers/turbulence_model/models.h>
#include <solvers/turbulence_model/sgs_flux.h>

#include <Kokkos_Profiling_ScopedRegion.hpp>

namespace alps::solver {

Vector3Field<Real***> FlowOverWaveSolver::calc_explicit_rhs() const
{
  ++step;

  const auto& mesh     = flow_field.mesh;
  auto        tmp_mesh = static_cast<Mesh>(mesh);
  tmp_mesh.hbar        = 1;
  Vector3Field<Real***> fu;

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

  invoke_convection_fluxes_functor<CurveConvectionFunctor<BottomWaveMesh>>(
    fluxes, flow_field);

  /* Process turbulence model */
  if (std::holds_alternative<std::monostate>(turbulence_model)) {
    (void)0; // Do nothing
  } else if (std::holds_alternative<AnisotropicMinimumDissipation>(
               turbulence_model)) {
    auto grad_u = grad(flow_field.u.x, mesh, CenterPt());
    auto grad_v = grad(flow_field.u.y, mesh, CenterPt());
    auto grad_w = grad(flow_field.u.z, mesh, NodePt());

    Tensor33Field<Real***> const grad_velocities(grad_u, grad_v, grad_w);
    const auto*                  model =
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
      fluxes, flow_field, flow_field.nu_t, grad_u, grad_v, grad_w);
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

  /* Process viscous flux */
  if (options.integrator == "ab2cn") {
    add_viscous_fluxes_from_u<true>(fluxes, flow_field, 1 / Real(options.Re));
  } else if (options.integrator == "ab2") {
    add_viscous_fluxes_from_u<false>(fluxes, flow_field, 1 / Real(options.Re));
  } else {
    throw std::runtime_error(
      "Unsupported integrator type for viscous flux calculation");
  }

  set_boundary_stress_flux(fluxes);

  ALPS_CHECK_LAST_DEVICE_ERROR();

  fluxes = {};

  // Reuse the divergence for rectangular mesh to calculate dFₖ/dξₖ
  fu.z   = div(flux_w, tmp_mesh, NodePt(), "div(w_flux)");
  flux_w = {};
  fu.y   = div(flux_v, tmp_mesh, CenterPt(), "div(v_flux)");
  flux_v = {};
  fu.x   = div(flux_u, tmp_mesh, CenterPt(), "div(u_flux)");
  flux_u = {};

  for (const auto& body_force : body_forces) {
    body_force->add_forces(fu, Kokkos::DefaultExecutionSpace());
  }
  Kokkos::DefaultExecutionSpace().fence();

  ALPS_CHECK_LAST_DEVICE_ERROR();

  return fu;
}

} // namespace alps::solver
