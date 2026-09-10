//
// Created by xuanx004 on 7/11/24.
//

#include "sgs_flux_curve-inl.h"
#include "sgs_scalar_flux_curve-inl.h"
#include <solvers/field/flow_over_wave_field.h>
#include <solvers/ns/curvilinear_common/add_boussinesq_buoyancy.h>
#include <solvers/turbulence_model/anisotropic_minimum_dissipation.h>
#include <solvers/turbulence_model/constant_schmidt.h>
#include <solvers/turbulence_model/dynamic_smagorinsky.h>
#include <solvers/turbulence_model/sgs_flux.h>
#include <solvers/turbulence_model/smagorinsky.h>

namespace alps::solver {

void calculate_eddy_viscosity(
  HaloView<Real***> const&                   nu,
  AnisotropicMinimumDissipation const&       sgs_model,
  Tensor33Field<Real***> const&              grad_u,
  FlowOverWaveField const&                   flow,
  BoussinesqForce<FlowOverWaveSolver> const* buoyancy)
{
  auto const flat_flow = flow.as_channel_flow_field();

  if (buoyancy != nullptr) {
    // calculate buoyancy force and its gradient
    auto const stream        = get_next_stream();
    auto const [nx, ny, nz]  = flow.mesh.extents();
    auto const buoyancy_flux = HaloView<Real***, default_memory_pool>(
      "buoyancy_flux", {0, nx - 1}, {0, ny - 1}, {-1, nz});
    Kokkos::deep_copy(stream, buoyancy_flux.view(), 0);
    // calculate the buoyance term then take
    add_buoyancy_force_without_ref_scalar_curve(buoyancy_flux,
                                                buoyancy->get_scalar()->array,
                                                buoyancy->get_Ri(),
                                                flow.mesh,
                                                stream);
    // scale the buoyance term by J
    auto const& J = flow.mesh.J;
    Kokkos::parallel_for(
      "scale buoyancy",
      LoopPolicy<3>(stream, {0, 0, 0}, {nx, ny, nz}),
      KOKKOS_LAMBDA(int i, int j, int k) {
        buoyancy_flux(i, j, k) *= J(i, j);
      });
    stream.fence();
    update_halo_lower_z(flow.mesh.grid.partition(), buoyancy_flux, 1);
    auto const grad_b = grad(buoyancy_flux, flow.mesh, NodePt());
    detail::calculate_eddy_viscosity_impl<true>(
      nu, sgs_model, grad_u, flow.mesh, flow.time, &grad_b);
  } else {
    detail::calculate_eddy_viscosity_impl<false>(
      nu, sgs_model, grad_u, flow.mesh, flow.time, nullptr);
  }
}

void calculate_eddy_viscosity(HaloView<Real***> const&          nu,
                              DynamicSmagorinsky const&         sgs_model,
                              SymmTensor33Field<Real***> const& Sij,
                              FlowOverWaveField const&          flow,
                              bool                              skip_update_C0)
{
  auto const tmp_flow = flow.as_channel_flow_field();

  calculate_eddy_viscosity(nu, sgs_model, Sij, tmp_flow, skip_update_C0);
}

void calculate_eddy_viscosity(HaloView<Real***> const&          nu,
                              ConstantSmagorinsky const&        sgs_model,
                              SymmTensor33Field<Real***> const& Sij,
                              FlowOverWaveField const&          flow)
{
  auto const tmp_flow = flow.as_channel_flow_field();

  calculate_eddy_viscosity(nu, sgs_model, Sij, tmp_flow);
}

void add_sgs_momentum_fluxes_from_gradu(Tensor33Field<Real***> const& fluxes,
                                        FlowOverWaveField const&      flow,
                                        HaloView<Real***> const&      nu_t,
                                        Vector3Field<Real***> const&  grad_u,
                                        Vector3Field<Real***> const&  grad_v,
                                        Vector3Field<Real***> const&  grad_w)
{
  detail::add_sgs_momentum_fluxes_from_gradu(
    fluxes, flow, nu_t, grad_u, grad_v, grad_w);
}

void add_sgs_momentum_fluxes_from_Sij(Tensor33Field<Real***> const&     fluxes,
                                      FlowOverWaveField const&          flow,
                                      HaloView<Real***> const&          nu_t,
                                      SymmTensor33Field<Real***> const& Sij)
{
  detail::add_sgs_momentum_fluxes_from_Sij(fluxes, flow, nu_t, Sij);
}

//=============================================================================
// Scalars
//=============================================================================
void calculate_eddy_diffusivity(HaloView<Real***> const&   nuD,
                                ConstantTurbulentSc const& sgs_model,
                                FlowOverWaveField const&   flow)
{
  auto const flat_flow = flow.as_channel_flow_field();

  calculate_eddy_diffusivity(nuD, sgs_model, flat_flow);
}

void calculate_eddy_diffusivity(
  HaloView<Real***> const&        nuD,
  DynamicSmagorinskyScalar const& scalar_sgs_model,
  ScalarField const&              scalar,
  Vector3Field<Real***> const&    grad_f,
  FlowOverWaveField const&        flow,
  bool                            skip_update_C0)
{
  auto const flat_flow = flow.as_channel_flow_field();

  calculate_eddy_diffusivity(
    nuD, scalar_sgs_model, scalar, grad_f, flat_flow, skip_update_C0);
}

void calculate_eddy_diffusivity(
  HaloView<Real***> const&                   De,
  AnisotropicMinimumDissipationScalar const& scalar_sgs_model,
  ScalarField const&                         f,
  Vector3Field<Real***> const&               grad_f,
  FlowOverWaveField const&                   flow)
{
  auto const flat_flow = flow.as_channel_flow_field();

  calculate_eddy_diffusivity(De, scalar_sgs_model, f, grad_f, flat_flow);
}

void add_SGS_and_molecular_diffusion_fluxes(Vector3Field<Real***> const& fluxes,
                                            HaloView<Real const***> const& f,
                                            FlowOverWaveField const&       flow,
                                            HaloView<Real***> const&       nu_D,
                                            Real                           D,
                                            Real gamma)
{
  detail::add_SGS_and_molecular_diffusion_fluxes_impl(
    fluxes, f, flow, nu_D, D * gamma);
}

void add_sgs_and_molecular_diffusion_fluxes_from_gradf(
  Vector3Field<Real***> const& fluxes,
  Vector3Field<Real***> const& grad_f,
  FlowOverWaveField const&     flow,
  HaloView<Real***> const&     nuD,
  Real                         D,
  Real                         gamma)
{
  detail::add_sgs_diffusion_fluxes(fluxes, grad_f, flow, nuD, D * gamma);
}

} // namespace alps::solver
