//
// Created by xuanx004 on 7/13/24.
//

#include "solver.h"

#include <common/base/logging.h>
#include <common/device/device_traits.h>
#include <common/device/devices.h>
#include <common/kokkos_abstraction/exec_policy.h>
#include <common/kokkos_abstraction/pool_space.h>
#include <common/runtime/async_utils.h>
#include <decomp/ghost_cell_exchange.h>
#include <solvers/ns/advection_invoke.h>
#include <solvers/ns/curvilinear_common/advection_functor.h>
#include <solvers/operators/div.h>
#include <solvers/operators/grad_curvilinear.h>
#include <solvers/operators/lap_curvilinear.h>
#include <solvers/turbulence_model/sgs_scalar_flux.h>
#include <spectral/spectral.h>

#include <Kokkos_Core.hpp>
#include <Kokkos_Profiling_ScopedRegion.hpp>

namespace alps::solver {

namespace {
void step_one_scalar(HaloView<Real***> const&     f,
                     MDView<Real const***> const& Rf,
                     MDView<Real const***> const& Rf0,
                     MDView<Real const**> const&  invJ0,
                     Real                         dt)
{
  auto constexpr tile_size = []() -> Kokkos::Array<int, 3> {
    using Device = Kokkos::DefaultExecutionSpace;
    if constexpr (is_cuda_execution_space_v<Device>) return {16, 8, 1};
    if constexpr (is_hip_execution_space_v<Device>) return {64, 4, 1};
    return {0, 0, 0};
  }();
  auto const policy = LoopPolicy<3>(get_next_stream(),
                                    {0, 0, 0},
                                    {Rf.extent(0), Rf.extent(1), Rf.extent(2)},
                                    {tile_size[0], tile_size[1], tile_size[2]});
  if (Rf0.is_allocated()) {
    auto const alpha = -dt / 2;
    auto const beta  = (Real)1.5 * dt;
    Kokkos::parallel_for(
      "step " + f.label(), policy, KOKKOS_LAMBDA(int i, int j, int k) {
        f(i, j, k) = Kokkos::fma(
          beta,
          Rf(i, j, k),
          Kokkos::fma(alpha, Rf0(i, j, k), f(i, j, k) * invJ0(i, j)));
      });
  } else {
    Kokkos::parallel_for(
      "step " + f.label(), policy, KOKKOS_LAMBDA(int i, int j, int k) {
        f(i, j, k) = Kokkos::fma(Rf(i, j, k), dt, f(i, j, k) * invJ0(i, j));
      });
  }
  policy.space().fence();
}
} // namespace

void FlowOverWaveSolver::advance_scalars_ab2()
{
  auto const& mesh = flow_field.mesh;

  auto tmp_mesh = static_cast<Mesh>(mesh);
  tmp_mesh.hbar = 1;

  for (std::size_t i_scalar = 0; i_scalar < flow_field.scalars.size();
       ++i_scalar) {
    auto const region =
      Kokkos::Profiling::ScopedRegion("scalar " + std::to_string(i_scalar));

    auto const& scalar         = flow_field.scalars.at(i_scalar);
    auto const  scalar_options = options.scalars.at(i_scalar);
    auto const  scalarD        = 1 / (options.Re * scalar_options.ScPr.value());

    Vector3Field<Real***, default_memory_pool> c_flux(
      Kokkos::view_alloc(scalar.label() + " flux", Kokkos::WithoutInitializing),
      mesh.extents(),
      {0, 0, 1});

    invoke_advection_fluxes_functor<CurveAdvectionFunctor<BottomWaveMesh>>(
      c_flux, scalar.array, flow_field);

    auto const& sgs_model = scalar_sgs_models.at(i_scalar);
    if (std::holds_alternative<std::monostate>(sgs_model)) {
      add_laplacian_fluxes(c_flux,
                           scalar.array,
                           mesh,
                           (Real)scalarD,
                           "diffusion " + scalar.array.label());
    } else if (std::holds_alternative<ConstantTurbulentSc>(sgs_model)) {
      auto const* model = std::get_if<ConstantTurbulentSc>(&sgs_model);

      calculate_eddy_diffusivity(scalar.nuD, *model, flow_field);

      add_SGS_and_molecular_diffusion_fluxes(
        c_flux, scalar.array, flow_field, scalar.nuD, (Real)scalarD);
    } else if (std::holds_alternative<DynamicSmagorinskyScalar>(sgs_model)) {
      auto const* model = std::get_if<DynamicSmagorinskyScalar>(&sgs_model);

      auto grad_f = grad(scalar.array, mesh, CenterPt());

      calculate_eddy_diffusivity(scalar.nuD,
                                 *model,
                                 scalar,
                                 grad_f,
                                 flow_field,
                                 (step % options.C0_update_frequency != 0));
      add_sgs_and_molecular_diffusion_fluxes_from_gradf(
        c_flux, grad_f, flow_field, scalar.nuD, (Real)scalarD);
    } else if (std::holds_alternative<AnisotropicMinimumDissipationScalar>(
                 sgs_model)) {
      auto const* model =
        std::get_if<AnisotropicMinimumDissipationScalar>(&sgs_model);

      auto grad_f = grad(scalar.array, mesh, CenterPt());

      calculate_eddy_diffusivity(
        scalar.nuD, *model, scalar, grad_f, flow_field);
      add_sgs_and_molecular_diffusion_fluxes_from_gradf(
        c_flux, grad_f, flow_field, scalar.nuD, (Real)scalarD);
    }

    set_boundary_scalar_flux(c_flux, i_scalar);

    ALPS_CHECK_LAST_DEVICE_ERROR();

    auto fc = div(c_flux, tmp_mesh, CenterPt());

    auto const stream = get_next_stream();
    for (auto const& source : scalar_sources.at(i_scalar)) {
      source->add_source(fc, stream);
    }
    stream.fence();

    c_flux = {};

    MDView<Real***>& Rf0 = Rc_saved.at(i_scalar);
    step_one_scalar(scalar.array, fc, Rf0, mesh.invJ, (Real)dt);

    std::swap(Rf0, fc);
    fc = {};

    spectral::dealias(
      create_inner_view(scalar.array).view(), mesh.grid, stream);
    stream.fence();

    ALPS_CHECK_LAST_DEVICE_ERROR();
  }
}

void FlowOverWaveSolver::scale_scalars_by_J() const
{
  auto const& mesh = flow_field.mesh;
  auto const& J    = mesh.J;

  LoopPolicy<3> policy(Kokkos::DefaultExecutionSpace(),
                       {0, 0, 0},
                       {mesh.extent(0), mesh.extent(1), mesh.extent(2)});
  for (std::size_t i_scalar = 0; i_scalar < flow_field.scalars.size();
       ++i_scalar) {
    auto const region = Kokkos::Profiling::ScopedRegion(
      "scale scalar " + std::to_string(i_scalar));

    auto const& scalar = flow_field.scalars.at(i_scalar);
    auto const& array  = scalar.array;

    auto const stream = get_next_stream();
    apply_scalar_bc(i_scalar, stream);
    stream.fence();

    ALPS_CHECK_LAST_DEVICE_ERROR();

    if (mesh.comm().rank() == 0) {
      logger->info("Scalar {} updated", scalar.label());
    }
    policy.m_space = stream;

    bool all_valid = true;
    Kokkos::parallel_reduce(
      "scale " + scalar.label() + " by J",
      policy,
      KOKKOS_LAMBDA(int i, int j, int k, bool& valid) {
        array(i, j, k) *= J(i, j);
        valid = valid && Kokkos::isfinite(array(i, j, k));
      },
      Kokkos::LAnd<bool>(all_valid));
    policy.space().fence();

    if (!all_valid) {
      throw std::runtime_error("Scalar field " + scalar.label()
                               + " contains NaN or Inf");
    }

    update_halo_z(mesh.partition(), array, (uint8_t)i_scalar);
  }
}

} // namespace alps::solver
