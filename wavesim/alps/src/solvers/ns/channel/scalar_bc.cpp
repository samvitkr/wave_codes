//
// Created by xuanx004 on 7/12/24.
//

#include "solver_ab2.h"

#include <common/base/logging.h>
#include <common/container/view_utils.h>
#include <common/kokkos_abstraction/exec_policy.h>
#include <common/runtime/async_utils.h>
#include <decomp/mdcomm.h>

namespace alps::solver {

void apply_top_bc(HaloView<Real***> const&             f,
                  ConstantGradientBC const&            bc,
                  FlowField const&                     flow,
                  Kokkos::DefaultExecutionSpace const& stream)
{
  auto const& mesh = flow.mesh;
  if (!flow.mesh.comm().is_last(2)) return;

  auto const  nz   = mesh.extent(2);
  auto const  grad = bc.grad_;
  auto const  hbar = mesh.hbar;
  auto const& dz   = mesh.dz;
  parallel_for(
    LoopPolicy<2>(stream, {0, 0}, {mesh.extent(0), mesh.extent(1)}),
    KOKKOS_LAMBDA(int i, int j) {
      auto alpha  = dz(nz - 2) + dz(nz - 3);
      auto beta   = alpha + dz(nz - 2);
      auto coeff1 = alpha * alpha / dz(nz - 3) / beta;
      auto coeff2 = -dz(nz - 2) * dz(nz - 2) / dz(nz - 3) / beta;
      auto coeff0 = dz(nz - 2) * alpha / beta * hbar;
      f(i, j, nz - 1) =
        coeff1 * f(i, j, nz - 2) + coeff2 * f(i, j, nz - 3) + coeff0 * grad;
    });
}

void apply_bottom_bc(HaloView<Real***> const&             f,
                     ConstantGradientBC const&            bc,
                     FlowField const&                     flow,
                     Kokkos::DefaultExecutionSpace const& stream)
{
  auto const& mesh = flow.mesh;
  if (!flow.mesh.comm().is_first(2)) return;

  auto const  grad = bc.grad_;
  auto const  hbar = mesh.hbar;
  auto const& dz   = mesh.dz;
  parallel_for(
    LoopPolicy<2>(stream, {0, 0}, {mesh.extent(0), mesh.extent(1)}),
    KOKKOS_LAMBDA(int i, int j) {
      auto alpha  = dz(0) + dz(1);
      auto beta   = alpha + dz(0);
      auto coeff1 = alpha * alpha / dz(1) / beta;
      auto coeff2 = -dz(0) * dz(0) / dz(1) / beta;
      auto coeff0 = -dz(0) * alpha / beta * hbar;
      f(i, j, 0)  = coeff1 * f(i, j, 1) + coeff2 * f(i, j, 2) + coeff0 * grad;
    });
}

void apply_top_bc(HaloView<Real***> const&             f,
                  ConstantDirichletBC const&           bc,
                  FlowField const&                     flow,
                  Kokkos::DefaultExecutionSpace const& stream)
{
  auto const& mesh = flow.mesh;
  if (!flow.mesh.comm().is_last(2)) return;

  auto const nz  = mesh.extent(2);
  auto const val = bc.value_;
  parallel_for(
    LoopPolicy<2>(stream, {0, 0}, {mesh.extent(0), mesh.extent(1)}),
    KOKKOS_LAMBDA(int i, int j) { f(i, j, nz - 1) = val; });
}

void apply_bottom_bc(HaloView<Real***> const&             f,
                     ConstantDirichletBC const&           bc,
                     FlowField const&                     flow,
                     Kokkos::DefaultExecutionSpace const& stream)
{
  auto const& mesh = flow.mesh;
  if (!flow.mesh.comm().is_first(2)) return;

  auto const val = bc.value_;
  parallel_for(
    LoopPolicy<2>(stream, {0, 0}, {mesh.extent(0), mesh.extent(1)}),
    KOKKOS_LAMBDA(int i, int j) { f(i, j, 0) = val; });
}

void apply_top_bc(HaloView<Real***> const&             f,
                  ConstantFluxBC const&                bc,
                  Real                                 D,
                  FlowField const&                     flow,
                  Kokkos::DefaultExecutionSpace const& stream)
{
  if (!flow.mesh.comm().is_last(2)) return;

  ConstantGradientBC grad_bc;
  grad_bc.grad_ = -bc.flux_ / D;
  apply_top_bc(f, grad_bc, flow, stream);
}

void apply_bottom_bc(HaloView<Real***> const&             f,
                     ConstantFluxBC const&                bc,
                     Real                                 D,
                     FlowField const&                     flow,
                     Kokkos::DefaultExecutionSpace const& stream)
{
  if (!flow.mesh.comm().is_first(2)) return;

  ConstantGradientBC grad_bc;
  grad_bc.grad_ = -bc.flux_ / D;
  apply_bottom_bc(f, grad_bc, flow, stream);
}

void ChannelFlowSolverAB2::apply_scalar_bc(
  std::size_t                          scalar_index,
  const Kokkos::DefaultExecutionSpace& space,
  WhichBoundary                        boundary) const
{
  auto const& stream1        = get_next_stream();
  auto const& scalars        = flow_field.scalars.storage();
  auto const& scalar         = scalars.at(scalar_index);
  auto const& scalar_options = options.scalars.at(scalar_index);

  if (boundary == WhichBoundary::TopBC || boundary == WhichBoundary::Both) {
    auto event = get_device_event();
    enqueue(event, space);
    wait_for(event, stream1);

    bool bc_set = false;
    if (auto const* bc =
          dynamic_cast<ConstantDirichletBC const*>(scalar.top_bc.get());
        bc != nullptr) {
      apply_top_bc(scalar.array, *bc, flow_field, stream1);
      bc_set = true;
    }
    if (auto const* bc =
          dynamic_cast<ConstantGradientBC const*>(scalar.top_bc.get());
        bc != nullptr) {
      apply_top_bc(scalar.array, *bc, flow_field, stream1);
      bc_set = true;
    }
    if (auto const* bc =
          dynamic_cast<ConstantFluxBC const*>(scalar.top_bc.get());
        bc != nullptr) {
      auto const scalarD = 1 / (options.Re * scalar_options.ScPr.value());
      apply_top_bc(scalar.array, *bc, (Real)scalarD, flow_field, stream1);
      bc_set = true;
    }

    if (!bc_set) {
      logger->warn("No top boundary condition set for scalar "
                   + scalar.label());
    }
  }

  if (boundary == WhichBoundary::BottomBC || boundary == WhichBoundary::Both) {
    bool bc_set = false;
    if (auto const* bc =
          dynamic_cast<ConstantDirichletBC const*>(scalar.bottom_bc.get());
        bc != nullptr) {
      apply_bottom_bc(scalar.array, *bc, flow_field, space);
      bc_set = true;
    }
    if (auto const* bc =
          dynamic_cast<ConstantGradientBC const*>(scalar.bottom_bc.get());
        bc != nullptr) {
      apply_bottom_bc(scalar.array, *bc, flow_field, space);
      bc_set = true;
    }
    if (auto const* bc =
          dynamic_cast<ConstantFluxBC const*>(scalar.bottom_bc.get());
        bc != nullptr) {
      auto const scalarD = 1 / (options.Re * scalar_options.ScPr.value());
      apply_bottom_bc(scalar.array, *bc, (Real)scalarD, flow_field, space);
      bc_set = true;
    }

    if (!bc_set) {
      logger->warn("No bottom boundary condition set for scalar "
                   + scalar.label());
    }
  }

  space.fence();
  stream1.fence();
}

void ChannelFlowSolverAB2::set_boundary_scalar_flux(
  Vector3Field<Real***> const& fluxes,
  std::size_t                  scalar_index,
  WhichBoundary                boundary) const
{
  auto const& stream  = get_next_stream();
  auto const& scalars = flow_field.scalars.storage();
  auto const& scalar  = scalars.at(scalar_index);

  auto const& mesh = flow_field.mesh;

  auto const is_top    = mesh.comm().is_last(2);
  auto const is_bottom = mesh.comm().is_first(2);

  if (is_top
      && (boundary == WhichBoundary::TopBC
          || boundary == WhichBoundary::Both)) {
    if (auto const* bc =
          dynamic_cast<ConstantFluxBC const*>(scalar.top_bc.get());
        bc != nullptr) {
      auto const& fz_top =
        subview(fluxes.z, ALL, ALL, mesh.extent(2) - 2).view();
      auto const flux = bc->flux_;
      Kokkos::parallel_for(
        LoopPolicy<2>(stream, {0, 0}, {mesh.extent(0), mesh.extent(1)}),
        KOKKOS_LAMBDA(int i, int j) { fz_top(i, j) = -flux; });
    }
  }

  if (is_bottom
      && (boundary == WhichBoundary::BottomBC
          || boundary == WhichBoundary::Both)) {
    if (auto const* bc =
          dynamic_cast<ConstantFluxBC const*>(scalar.bottom_bc.get());
        bc != nullptr) {
      auto const& fz_bottom = subview(fluxes.z, ALL, ALL, 0).view();
      auto const  flux      = bc->flux_;
      Kokkos::parallel_for(
        LoopPolicy<2>(stream, {0, 0}, {mesh.extent(0), mesh.extent(1)}),
        KOKKOS_LAMBDA(int i, int j) { fz_bottom(i, j) = -flux; });
    }
  }

  stream.fence();
}

} // namespace alps::solver
