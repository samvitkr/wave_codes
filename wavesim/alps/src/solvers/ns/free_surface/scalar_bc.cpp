//
// Created by xuanx004 on 7/13/24.
//

#include "solver.h"

#include <common/base/logging.h>
#include <common/kokkos_abstraction/exec_policy.h>
#include <common/runtime/async_utils.h>
#include <decomp/mdcomm.h>
#include <solvers/field/scalar_bc_types.h>

#include <solvers/ns/curvilinear_common/scalar_bc.h>

namespace alps::solver {

void FreeSurfaceSolver::apply_scalar_bc(
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
      apply_top_bc(scalar.array, *bc, flow_field.mesh, stream1);
      bc_set = true;
    }
    if (auto const* bc =
          dynamic_cast<ConstantGradientBC const*>(scalar.top_bc.get());
        bc != nullptr) {
      apply_top_bc(scalar.array, *bc, flow_field.mesh, stream1);
      bc_set = true;
    }
    if (auto const* bc =
          dynamic_cast<ConstantFluxBC const*>(scalar.top_bc.get());
        bc != nullptr) {
      auto const scalarD = 1 / (options.Re * scalar_options.ScPr.value());
      apply_top_bc(scalar.array, *bc, (Real)scalarD, flow_field.mesh, stream1);
      bc_set = true;
    }

    if (!bc_set) {
      logger->warn("No top boundary condition set for scalar {}",
                   scalar.label());
    }
  }

  if (boundary == WhichBoundary::BottomBC || boundary == WhichBoundary::Both) {
    bool bc_set = false;
    if (auto const* bc =
          dynamic_cast<ConstantDirichletBC const*>(scalar.bottom_bc.get());
        bc != nullptr) {
      solver::apply_bottom_bc(scalar.array, *bc, flow_field.mesh, space);
      bc_set = true;
    }
    if (auto const* bc =
          dynamic_cast<ConstantGradientBC const*>(scalar.bottom_bc.get());
        bc != nullptr) {
      solver::apply_bottom_bc(scalar.array, *bc, flow_field.mesh, space);
      bc_set = true;
    }
    if (auto const* bc =
          dynamic_cast<ConstantFluxBC const*>(scalar.bottom_bc.get());
        bc != nullptr) {
      auto const scalarD = 1 / (options.Re * scalar_options.ScPr.value());
      solver::apply_bottom_bc(
        scalar.array, *bc, (Real)scalarD, flow_field.mesh, space);
      bc_set = true;
    }

    if (!bc_set) {
      logger->warn("No bottom boundary condition set for scalar {}",
                   scalar.label());
    }
  }

  space.fence();
  stream1.fence();
}

void FreeSurfaceSolver::set_boundary_scalar_flux(
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
        subview(fluxes.z, Kokkos::ALL, Kokkos::ALL, mesh.extent(2) - 2).view();
      auto const flux = bc->flux_;
      Kokkos::parallel_for(
        LoopPolicy<2>(stream, {0, 0}, {mesh.extent(0), mesh.extent(1)}),
        KOKKOS_LAMBDA(int i, int j) {
          // the specified flux is the flux across the boundary over unit dξd𝜓
          fz_top(i, j) = -flux;
        });
    }
  }

  if (is_bottom
      && (boundary == WhichBoundary::BottomBC
          || boundary == WhichBoundary::Both)) {
    if (auto const* bc =
          dynamic_cast<ConstantFluxBC const*>(scalar.bottom_bc.get());
        bc != nullptr) {
      auto const& fz_bottom =
        subview(fluxes.z, Kokkos::ALL, Kokkos::ALL, 0).view();
      auto const flux = bc->flux_;
      Kokkos::parallel_for(
        LoopPolicy<2>(stream, {0, 0}, {mesh.extent(0), mesh.extent(1)}),
        KOKKOS_LAMBDA(int i, int j) { fz_bottom(i, j) = -flux; });
    }
  }

  stream.fence();
}

} // namespace alps::solver
