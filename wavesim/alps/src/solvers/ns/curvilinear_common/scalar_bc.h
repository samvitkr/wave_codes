//
// Created by xuanx004 on 7/13/24.
//

#pragma once

#include <common/container/view_types.h>
#include <solvers/field/scalar_bc_types.h>

namespace alps::solver {

class CurvilinearMesh;

void apply_top_bc(HaloView<Real***> const&             f,
                  ConstantDirichletBC const&           bc,
                  CurvilinearMesh const&               mesh,
                  Kokkos::DefaultExecutionSpace const& stream);

void apply_bottom_bc(HaloView<Real***> const&             f,
                     ConstantDirichletBC const&           bc,
                     CurvilinearMesh const&               mesh,
                     Kokkos::DefaultExecutionSpace const& stream);

void apply_top_bc(HaloView<Real***> const&             f,
                  ConstantGradientBC const&            bc,
                  CurvilinearMesh const&               mesh,
                  const Kokkos::DefaultExecutionSpace& stream);

void apply_bottom_bc(HaloView<Real***> const&             f,
                     ConstantGradientBC const&            bc,
                     CurvilinearMesh const&               mesh,
                     const Kokkos::DefaultExecutionSpace& stream);

void apply_top_bc(HaloView<Real***> const&             f,
                  ConstantFluxBC const&                bc,
                  Real                                 D,
                  CurvilinearMesh const&               mesh,
                  Kokkos::DefaultExecutionSpace const& stream);

void apply_bottom_bc(HaloView<Real***> const&             f,
                     ConstantFluxBC const&                bc,
                     Real                                 D,
                     CurvilinearMesh const&               mesh,
                     Kokkos::DefaultExecutionSpace const& stream);

} // namespace alps::solver
