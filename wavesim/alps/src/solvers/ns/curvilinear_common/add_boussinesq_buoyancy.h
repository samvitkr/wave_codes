//
// Created by xuanx004 on 7/24/24.
//

#pragma once

#include <common/container/vector_field.h>
#include <common/real_type.h>
#include <solvers/mesh/curvilinear_mesh.h>

namespace alps::solver {

void add_buoyancy_force_with_ref_scalar_curve(
  HaloView<Real***> const&             fz,
  HaloView<Real const***> const&       theta,
  Real                                 theta0,
  Real                                 beta,
  CurvilinearMesh const&               mesh,
  Kokkos::DefaultExecutionSpace const& stream);

void add_buoyancy_force_without_ref_scalar_curve(
  HaloView<Real***> const&             fz,
  HaloView<Real const***> const&       theta,
  Real                                 beta,
  CurvilinearMesh const&               mesh,
  Kokkos::DefaultExecutionSpace const& stream);

} // namespace alps::solver
