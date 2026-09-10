//
// Created by xuanx004 on 7/16/24.
//

#pragma once

#include <common/container/vector_field.h>
#include <common/real_type.h>
#include <solvers/mesh/curvilinear_mesh.h>

namespace alps::solver {
void add_coriolis_forces(Vector3Field<Real***> const&         fu,
                         HaloView<Real const***> const&       u,
                         HaloView<Real const***> const&       v,
                         Real                                 fz,
                         CurvilinearMesh const&               mesh,
                         Kokkos::DefaultExecutionSpace const& space);
} // namespace alps::solver
