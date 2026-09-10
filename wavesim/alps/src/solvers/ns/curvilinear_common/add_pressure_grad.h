//
// Created by xuanx004 on 7/16/24.
//

#pragma once

#include <common/container/vector_field.h>
#include <solvers/mesh/curvilinear_mesh.h>

#include <array>

namespace alps::solver {
void add_pressure_gradient(const Vector3Field<Real***>& source,
                           const CurvilinearMesh&       mesh,
                           std::array<Real, 2>          gradients);
} // namespace alps::solver
