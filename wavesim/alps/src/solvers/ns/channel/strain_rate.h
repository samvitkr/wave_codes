#pragma once

#include <common/container/matrix_field.h>
#include <common/real_type.h>
#include <solvers/mesh/mesh_fwd.h>

namespace alps::solver {

void calculate_strain_rate(const SymmTensor33Field<Real***>& Sij,
                           const Vector3Field<Real***>&      u_vec,
                           const Mesh&                       mesh);

} // namespace alps::solver
