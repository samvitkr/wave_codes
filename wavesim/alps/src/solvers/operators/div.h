#pragma once

#include <common/container/vector_field.h>
#include <common/container/view_types.h>
#include <solvers/mesh/mesh.h>

#include <string_view>

namespace alps {
namespace solver {

/// Calculate the divergence of a vector field u to the cell center on a
/// staggered mesh
/** @note It is expected that vec_u.x and vec_u.y defined on the cell center,
 * and vec_u.z is defined on cell nodes
 */
[[nodiscard]] MDView<Real***> div(const Vector3Field<Real***>& vec_u,
                                  const Mesh&                  mesh,
                                  CenterPt,
                                  std::string_view label = "div_u");

/// Calculate the divergence of a vector field u to the cell nodes on a
/// staggered mesh
/** @note It is expected that vec_u.x and vec_u.y defined on the cell nodes, and
 * vec_u.z is defined on cell centers.
 */
[[nodiscard]] MDView<Real***> div(const Vector3Field<Real***>& vec_u,
                                  const Mesh&                  mesh,
                                  NodePt,
                                  std::string_view label = "div_u");

/// Calculate the maximum divergence of u (halo of u must be updated beforehand)
[[nodiscard]] std::pair<Real, std::array<int, 3>>
max_div(const Vector3Field<Real***>& vec_u, const Mesh& mesh, CenterPt);

} // namespace solver
} // namespace alps
