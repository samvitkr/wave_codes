#pragma once

#include <common/container/vector_field.h>
#include <common/container/view_types.h>
#include <solvers/mesh/curvilinear_mesh.h>

namespace alps::solver {

/// Calculate the grad of f, where f at cell center
[[nodiscard]] Vector3Field<Real***>
grad(const HaloView<Real***>& f, const BottomWaveMesh& mesh, CenterPt tag);

/// Calculate the grad of f, where f at cell center
[[nodiscard]] Vector3Field<Real***>
grad(const HaloView<Real***>& f, const TopWaveMesh& mesh, CenterPt tag);

/// Calculate the grad of f, where f at cell node
[[nodiscard]] Vector3Field<Real***>
grad(const HaloView<Real***>& f, const BottomWaveMesh& mesh, NodePt tag);

/// Calculate the grad of f, where f at cell node
[[nodiscard]] Vector3Field<Real***>
grad(const HaloView<Real***>& f, const TopWaveMesh& mesh, NodePt tag);

/** @brief Calculate the gradient of a field f, where f is defined on cell
 * centers
 *
 * ∂u/∂x and ∂u/∂y are stored in grad_f.x and grad_f.y and are defined on cell
 * centers, ∂u/∂z is stored in grad_f.z and is defined on cell nodes.
 */
void grad(const Vector3Field<Real***>& grad_f,
          const HaloView<Real***>&     f,
          const BottomWaveMesh&        mesh,
          CenterPt                     tag);

/** @brief Calculate the gradient of a field f, where f is defined on cell
 * centers
 *
 * ∂u/∂x and ∂u/∂y are stored in grad_f.x and grad_f.y and are defined on cell
 * centers, ∂u/∂z is stored in grad_f.z and is defined on cell nodes.
 */
void grad(const Vector3Field<Real***>& grad_f,
          const HaloView<Real***>&     f,
          const TopWaveMesh&           mesh,
          CenterPt                     tag);

/** @brief Calculate the gradient of a field f, where f is defined on cell nodes
 *
 * ∂w/∂x and ∂w/∂y are stored in grad_f.x and grad_f.y and are defined on cell
 * nodes, ∂u/∂z is stored in grad_f.z and is defined on cell centers.
 */
void grad(const Vector3Field<Real***>& grad_f,
          const HaloView<Real***>&     f,
          const BottomWaveMesh&        mesh,
          NodePt                       tag);

/** @brief Calculate the gradient of a field f, where f is defined on cell nodes
 *
 * ∂w/∂x and ∂w/∂y are stored in grad_f.x and grad_f.y and are defined on cell
 * nodes, ∂u/∂z is stored in grad_f.z and is defined on cell centers.
 */
void grad(const Vector3Field<Real***>& grad_f,
          const HaloView<Real***>&     f,
          const TopWaveMesh&           mesh,
          NodePt                       tag);

} // namespace alps::solver
