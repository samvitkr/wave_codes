#pragma once

#include <common/container/vector_field.h>
#include <common/container/view_types.h>
#include <solvers/mesh/curvilinear_mesh.h>

#include <string_view>

namespace alps::solver {

/** @brief Add curvilinear Laplacian fluxes for a cell-centered field.
 *
 * Calculates and adds F_j = J^{-1} g^{ij} df/dξ_i to fluxes. Inputs fluxes.x
 * and fluxes.y are defined on cell centers, and fluxes.z is defined on cell
 * nodes.
 */
void add_laplacian_fluxes(const Vector3Field<Real***>&   fluxes,
                          const HaloView<Real const***>& f,
                          const BottomWaveMesh&          mesh,
                          Real                           coeff,
                          std::string_view               label = "lap flux");

/** @brief Add curvilinear Laplacian fluxes for a cell-centered field, excluding
 * g33 contributions.
 */
void add_laplacian_fluxes_no_g33(const Vector3Field<Real***>&   fluxes,
                                 const HaloView<Real const***>& f,
                                 const BottomWaveMesh&          mesh,
                                 Real                           coeff,
                                 std::string_view label = "lap flux");

/** @brief Add curvilinear Laplacian fluxes for a cell-centered field.
 *
 * Calculates and adds F_j = J^{-1} g^{ij} df/dξ_i to fluxes. Inputs fluxes.x
 * and fluxes.y are defined on cell centers, and fluxes.z is defined on cell
 * nodes.
 */
void add_laplacian_fluxes(const Vector3Field<Real***>&   fluxes,
                          const HaloView<Real const***>& f,
                          const TopWaveMesh&             mesh,
                          Real                           coeff,
                          std::string_view               label = "lap flux");

} // namespace alps::solver
