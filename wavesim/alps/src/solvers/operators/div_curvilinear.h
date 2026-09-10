#pragma once

#include <common/container/vector_field.h>
#include <common/container/view_types.h>
#include <solvers/mesh/curvilinear_mesh.h>

#include <string_view>

namespace alps {
namespace solver {

/** @brief Calculate the divergence of u, with input J^{-1}u
 *
 * Calculate ∂[J^{-1}(∂ξ_k/∂x_i)u_i]/∂ξ_k, which is actually J∇·u
 * @note Inputs vec_invJu.x and vec_invJu.y should be defined on the cell
 * center, and vec_invJu.z is defined on cell nodes
 */
[[nodiscard]] MDView<Real***> div(const Vector3Field<Real***>& vec_invJu,
                                  const BottomWaveMesh&        mesh,
                                  CenterPt,
                                  std::string_view label = "div_u");

/** @brief Calculate the divergence of u, with input J^{-1}u
 *
 * Calculate ∂[J^{-1}(∂ξ_k/∂x_i)u_i]/∂ξ_k, which is actually J∇·u
 * @note Inputs vec_invJu.x and vec_invJu.y should be defined on the cell
 * center, and vec_invJu.z is defined on cell nodes
 */
[[nodiscard]] MDView<Real***> div(const Vector3Field<Real***>& vec_invJu,
                                  const TopWaveMesh&           mesh,
                                  CenterPt,
                                  std::string_view label = "div_u");

/** @brief Calculate the maximum divergence of u, as well as the local indices
 * where the maximum is (halo of u must be updated beforehand)
 *
 * Calculate J^{-1}∇·u = ∂[J^{-1}(∂ξ_k/∂x_i)u_i]/∂ξ_k
 * @note Inputs vec_u.x and vec_u.y should be defined on the cell center,
 * and vec_u.z is defined on cell nodes
 */
[[nodiscard]] std::pair<Real, std::array<int, 3>>
max_div(const Vector3Field<Real***>& vec_u,
        const BottomWaveMesh&        mesh,
        CenterPt /*unused*/);

/** @brief Calculate the maximum divergence of u, as well as the local indices
 * where the maximum is (halo of u must be updated beforehand)
 *
 * Calculate J^{-1}∇·u = ∂[J^{-1}(∂ξ_k/∂x_i)u_i]/∂ξ_k
 * @note Inputs vec_u.x and vec_u.y should be defined on the cell center,
 * and vec_u.z is defined on cell nodes
 */
[[nodiscard]] std::pair<Real, std::array<int, 3>>
max_div(const Vector3Field<Real***>& vec_u,
        const TopWaveMesh&           mesh,
        CenterPt /*unused*/);

} // namespace solver
} // namespace alps
