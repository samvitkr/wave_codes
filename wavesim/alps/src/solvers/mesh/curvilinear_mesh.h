#pragma once

#include "mesh.h"
#include <common/container/view_types.h>

#include <Kokkos_Macros.hpp>

namespace alps {
namespace solver {

/** @brief Represents a 3D curvilinear computational mesh.
 * The mesh is transformed in the vertical direction
 *
 * @note CurvilinearMesh is not default-constructible
 */
class CurvilinearMesh : public Mesh
{
 public:
  MDView<Real**> ex; ///< eta_x
  MDView<Real**> ey; ///< eta_y

  MDView<Real**> J;    ///< Jacobian
  MDView<Real**> invJ; ///< Inverse Jacobian J^{-1}

  MDView<Real**> exr; ///< eta_x*J
  MDView<Real**> eyr; ///< eta_y*J

  MDView<Real**> et; ///< eta_t

  CurvilinearMesh() = delete;

  CurvilinearMesh(Grid const& grid_, Real Lz, int n_ghost = 1);
};

/// @brief A 3D curvilinear computational mesh with a bottom wave
class BottomWaveMesh : public CurvilinearMesh
{
 public:
  using CurvilinearMesh::CurvilinearMesh;

  /// Update the metric coefficients in the curvilinear mesh given eta.
  void
  update_metric_coefficients(MDView<Real const**> const&          eta,
                             Kokkos::DefaultExecutionSpace const& space) const;

  /// Compute ζₓ
  template<class T, class S>
  static KOKKOS_FORCEINLINE_FUNCTION auto zeta_x(T zeta, S exr)
  {
    return (zeta - 1) * exr;
  }

  /// Compute ζ_y
  template<class T, class S>
  static KOKKOS_FORCEINLINE_FUNCTION auto zeta_y(T zeta, S eyr)
  {
    return zeta_x(zeta, eyr);
  }

  /// Compute ζ_z
  template<class T, class S>
  static KOKKOS_FORCEINLINE_FUNCTION auto zeta_z(T /*zeta*/, S invJ)
  {
    return 1 / invJ;
  }

  /// Compute J^{-1} ζ_x
  template<class T, class S>
  static KOKKOS_FORCEINLINE_FUNCTION auto invJ_zeta_x(T zeta, S eta_x)
  {
    return (zeta - 1) * eta_x;
  }

  /// @brief Compute J^{-1} ζ_y
  template<class T, class S>
  static KOKKOS_FORCEINLINE_FUNCTION auto invJ_zeta_y(T zeta, S eta_y)
  {
    return invJ_zeta_x(zeta, eta_y);
  }

  /// @brief Compute J^{-1} ζ_z
  static KOKKOS_FORCEINLINE_FUNCTION auto invJ_zeta_z() { return 1; }

  /// @brief Compute J^{-1} W, where W is the contravariant velocity in the ζ
  /// direction
  template<class T, class S>
  static KOKKOS_FORCEINLINE_FUNCTION auto invJWg(T zeta, S eta_t)
  {
    return (1 - zeta) * eta_t;
  }

  template<class T, class S>
  static KOKKOS_FORCEINLINE_FUNCTION auto g33(T zeta, S exr, S eyr, S invJ)
  {
    auto zx  = zeta_x(zeta, exr);
    auto zy  = zeta_y(zeta, eyr);
    auto ztz = zeta_z(zeta, invJ);
    return zx * zx + zy * zy + ztz * ztz;
  }
};

/// @brief A 3D curvilinear computational mesh with a top wave
class TopWaveMesh : public CurvilinearMesh
{
 public:
  using CurvilinearMesh::CurvilinearMesh;

  /// Update the metric coefficients in the curvilinear mesh given eta.
  void
  update_metric_coefficients(MDView<Real const**> const&          eta,
                             Kokkos::DefaultExecutionSpace const& space) const;

  /// Compute ζₓ
  template<class T, class S>
  static KOKKOS_FORCEINLINE_FUNCTION auto zeta_x(T zeta, S exr)
  {
    return -zeta * exr;
  }

  /// Compute ζ_y
  template<class T, class S>
  static KOKKOS_FORCEINLINE_FUNCTION auto zeta_y(T zeta, S eyr)
  {
    return zeta_x(zeta, eyr);
  }

  /// Compute ζ_z
  template<class T, class S>
  static KOKKOS_FORCEINLINE_FUNCTION auto zeta_z(T /*zeta*/, S invJ)
  {
    return 1 / invJ;
  }

  /// Compute J^{-1} ζ_x
  template<class T, class S>
  static KOKKOS_FORCEINLINE_FUNCTION auto invJ_zeta_x(T zeta, S eta_x)
  {
    return -zeta * eta_x;
  }

  /// @brief Compute J^{-1} ζ_y
  template<class T, class S>
  static KOKKOS_FORCEINLINE_FUNCTION auto invJ_zeta_y(T zeta, S eta_y)
  {
    return invJ_zeta_x(zeta, eta_y);
  }

  /// @brief Compute J^{-1} ζ_z
  static KOKKOS_FORCEINLINE_FUNCTION auto invJ_zeta_z() { return 1; }

  /// @brief Compute J^{-1} W, where W is the contravariant velocity in the ζ
  /// direction
  template<class T, class S>
  static KOKKOS_FORCEINLINE_FUNCTION auto invJWg(T zeta, S eta_t)
  {
    return zeta * eta_t;
  }

  template<class T, class S>
  static KOKKOS_FORCEINLINE_FUNCTION auto g33(T zeta, S exr, S eyr, S invJ)
  {
    auto zx  = zeta_x(zeta, exr);
    auto zy  = zeta_y(zeta, eyr);
    auto ztz = zeta_z(zeta, invJ);
    return zx * zx + zy * zy + ztz * ztz;
  }
};

} // namespace solver
} // namespace alps
