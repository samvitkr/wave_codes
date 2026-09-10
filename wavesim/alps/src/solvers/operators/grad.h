#pragma once

#include <common/container/vector_field.h>
#include <common/container/view_types.h>
#include <solvers/mesh/mesh.h>

#include <Kokkos_Core.hpp>

namespace alps {
namespace solver {

/// Calculate the gradient of a field u, where u is defined on cell centers
/** (see the in-place version for details) */
[[nodiscard]] Vector3Field<Real***>
grad(const HaloView<Real***>& u, const Mesh& mesh, CenterPt);

/// Calculate the gradient of a field w, where w is defined on cell nodes
/** (see the in-place version for details) */
[[nodiscard]] Vector3Field<Real***>
grad(const HaloView<Real***>& w, const Mesh& mesh, NodePt);

/// Calculate the gradient of a field u, where u is defined on cell centers
/**
 * ∂u/∂x and ∂u/∂y are stored in grad_u.x and grad_u.y and are defined on cell
 * centers, ∂u/∂z is stored in grad_u.z and is defined on cell nodes.
 */
void grad(const Vector3Field<Real***>& grad_u,
          const HaloView<Real***>&     u,
          const Mesh&                  mesh,
          CenterPt /*tag*/);

/// Calculate the gradient of a field w, where w is defined on cell nodes
/**
 * ∂w/∂x and ∂w/∂y are stored in grad_u.x and grad_u.y and are defined on cell
 * nodes, ∂u/∂z is stored in grad_u.z and is defined on cell centers.
 */
void grad(const Vector3Field<Real***>& grad_u,
          const HaloView<Real***>&     w,
          const Mesh&                  mesh,
          NodePt /*tag*/);

struct LeftBoundary
{};
struct RightBoundary
{};

template<class T, class ExecSpace>
struct FDBoundaryFunctor
{
  struct Add
  {};

  FDBoundaryFunctor(MDView<T**, ExecSpace>                output,
                    MDView<T const** [3], ExecSpace>      input_stencil,
                    MDView<T const[2], Kokkos::HostSpace> delta,
                    T                                     Lz,
                    LeftBoundary                          tag);

  FDBoundaryFunctor(MDView<T**, ExecSpace>                output,
                    MDView<T const** [3], ExecSpace>      input_stencil,
                    MDView<T const[2], Kokkos::HostSpace> delta,
                    T                                     Lz,
                    RightBoundary                         tag);

  MDView<T**, ExecSpace>           fz;
  MDView<T const** [3], ExecSpace> f;
  T                                coeff[3];

  KOKKOS_INLINE_FUNCTION void operator()(int i, int j) const
  {
    fz(i, j) =
      (coeff[0] * f(i, j, 0) + coeff[1] * f(i, j, 1) + coeff[2] * f(i, j, 2));
  }

  KOKKOS_INLINE_FUNCTION void operator()(const Add& /*tag*/, int i, int j) const
  {
    fz(i, j) +=
      (coeff[0] * f(i, j, 0) + coeff[1] * f(i, j, 1) + coeff[2] * f(i, j, 2));
  }
};

} // namespace solver
} // namespace alps
