//
// Created by xuanx004 on 7/22/24.
//

#include <solvers/source_terms/boussinesq.h>

#include "solver_ab2.h"

namespace alps::solver {

template<>
void BoussinesqForce<ChannelFlowSolverAB2>::add_forces(
  const Vector3Field<Real***>&         fu,
  const Kokkos::DefaultExecutionSpace& space) const
{
  if (std::abs(Ri) < std::numeric_limits<Real>::epsilon()) return;

  if (theta0 == std::numeric_limits<double>::max()) {
    add_buoyancy_force_without_ref_scalar(
      fu.z, scalar->array, Ri, solver.flow_field.mesh, space);
  } else {
    add_buoyancy_force_with_ref_scalar(
      fu.z, scalar->array, (Real)theta0, Ri, solver.flow_field.mesh, space);
  }
}

} // namespace alps::solver

// Explicit instantiation
namespace alps::solver {
template class BoussinesqForce<ChannelFlowSolverAB2>;
} // namespace alps::solver
