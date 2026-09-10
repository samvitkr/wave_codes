//
// Created by xuanx004 on 12/30/23.
//

#include <solvers/source_terms/coriolis.h>

#include <solvers/ns/curvilinear_common/add_coriolis_force.h>
#include <solvers/ns/free_surface/solver.h>

namespace alps::solver {

template<>
void CoriolisForce<FreeSurfaceSolver>::add_forces(
  const Vector3Field<Real***>&         fu,
  const Kokkos::DefaultExecutionSpace& space) const
{
  add_coriolis_forces(fu,
                      solver.flow_field.u.x,
                      solver.flow_field.u.y,
                      this->fz,
                      solver.flow_field.mesh,
                      space);
}

} // namespace alps::solver
