#include <solvers/source_terms/pressure_grad.h>

#include <solvers/ns/curvilinear_common/add_pressure_grad.h>
#include <solvers/ns/flow_over_wave/solver.h>

namespace alps {
namespace solver {

template<>
void ConstantPressureGradient<FlowOverWaveSolver>::add_forces(
  const Vector3Field<Real***>&         fu,
  const Kokkos::DefaultExecutionSpace& space) const
{
  space.fence();
  add_pressure_gradient(fu, solver_.flow_field.mesh, gradients_);
}

} // namespace solver
} // namespace alps
