#include <solvers/ns/free_surface/solver.h>

#include <solvers/ns/curvilinear_common/rayleigh_damp-inl.h>

namespace alps::solver {

template void RayleighDamp<FreeSurfaceSolver>::add_forces(
  const Vector3Field<Real***>&         fu,
  const Kokkos::DefaultExecutionSpace& space) const;

template void RayleighDampScalar<FreeSurfaceSolver>::add_source(
  const Kokkos::View<Real***, Kokkos::LayoutLeft>& fc,
  const Kokkos::DefaultExecutionSpace&             space) const;

} // namespace alps::solver
