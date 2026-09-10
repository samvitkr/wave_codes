//
// Created by xuanx004 on 12/30/23.
//

#include <solvers/source_terms/coriolis.h>

#include <common/container/view_utils.h>
#include <common/kokkos_abstraction/exec_policy.h>
#include <solvers/ns/channel/solver.h>

namespace alps::solver {

namespace {
void add_coriolis_forces(const Vector3Field<Real***>&         fu,
                         const HaloView<const Real***>&       u,
                         const HaloView<const Real***>&       v,
                         Real                                 fz,
                         const Kokkos::DefaultExecutionSpace& space)
{
  auto const& fux = fu.x;
  auto const& fuy = fu.y;

  const auto f = fz;
  Kokkos::parallel_for(
    "add CoriolisForce",
    LoopPolicy<3>(space, local_begins(fux), local_ends(fux)),
    KOKKOS_LAMBDA(int i, int j, int k) {
      fux(i, j, k) += f * v(i, j, k);
      fuy(i, j, k) -= f * u(i, j, k);
    });
}
} // anonymous namespace

template<>
void CoriolisForce<ChannelFlowSolverAB2>::add_forces(
  const Vector3Field<Real***>&         fu,
  const Kokkos::DefaultExecutionSpace& space) const
{
  if (std::abs(this->fz) < std::numeric_limits<Real>::epsilon()) return;

  add_coriolis_forces(
    fu, solver.flow_field.u.x, solver.flow_field.u.y, this->fz, space);
}

template<>
void CoriolisForce<ChannelFlowSolverAB2CN>::add_forces(
  const Vector3Field<Real***>&         fu,
  const Kokkos::DefaultExecutionSpace& space) const
{
  if (std::abs(this->fz) < std::numeric_limits<Real>::epsilon()) return;

  add_coriolis_forces(
    fu, solver.flow_field.u.x, solver.flow_field.u.y, this->fz, space);
}
} // namespace alps::solver
