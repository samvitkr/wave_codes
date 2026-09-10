//
// Created by xuanx004 on 7/16/24.
//

#include "add_coriolis_force.h"

#include <common/container/view_utils.h>
#include <common/kokkos_abstraction/exec_policy.h>

namespace alps::solver {

void add_coriolis_forces(Vector3Field<Real***> const&         fu,
                         HaloView<Real const***> const&       u,
                         HaloView<Real const***> const&       v,
                         Real const                           fz,
                         CurvilinearMesh const&               mesh,
                         const Kokkos::DefaultExecutionSpace& space)
{
  if (std::abs(fz) < std::numeric_limits<Real>::epsilon()) return;

  auto const& fux = fu.x;
  auto const& fuy = fu.y;

  const auto& invJ = mesh.invJ;
  Kokkos::parallel_for(
    "add CoriolisForce",
    LoopPolicy<3>(space, local_begins(fux), local_ends(fux)),
    KOKKOS_LAMBDA(int i, int j, int k) {
      fux(i, j, k) += fz * v(i, j, k) * invJ(i, j);
      fuy(i, j, k) -= fz * u(i, j, k) * invJ(i, j);
    });
}

} // namespace alps::solver
