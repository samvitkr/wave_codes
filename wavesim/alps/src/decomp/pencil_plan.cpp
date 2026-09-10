#include "pencil_plan.h"

#include "block_partition.h"
#include "mdcomm.h"
#include <common/real_type.h>
#include <transpose/transposer_pool.h>

#include <array>
#include <vector>

namespace alps {

PencilPlan::PencilPlan(const MPIComm3D& mdComm, std::vector<int> grid_size)
  : x_pencil(mdComm, grid_size)
  , y_pencil(mdComm, {grid_size[1], grid_size[0], grid_size[2]})
{
  // Create transposers of common sizes
  auto& pool =
    transpose::TransposerPool<Real,
                              Kokkos::DefaultExecutionSpace>::get_instance();

  pool.get_transposer(x_pencil.comm.axis_comm[1],
                      x_pencil.global_extents[0],
                      x_pencil.global_extents[1],
                      extent(2));
  pool.get_transposer(y_pencil.comm.axis_comm[1],
                      y_pencil.global_extents[0],
                      y_pencil.global_extents[1],
                      extent(2));
}

Kokkos::LayoutLeft create_local_layout(const PencilPlan&  plan,
                                       std::array<int, 3> n_halos,
                                       Pencil             which)
{
  auto local_extents = plan.extents(which);
  return Kokkos::LayoutLeft(local_extents[0] + 2 * n_halos[0],
                            local_extents[1] + 2 * n_halos[1],
                            local_extents[2] + 2 * n_halos[2]);
}

Kokkos::LayoutLeft create_local_layout(const PencilPlan& plan, Pencil which)
{
  return create_local_layout(plan, {0, 0, 0}, which);
}
} // namespace alps
