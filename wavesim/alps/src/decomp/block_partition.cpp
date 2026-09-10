#include "block_partition.h"

#include <Kokkos_Array.hpp>
#include <mpipp/collectives.h>

#include <stdexcept>
#include <vector>

namespace alps {

namespace {
/// Return a distribution of dim over n blocks
std::vector<int> distribute(int dim, int n)
{
  if (n <= 0) {
    throw std::invalid_argument("Number of blocks must be positive.");
  }
  if (dim <= 0) {
    throw std::invalid_argument("The dimension size must be positive.");
  }

  std::vector<int> distribution(n, dim / n);
  for (int i = 0; i < dim % n; ++i) {
    distribution[i] += 1;
  }
  return distribution;
}

int calc_offset(std::vector<int> distribution, int n)
{
  int offset = 0;
  for (int i = 0; i < n; ++i) {
    offset += distribution[i];
  }
  return offset;
}
} // anonymous namespace

BlockPartition::BlockPartition(const CommType&  mdComm,
                               std::vector<int> grid_size)
  : comm{mdComm}
  , extents{distribute(grid_size.at(0), mdComm.dims[0])[mdComm.coords[0]],
            distribute(grid_size.at(1), mdComm.dims[1])[mdComm.coords[1]],
            distribute(grid_size.at(2), mdComm.dims[2])[mdComm.coords[2]]}
  , offsets{calc_offset(distribute(grid_size.at(0), mdComm.dims[0]),
                        mdComm.coords[0]),
            calc_offset(distribute(grid_size.at(1), mdComm.dims[1]),
                        mdComm.coords[1]),
            calc_offset(distribute(grid_size.at(2), mdComm.dims[2]),
                        mdComm.coords[2])}
  , global_extents{grid_size.at(0), grid_size.at(1), grid_size.at(2)}
{}

std::vector<int> BlockPartition::get_distribution(int axis) const
{
  const auto&      axisComm = comm.axis_comm[axis];
  std::vector<int> distribution(axisComm.size());
  mpipp::allgather(extents[axis], distribution.data(), axisComm);
  return distribution;
}

std::vector<int> BlockPartition::get_offsets(int axis) const
{
  const auto&      axisComm = comm.axis_comm[axis];
  std::vector<int> axis_offsets(axisComm.size());
  mpipp::allgather(offsets[axis], axis_offsets.data(), axisComm);
  return axis_offsets;
}

} // namespace alps
