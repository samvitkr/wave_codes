#pragma once

#include "mdcomm.h"

#include <Kokkos_Array.hpp>

#include <vector>

namespace alps {

/** \brief Describes the block partition of a 3D grid over a 3D process grid
 */
struct BlockPartition
{
 private:
  using dims_t   = Kokkos::Array<int, 3>;
  using CommType = MPIComm3D;

 public:
  CommType comm;

  dims_t extents{0};

  dims_t offsets{0};

  dims_t global_extents{0};

  BlockPartition(const CommType& mdComm, std::vector<int> grid_size);

  BlockPartition(const BlockPartition&) = default;

  /// Get how the each dimension is distributed along the specified axis
  [[nodiscard]] std::vector<int> get_distribution(int axis) const;

  /// Get how all the offsets in each dimension
  [[nodiscard]] std::vector<int> get_offsets(int axis) const;
};

} // namespace alps
