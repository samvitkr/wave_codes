#pragma once

#include "block_partition.h"
#include "mdcomm.h"

#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE
#include <Kokkos_Layout.hpp>
#undef KOKKOS_IMPL_PUBLIC_INCLUDE
#endif

#include <array>
#include <vector>

namespace alps {

enum class Pencil
{
  X,
  Y
};

/// @brief An aggregate of different pencil partitions
class PencilPlan
{
 public:
  PencilPlan(const MPIComm3D& mdComm, std::vector<int> grid_size);

  /// Return the corresponding pencil partition
  const auto& partition(Pencil which = Pencil::X) const noexcept
  {
    return which == Pencil::Y ? y_pencil : x_pencil;
  }

  /// Size of one global dimension
  auto global_extent(int axis, Pencil which = Pencil::X) const noexcept
  {
    return partition(which).global_extents[axis];
  }

  /// Size of one dimension of the local block
  auto extent(int axis, Pencil which = Pencil::X) const noexcept
  {
    return partition(which).extents[axis];
  }

  /// Sizes of the local block in all dimensions
  const auto& extents(Pencil which = Pencil::X) const noexcept
  {
    return partition(which).extents;
  }

  /// Global index offset of the local block for a specified axis
  auto offset(int axis, Pencil which = Pencil::X) const noexcept
  {
    return partition(which).offsets[axis];
  }

  /// Global index offsets of the local block in all dimensions
  const auto& offsets(Pencil which = Pencil::X) const noexcept
  {
    return partition(which).offsets;
  }

  const auto& comm() const noexcept { return x_pencil.comm; }

  // Partition of the x-pencil and y-pencil
  BlockPartition x_pencil;
  BlockPartition y_pencil;
};

Kokkos::LayoutLeft create_local_layout(const PencilPlan&  plan,
                                       std::array<int, 3> n_halos,
                                       Pencil             which = Pencil::X);

Kokkos::LayoutLeft create_local_layout(const PencilPlan& plan,
                                       Pencil            which = Pencil::X);

} // namespace alps
