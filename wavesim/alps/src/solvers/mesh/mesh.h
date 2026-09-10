#pragma once

#include "mesh_fwd.h"
#include <common/container/view_types.h>
#include <common/math.h>
#include <common/real_type.h>
#include <io/io_fwd.h>
#include <spectral/spectral_fwd.h>

#include <Kokkos_Macros.hpp>
#include <nonstd/span.hpp>

#include <filesystem>

namespace alps {
struct MPIComm3D;
struct BlockPartition;
} // namespace alps

namespace alps::solver {

/**
 * @brief Represents a 3D computational mesh.
 *
 * This class represents a 3D mesh, defined on a given SpectralGrid. The
 * vertical coordinates distributed on each processor and the ghost cells are
 * stored. It provides methods to load and save the mesh to a file, and to
 * access information about the mesh, such as its size and partitioning.
 *
 * @note This class is not default-constructible, and it must be initialized
 * with a grid and optional parameters.
 */
class Mesh
{
 public:
  Grid const& grid;

  /// @brief The base wavenumber of the mesh in the x and y directions
  Real pex, pey;

  /// @brief The length of the mesh in the z direction
  Real hbar;

  /// @brief The computational coordinates ζ defined on cell centers
  HaloView<Real*> zz;

  /// @brief The spacings of computational coordinates ζ between cell centers
  HaloView<Real*> dz;

  /// @brief The computational coordinates ζ defined on cell nodes
  HaloView<Real*> zw;

  /// @brief The spacings of computational coordinates ζ between cell nodes
  HaloView<Real*> dzw;

  /// @brief The computational coordinates on the host
  HaloView<Real*, Kokkos::HostSpace> zz_h, dz_h, zw_h, dzw_h;

  Mesh() = delete;

  Mesh(Grid const& grid_, Real Lz, int n_ghost = 1);

  /// Load a global mesh, the input should be the global mesh
  void load(nonstd::span<Real const> z_node) const;

  void save_to_file(std::filesystem::path filename) const;

  void save_to_file(HighFive::File& hdf_file) const;

  /// @brief Read mesh from file by invoking @ref read_and_interp
  void read_from_file(std::filesystem::path filename) const;

  /// @brief Read in `zw` from a file and populate other arrays of the mesh
  void read_and_interp(HighFive::File const& hdf_file) const;

  const MPIComm3D& comm() const;

  const BlockPartition& partition() const;

  /// @brief The size of the grid block locally on the processor in each axis
  /// (alias of local_extent)
  int extent(int axis) const;

  /// @brief The size of the grid block locally on the processor in each axis
  int local_extent(int axis) const;

  /// @brief The global grid size in each axis
  int global_extent(int axis) const;

  /// @brief The shape of the grid block locally on the processor (alias of
  /// local_extents)
  Kokkos::Array<int, 3> extents() const;

  /// @brief The shape of the grid block locally on the processor
  Kokkos::Array<int, 3> local_extents() const;

  /// @brief The shape of the global grid
  Kokkos::Array<int, 3> global_extents() const;
};

/// @brief Interpolate the cell-centered values v0 and v1 to the node
/// in-between, where d1 and d2 are the spacings of the cells where v0 and v1
/// are located.
template<class T, class S>
KOKKOS_FORCEINLINE_FUNCTION auto itp2node(T v0, T v1, S d0, S d1)
{
  auto beta = d0 / (d0 + d1);
  return Kokkos::fma(beta, v1, Kokkos::fma(-beta, v0, v0));
}

/// @brief Interpolate the cell-centered values v0 and v1 to the node
/// in-between, where ratio is the ratio of the cell spacing of v0 to the
/// spacing of v1.
/** v = v0*(1-ratio) + v1*ratio */
template<class T, class S>
KOKKOS_FORCEINLINE_FUNCTION auto itp2node(T v0, T v1, S ratio)
{
  return Kokkos::fma(ratio, v1, Kokkos::fma(-ratio, v0, v0));
}

/// @brief Interpolate the node-centered values v0 and v1 to the cell center
template<class T>
KOKKOS_FORCEINLINE_FUNCTION auto itp2center(T v0, T v1)
{
  return (v0 + v1) / 2;
}
} // namespace alps::solver
