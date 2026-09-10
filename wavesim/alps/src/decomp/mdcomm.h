#pragma once

#include <mpipp/comm.h>

#include <array>
#include <vector>

namespace alps {

/** \brief Three-dimensional Cartesian communicator
 *
 * An MPI communicator representing a processor grid in Fortran ordering, i.e.
 * the rank stride of the leftmost coordinate is 1.
 * Design notes: copying the struct shallow copies the underlying MPI
 * communicators. The communicators are not destroyed explicitly, which should
 * be fine because the resources are cleaned up by MPI_Finalize.
 */
struct MPIComm3D
{
 private:
  template<class T>
  using Array = std::array<T, 3>;

 public:
  /// Cartesian communicator
  mpipp::communicator comm;

  /// Axis communicators for the process where the class is constructed
  Array<mpipp::communicator> axis_comm;

  /// Number of processes in each axis
  Array<int> dims;

  /// Coordinates along each axis for the process where the class is constructed
  Array<int> coords;

  /// Rank of the next process of the calling process along each axis
  /** Stores the destination rank obtained from MPI_Cart_shift with disp=1 */
  Array<int> next_proc_on_axis;

  /// Rank of the previous process of the calling process along each axis
  /** Stores the source rank obtained from MPI_Cart_shift with disp=1 */
  Array<int> prev_proc_on_axis;

  /// Constructor with communicator and axis sizes
  MPIComm3D(mpipp::communicator const& base_comm,
            std::vector<int>           comm_dims,
            std::vector<int>           periodic);

  [[nodiscard]] int rank() const noexcept { return comm.rank(); }

  [[nodiscard]] int size() const noexcept { return comm.size(); }

  MPI_Comm raw_handle() const noexcept { return comm.raw_handle(); }

  operator mpipp::communicator&() { return comm; }

  operator mpipp::communicator const&() const { return comm; }

  /// Query if the calling process is the first along an axis
  bool is_first(int axis) const { return coords.at(axis) == 0; }

  /// Query if the calling process is the last along an axis
  bool is_last(int axis) const { return coords.at(axis) == dims.at(axis) - 1; }

  /// Get the global rank in the communicator with a given coordinate
  int to_rank(std::vector<int> proc_coords) const noexcept;

  /// Get the coordinate in the communicator with a given global rank
  Array<int> to_coordinates(int rk) const noexcept;

  /// Get the shifted index of the calling process along a given axis
  std::pair<int, int> shifted_ranks(int dim, int disp) const noexcept;
};
} // namespace alps
