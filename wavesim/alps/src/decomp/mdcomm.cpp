#include "mdcomm.h"

#include <fmt/format.h>
#include <mpipp/comm.h>
#include <mpipp/environment.h>

#include <stdexcept>

namespace alps {
/**
 * @brief Compute the strides of an array size in Fortran ordering
 */
std::vector<int> compute_strides(std::vector<int> dimensions);

mpipp::communicator slice_communicator(mpipp::communicator const& comm,
                                       std::vector<bool>          remain_dims);

MPIComm3D::MPIComm3D(const mpipp::communicator& base_comm,
                     const std::vector<int>     comm_dims,
                     const std::vector<int>     periodic)
  : comm{[base_comm, &comm_dims, &periodic]() -> mpipp::communicator {
    if (comm_dims.size() != 3 || periodic.size() != 3) {
      throw std::invalid_argument("Invalid processor grid size: incorrect "
                                  "number of axes.");
    }
    Array<int> pdims{};
    Array<int> periodic_flags{};
    int        n_processes{1};
    for (int i = 0; i < 3; ++i) {
      pdims[i] = comm_dims[3 - i - 1];
      n_processes *= pdims[i];
      periodic_flags[i] = periodic[3 - i - 1];
    } // reverse dims and periodic to C-order

    if (base_comm.size() != n_processes) {
      auto msg = fmt::format("Invalid processor grid size: cannot build a grid "
                             "of {}x{}x{} with {} processes.",
                             comm_dims[0],
                             comm_dims[1],
                             comm_dims[2],
                             base_comm.size());
      throw std::invalid_argument(msg);
    }

    MPI_Comm new_comm{MPI_COMM_NULL};
    MPI_Cart_create(base_comm.raw_handle(),
                    int(pdims.size()),
                    pdims.data(),
                    periodic_flags.data(),
                    0,
                    &new_comm);

    return mpipp::communicator::create_by_take_over(new_comm);
  }()}
  , axis_comm{slice_communicator(comm, {true, false, false}),
              slice_communicator(comm, {false, true, false}),
              slice_communicator(comm, {false, false, true})}
  , dims{comm_dims[0], comm_dims[1], comm_dims[2]}
  , coords{to_coordinates(comm.rank())}
  , next_proc_on_axis{shifted_ranks(0, 1).second,
                      shifted_ranks(1, 1).second,
                      shifted_ranks(2, 1).second}
  , prev_proc_on_axis{shifted_ranks(0, 1).first,
                      shifted_ranks(1, 1).first,
                      shifted_ranks(2, 1).first}
{}

int MPIComm3D::to_rank(const std::vector<int> proc_coords) const noexcept
{
  Array<int> c_coords;
  for (int i = 0; i < 3; ++i) {
    c_coords[3 - i - 1] = proc_coords[i];
  }

  int r = 0;
  MPI_Cart_rank(comm.raw_handle(), c_coords.data(), &r);
  return r;
}

MPIComm3D::Array<int> MPIComm3D::to_coordinates(int rk) const noexcept
{
  Array<int> c_coords, proc_coords;
  MPI_Cart_coords(comm.raw_handle(), rk, 3, c_coords.data());
  for (int i = 0; i < 3; ++i) {
    proc_coords[3 - i - 1] = c_coords[i];
  }
  return proc_coords;
}

std::vector<int> compute_strides(const std::vector<int> dimensions)
{
  auto             n = dimensions.size();
  std::vector<int> strides(n);

  strides.at(0) = 1;
  for (decltype(n) axis = 1; axis < n; ++axis) {
    strides[axis] = strides[axis - 1] * dimensions[axis - 1];
  }
  return strides;
}

mpipp::communicator slice_communicator(mpipp::communicator const& comm,
                                       std::vector<bool>          remain_dims)
{
  MPI_Comm         new_sub_comm{MPI_COMM_NULL};
  std::vector<int> c_remain_dims(3, int(false));
  for (int i = 0; i < 3; ++i) {
    c_remain_dims[i] = int(remain_dims[3 - i - 1]);
  }
  MPI_Cart_sub(comm.raw_handle(), c_remain_dims.data(), &new_sub_comm);
  if (new_sub_comm == MPI_COMM_NULL) {
    throw std::runtime_error("Failed to create axis communicators.");
  }
  return mpipp::communicator::create_by_take_over(new_sub_comm);
}

static_assert(std::is_copy_constructible_v<MPIComm3D>);
static_assert(std::is_move_constructible_v<MPIComm3D>);

std::pair<int, int> MPIComm3D::shifted_ranks(int dim, int disp) const noexcept
{
  std::pair<int, int> r(MPI_PROC_NULL, MPI_PROC_NULL);
  MPI_Cart_shift(comm.raw_handle(), 3 - dim - 1, disp, &r.first, &r.second);
  return r;
}

} // namespace alps
