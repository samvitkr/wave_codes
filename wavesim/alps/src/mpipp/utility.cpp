//
// Created by xuanx004 on 3/14/24.
//

#include "utility.h"

#include <numeric>

#include <fmt/format.h>
#include <mpipp/comm.h>

namespace mpipp {

std::vector<int> generate_displacements(nonstd::span<const int> counts)
{
  std::vector<int> displacements(counts.size(), 0);
  std::partial_sum(counts.begin(), counts.end() - 1, displacements.begin() + 1);
  return displacements;
}

MPI_Aint size_t_to_mpi_aint(std::size_t size)
{
  if (size > static_cast<std::size_t>(std::numeric_limits<MPI_Aint>::max())) {
    throw std::invalid_argument(
      fmt::format("size {} exceeds MPI_Aint range", size));
  }
  return static_cast<MPI_Aint>(size);
}

namespace detail {
void check_root(int root, const communicator& comm)
{
  if (root < 0 || root >= comm.size()) {
    throw std::invalid_argument("mpi root is invalid");
  }
  if (comm.rank() != root) {
    throw std::invalid_argument("mpi root is not the calling rank");
  }
}

void check_non_root(int root, const communicator& comm)
{
  if (root < 0 || root >= comm.size()) {
    throw std::invalid_argument("mpi root is invalid");
  }
  if (comm.rank() == root) {
    throw std::invalid_argument("mpi root is the calling rank");
  }
}

void check_gather_size(std::size_t         send_count,
                       std::size_t         recv_size,
                       const communicator& comm)
{
  if (send_count * comm.size() > recv_size) {
    throw std::invalid_argument(
      "mpi gather receive size is less than required");
  }
}

} // namespace detail
} // namespace mpipp
