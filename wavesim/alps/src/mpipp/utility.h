//
// Created by xuanx004 on 3/14/24.
//

#pragma once

#include <limits>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include <mpi.h>
#include <nonstd/span.hpp>

namespace mpipp {

class communicator;

template<typename T, std::enable_if_t<std::is_signed_v<T>, bool> = true>
constexpr int to_int_size(T size)
{
  if (size >= 0 && size <= std::numeric_limits<int>::max()) {
    return static_cast<int>(size);
  }
  throw std::invalid_argument("mpi to_int_size: size is too large");
}

template<typename T, std::enable_if_t<std::is_unsigned_v<T>, bool> = true>
constexpr int to_int_size(T size)
{
  if (size <= std::numeric_limits<int>::max()) {
    return static_cast<int>(size);
  }
  throw std::invalid_argument("mpi to_int_size: size is too large");
}

std::vector<int> generate_displacements(nonstd::span<const int> counts);

MPI_Aint size_t_to_mpi_aint(std::size_t size);

namespace detail {
[[maybe_unused]] void check_root(int root, const communicator& comm);

[[maybe_unused]] void check_non_root(int root, const communicator& comm);

[[maybe_unused]] void check_gather_size(std::size_t         send_count,
                                        std::size_t         recv_size,
                                        const communicator& comm);

} // namespace detail

} // namespace mpipp
