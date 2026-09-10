#pragma once

#include <Kokkos_Core.hpp>

namespace alps {

template<int Rank, typename... Properties>
using LoopPolicy = std::conditional_t<
  Rank == 1,
  Kokkos::RangePolicy<Kokkos::IndexType<int>, Properties...>,
  Kokkos::MDRangePolicy<
    Kokkos::Rank<Rank, Kokkos::Iterate::Left, Kokkos::Iterate::Left>,
    Kokkos::IndexType<int>,
    Properties...>>;

template<typename... Properties>
using GridPolicy = Kokkos::TeamPolicy<Kokkos::IndexType<int>, Properties...>;

} // namespace alps
