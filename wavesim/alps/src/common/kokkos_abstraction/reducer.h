#pragma once

#include <Kokkos_Array.hpp>

namespace Kokkos {

template<class T>
struct reduction_identity;

/// Specialization of reduction_identity for Kokkos::Array type to be able to
/// use as loc index with MaxLoc or MinLoc
template<>
struct reduction_identity<Kokkos::Array<int, 3>>
{
  KOKKOS_FORCEINLINE_FUNCTION static Kokkos::Array<int, 3> min()
  {
    return {0, 0, 0};
  }
};

} // namespace Kokkos
