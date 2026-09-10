#pragma once

#include <common/container/view_types.h>
#include <common/kokkos_abstraction/exec_policy.h>

namespace alps {

/// Calculate the maximum absolute value of a 3D array.
template<class T>
[[nodiscard]] extern T nrminf_3d(MDView<T const***> const&            x,
                                 Kokkos::DefaultExecutionSpace const& stream);

/// Calculate the L-inf distance between two 3D arrays.
template<class T>
[[nodiscard]] extern T nrminf_distance_3d(MDView<T const***> const& x,
                                          MDView<T const***> const& y);

/// Calculate the L-inf distance between two 3D arrays given a parallel policy.
template<class T>
[[nodiscard]] extern T nrminf_distance_3d(MDView<T const***> const& x,
                                          MDView<T const***> const& y,
                                          LoopPolicy<3> const&      policy);

} // namespace alps
