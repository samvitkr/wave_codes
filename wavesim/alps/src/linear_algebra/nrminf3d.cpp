#include "nrminf3d.h"

#include <common/container/view_types.h>
#include <common/container/view_utils.h>
#include <common/kokkos_abstraction/exec_policy.h>

#include <Kokkos_Core.hpp>

namespace alps {
namespace detail {
template<class T>
struct nrminf_3d_kernel
{
  nrminf_3d_kernel(MDView<T const***> const& x)
    : x_{x}
  {}

  KOKKOS_INLINE_FUNCTION void operator()(int i, int j, int k, T& max) const
  {
    auto val = Kokkos::abs(x_(i, j, k));
    if (val > max) max = val;
  }

  MDView<T const***> x_;
};
} // namespace detail

template<class T>
[[nodiscard]] T nrminf_3d(MDView<T const***> const&            x,
                          Kokkos::DefaultExecutionSpace const& stream)
{
  T             r;
  LoopPolicy<3> policy(stream, begins(x), ends(x));
  Kokkos::parallel_reduce(
    "nrminf3d", policy, detail::nrminf_3d_kernel<T>(x), Kokkos::Max<T>(r));

  policy.space().fence();
  return r;
}

template float  nrminf_3d(MDView<float const***> const&,
                          Kokkos::DefaultExecutionSpace const&);
template double nrminf_3d(MDView<double const***> const&,
                          Kokkos::DefaultExecutionSpace const&);

template<class T>
[[nodiscard]] T nrminf_distance_3d(MDView<T const***> const& x,
                                   MDView<T const***> const& y,
                                   LoopPolicy<3> const&      policy)
{
  T r;
  Kokkos::parallel_reduce(
    "nrminf3d",
    policy,
    KOKKOS_LAMBDA(int i, int j, int k, T& max) {
      auto val = Kokkos::abs(x(i, j, k) - y(i, j, k));
      if (val > max) max = val;
    },
    Kokkos::Max<T>(r));

  policy.space().fence();
  return r;
}

template<class T>
[[nodiscard]] T nrminf_distance_3d(MDView<T const***> const& x,
                                   MDView<T const***> const& y)
{
  return nrminf_distance_3d(x, y, LoopPolicy<3>(begins(x), ends(x)));
}

template float  nrminf_distance_3d(MDView<float const***> const&,
                                   MDView<float const***> const&);
template double nrminf_distance_3d(MDView<double const***> const&,
                                   MDView<double const***> const&);
template float  nrminf_distance_3d(MDView<float const***> const&,
                                   MDView<float const***> const&,
                                   LoopPolicy<3> const&);
template double nrminf_distance_3d(MDView<double const***> const&,
                                   MDView<double const***> const&,
                                   LoopPolicy<3> const&);
} // namespace alps
