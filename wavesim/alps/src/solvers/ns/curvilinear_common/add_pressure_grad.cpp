//
// Created by xuanx004 on 7/16/24.
//

#include "add_pressure_grad.h"

#include <common/container/view_utils.h>
#include <common/device/device_traits.h>
#include <common/kokkos_abstraction/exec_policy.h>
#include <common/runtime/async_utils.h>

namespace alps::solver {
void add_pressure_gradient(const Vector3Field<Real***>& source,
                           const CurvilinearMesh&       mesh,
                           std::array<Real, 2> const    gradients)
{
  std::vector<decltype(get_next_stream())> streams;

  auto const dpx = gradients[0];
  auto const dpy = gradients[1];

  if (std::abs(dpx) < std::numeric_limits<Real>::epsilon()
      && std::abs(dpy) < std::numeric_limits<Real>::epsilon())
    return;

  auto const& invJ = mesh.invJ;
  auto const& fx   = source.x;
  auto const& fy   = source.y;

  auto constexpr tile_size = []() -> Kokkos::Array<long, 3> {
    using Device = Kokkos::DefaultExecutionSpace;
    if constexpr (is_cuda_execution_space_v<Device>) return {16, 4, 2};
    if constexpr (is_hip_execution_space_v<Device>) return {64, 1, 4};
    return {0, 0, 0};
  }();

  auto const stream = get_next_stream();
  auto const policy =
    LoopPolicy<3>(stream, local_begins(fx), local_ends(fx), tile_size);

  if (std::abs(dpx) < std::numeric_limits<Real>::epsilon()) {
    Kokkos::parallel_for(
      "add dp/dx", policy, KOKKOS_LAMBDA(int i, int j, int k) {
        fy(i, j, k) += dpy * invJ(i, j);
      });
  } else if (std::abs(dpy) < std::numeric_limits<Real>::epsilon()) {
    Kokkos::parallel_for(
      "add dp/dx", policy, KOKKOS_LAMBDA(int i, int j, int k) {
        fx(i, j, k) += dpx * invJ(i, j);
      });
  } else {
    Kokkos::parallel_for(
      "add dp/dx", policy, KOKKOS_LAMBDA(int i, int j, int k) {
        fx(i, j, k) += dpx * invJ(i, j);
        fy(i, j, k) += dpy * invJ(i, j);
      });
  }
  stream.fence();
}
} // namespace alps::solver
