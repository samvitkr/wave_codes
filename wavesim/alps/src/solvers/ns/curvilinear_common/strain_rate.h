#pragma once

#include <common/container/matrix_field.h>
#include <common/container/vector_field.h>
#include <common/container/view_utils.h>
#include <common/device/device_traits.h>
#include <common/kokkos_abstraction/exec_policy.h>
#include <common/kokkos_abstraction/pool_space.h>
#include <common/real_type.h>
#include <common/runtime/async_utils.h>
#include <solvers/operators/grad_curvilinear.h>

namespace alps::solver {

template<typename MT>
void calculate_strain_rate(SymmTensor33Field<Real***> const& Sij,
                           Vector3Field<Real***> const&      u_vec,
                           MT const&                         mesh)
{
  using Kokkos::parallel_for;

  const auto [nx, ny, nz] = mesh.extents();

  const auto     stream1   = get_next_stream();
  constexpr auto tile_size = []() -> Kokkos::Array<range_index_t, 3> {
    using Device = Kokkos::DefaultExecutionSpace;
    if constexpr (is_cuda_execution_space_v<Device>) return {16, 1, 4};
    if constexpr (is_hip_execution_space_v<Device>) return {64, 1, 4};
    return {0, 0, 0};
  }();
  LoopPolicy<3> const policy1(stream1, {0, 0, 0}, {nx, ny, nz}, tile_size);

  // du/dx, du/dy, du/dz
  grad({Sij.xx, Sij.xy, Sij.xz}, u_vec.x, mesh, CenterPt());

  // calculate dv/dx, dv/dy, dv/dz
  // use Sij.zz as temporary storage for dv/dx
  grad({Sij.zz, Sij.yy, Sij.yz}, u_vec.y, mesh, CenterPt());

  const auto& vx = Sij.zz;
  const auto& uy = Sij.xy;
  parallel_for(
    "S12", policy1, KOKKOS_LAMBDA(int i, int j, int k) {
      uy(i, j, k) = (uy(i, j, k) + vx(i, j, k)) / 2;
    });
  policy1.space().fence();

  // dw/dx and dw/dy need temporay storage
  HaloView<Real***, default_memory_pool> wx(
    Kokkos::view_alloc("wx", Kokkos::WithoutInitializing),
    Sij.xx.layout(),
    {Sij.xx.begin(0), Sij.xx.begin(1), Sij.xx.begin(2)});
  HaloView<Real***, default_memory_pool> wy(
    Kokkos::view_alloc("wy", Kokkos::WithoutInitializing),
    Sij.xx.layout(),
    {Sij.xx.begin(0), Sij.xx.begin(1), Sij.xx.begin(2)});

  // calculate dw/dx, dw/dy, dw/dz
  grad({wx, wy, Sij.zz}, u_vec.z, mesh, NodePt());

  const auto& uz = Sij.xz;
  const auto& vz = Sij.yz;
  parallel_for(
    "S13 and S23", policy1, KOKKOS_LAMBDA(int i, int j, int k) {
      uz(i, j, k) = (uz(i, j, k) + wx(i, j, k)) / 2;
      vz(i, j, k) = (vz(i, j, k) + wy(i, j, k)) / 2;
    });
  policy1.space().fence();
}

} // namespace alps::solver
