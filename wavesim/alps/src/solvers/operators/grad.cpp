#include "grad.h"

#include <common/container/view_types.h>
#include <common/container/view_utils.h>
#include <common/kokkos_abstraction/exec_policy.h>
#include <common/kokkos_abstraction/pool_space.h>
#include <common/runtime/async_utils.h>
#include <spectral/spectral.h>

#include <Kokkos_Profiling_ScopedRegion.hpp>

namespace alps {
namespace solver {

using Kokkos::parallel_for;

template<class T, class ExecSpace>
FDBoundaryFunctor<T, ExecSpace>::FDBoundaryFunctor(
  MDView<T**, ExecSpace>                output,
  MDView<T const** [3], ExecSpace>      input_stencil,
  MDView<T const[2], Kokkos::HostSpace> delta,
  T                                     Lz,
  LeftBoundary /*tag*/)
  : fz{output}
  , f{input_stencil}
{
  const auto dz0   = delta(0);
  const auto dz1   = delta(1);
  const auto alpha = dz0 / (dz0 + dz1);
  coeff[0]         = -(1 + alpha) / dz0 / Lz;
  coeff[1]         = 1 / alpha / dz1 / Lz;
  coeff[2]         = -alpha / dz1 / Lz;
}

template<class T, class ExecSpace>
FDBoundaryFunctor<T, ExecSpace>::FDBoundaryFunctor(
  MDView<T**, ExecSpace>                output,
  MDView<T const** [3], ExecSpace>      input_stencil,
  MDView<T const[2], Kokkos::HostSpace> delta,
  T                                     Lz,
  RightBoundary /*tag*/)
  : fz{output}
  , f{input_stencil}
{
  const auto dz0   = delta(1);
  const auto dz1   = delta(0);
  const auto alpha = dz0 / (dz0 + dz1);
  coeff[2]         = (1 + alpha) / dz0 / Lz;
  coeff[1]         = -1 / alpha / dz1 / Lz;
  coeff[0]         = alpha / dz1 / Lz;
}

template struct FDBoundaryFunctor<Real, Kokkos::DefaultExecutionSpace>;

void grad(const Vector3Field<Real***>& grad_u,
          const HaloView<Real***>&     u,
          const Mesh&                  mesh,
          CenterPt /*tag*/)
{
  auto const region =
    Kokkos::Profiling::ScopedRegion("grad(" + u.label() + ")");

  using Kokkos::parallel_for;

  const auto& grid      = mesh.grid;
  const auto  hbar      = mesh.hbar;
  const auto& dz        = mesh.dz;
  const auto  is_top    = grid.comm().is_last(2);
  const auto  is_bottom = grid.comm().is_first(2);

  auto nx = local_extent(u, 0);
  auto ny = local_extent(u, 1);
  auto nz = local_extent(u, 2);

  const auto stream1 = get_next_stream();
  const auto stream2 = get_next_stream();

  auto u_inner  = create_inner_view(u).view();
  auto ux_inner = create_inner_view(grad_u.x).view();
  auto uy_inner = create_inner_view(grad_u.y).view();
  spectral::ddx(ux_inner, u_inner, grid, stream1);
  spectral::ddy(uy_inner, u_inner, grid, stream2);

  const auto& uz = grad_u.z;

  auto z_begin = is_bottom ? 1 : 0;
  auto z_end   = is_top ? nz - 2 : nz;
  if (is_bottom) {
    FDBoundaryFunctor<Real, std::decay_t<decltype(stream2)>> const functor(
      subview(uz, ALL, ALL, 0).view(),
      subview(u, ALL, ALL, index_range(0, 3)).view(),
      subview(mesh.dz_h, index_range(0, 2)).view(),
      mesh.hbar,
      LeftBoundary());
    parallel_for(LoopPolicy<2>(stream2, {0, 0}, {nx, ny}), functor);
  }
  if (is_top) {
    FDBoundaryFunctor<Real, std::decay_t<decltype(stream2)>> const functor(
      subview(uz, ALL, ALL, nz - 2).view(),
      subview(u, ALL, ALL, index_range(nz - 3, nz)).view(),
      subview(mesh.dz_h, index_range(nz - 3, nz - 1)).view(),
      mesh.hbar,
      RightBoundary());
    parallel_for(LoopPolicy<2>(stream2, {0, 0}, {nx, ny}), functor);
  }

  constexpr auto tile_size = []() -> Kokkos::Array<range_index_t, 3> {
    using Device = Kokkos::DefaultExecutionSpace;
    if constexpr (is_cuda_execution_space_v<Device>) return {16, 1, 8};
    if constexpr (is_hip_execution_space_v<Device>) return {64, 1, 4};
    return {0, 0, 0};
  }();
  LoopPolicy<3> const policy1(
    stream1, {0, 0, z_begin}, {nx, ny, z_end}, tile_size);
  parallel_for(
    "grad_ddz", policy1, KOKKOS_LAMBDA(int i, int j, int k) {
      auto alpha  = 1 / (dz(k) * hbar);
      uz(i, j, k) = alpha * (u(i, j, k + 1) - u(i, j, k));
    });

  stream1.fence();
  stream2.fence();
}

void grad(const Vector3Field<Real***>& grad_u,
          const HaloView<Real***>&     w,
          const Mesh&                  mesh,
          NodePt /*tag*/)
{
  auto const region =
    Kokkos::Profiling::ScopedRegion("grad(" + w.label() + ")");

  using Kokkos::parallel_for;

  const auto& grid      = mesh.grid;
  auto        hbar      = mesh.hbar;
  const auto& dzw       = mesh.dzw;
  auto        is_top    = grid.comm().is_last(2);
  auto        is_bottom = grid.comm().is_first(2);

  auto nx = local_extent(w, 0);
  auto ny = local_extent(w, 1);
  auto nz = local_extent(w, 2);

  const auto stream1 = get_next_stream();
  const auto stream2 = get_next_stream();

  auto w_inner  = create_inner_view(w).view();
  auto wx_inner = create_inner_view(grad_u.x).view();
  auto wy_inner = create_inner_view(grad_u.y).view();
  spectral::ddx(wx_inner, w_inner, grid, stream1);
  spectral::ddy(wy_inner, w_inner, grid, stream2);

  const auto& wz = grad_u.z;

  auto z_begin = is_bottom ? 1 : 0;
  auto z_end   = is_top ? nz - 1 : nz;
  if (is_bottom) {
    FDBoundaryFunctor<Real, std::decay_t<decltype(stream2)>> const functor(
      subview(wz, ALL, ALL, 0).view(),
      subview(w_inner, ALL, ALL, index_range(0, 3)),
      subview(mesh.dzw_h, index_range(0, 2)).view(),
      mesh.hbar,
      LeftBoundary());
    parallel_for(LoopPolicy<2>(stream2, {0, 0}, {nx, ny}), functor);
  }
  if (is_top) {
    FDBoundaryFunctor<Real, std::decay_t<decltype(stream2)>> const functor(
      subview(wz, ALL, ALL, nz - 1).view(),
      subview(w_inner, ALL, ALL, index_range(nz - 4, nz - 1)),
      subview(mesh.dzw_h, index_range(nz - 4, nz - 2)).view(),
      mesh.hbar,
      RightBoundary());
    parallel_for(LoopPolicy<2>(stream2, {0, 0}, {nx, ny}), functor);
  }

  constexpr auto tile_size = []() -> Kokkos::Array<range_index_t, 3> {
    using Device = Kokkos::DefaultExecutionSpace;
    if constexpr (is_cuda_execution_space_v<Device>) return {16, 1, 8};
    if constexpr (is_hip_execution_space_v<Device>) return {64, 1, 4};
    return {0, 0, 0};
  }();
  LoopPolicy<3> const policy(
    stream1, {0, 0, z_begin}, {nx, ny, z_end}, tile_size);
  parallel_for(
    "grad_ddz", policy, KOKKOS_LAMBDA(int i, int j, int k) {
      auto beta   = 1 / (dzw(k - 1) * hbar);
      wz(i, j, k) = beta * (w(i, j, k) - w(i, j, k - 1));
    });

  stream1.fence();
  stream2.fence();
}

Vector3Field<Real***>
grad(const HaloView<Real***>& u, const Mesh& mesh, CenterPt /*unused*/)
{
  Vector3Field<Real***, default_memory_pool> grad_u(
    Kokkos::view_alloc("grad(" + u.label() + ")", Kokkos::WithoutInitializing),
    local_extents(u),
    {0, 0, 1});

  grad(grad_u, u, mesh, CenterPt());

  return {grad_u.x, grad_u.y, grad_u.z};
}

Vector3Field<Real***>
grad(const HaloView<Real***>& w, const Mesh& mesh, NodePt /*unused*/)
{
  Vector3Field<Real***, default_memory_pool> grad_u(
    Kokkos::view_alloc("grad(" + w.label() + ")", Kokkos::WithoutInitializing),
    local_extents(w),
    {0, 0, 1});

  grad(grad_u, w, mesh, NodePt());

  return {grad_u.x, grad_u.y, grad_u.z};
}
} // namespace solver
} // namespace alps
