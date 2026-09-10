//
// Created by xuanx004 on 10/4/22.
//

#pragma once

#include <common/container/vector_field.h>
#include <common/container/view_utils.h>
#include <common/kokkos_abstraction/exec_policy.h>
#include <common/kokkos_abstraction/pool_space.h>
#include <common/real_type.h>
#include <common/runtime/async_utils.h>
#include <decomp/ghost_cell_exchange.h>
#include <solvers/mesh/curvilinear_mesh.h>
#include <solvers/operators/grad.h>
#include <spectral/spectral.h>

#include <Kokkos_Profiling_ScopedRegion.hpp>

namespace alps::solver::detail {

#define M_INVJ_G11     (invJ.data()[offset_ij])
#define M_INVJ_G22     (invJ.data()[offset_ij])
#define M_INVJ_G13(_z) (MT::invJ_zeta_x(_z, eta_x.data()[offset_ij]))
#define M_INVJ_G23(_z) (MT::invJ_zeta_y(_z, eta_y.data()[offset_ij]))
#define M_INVJ_G33(_z)                            \
  (MT::invJ_zeta_x(_z, eta_x.data()[offset_ij])   \
     * MT::zeta_x(_z, exr.data()[offset_ij])      \
   + MT::invJ_zeta_y(_z, eta_y.data()[offset_ij]) \
       * MT::zeta_y(_z, eyr.data()[offset_ij])    \
   + MT::invJ_zeta_z() * J.data()[offset_ij])

template<typename FieldType>
void add_SGS_and_molecular_diffusion_fluxes_impl(
  Vector3Field<Real***> const&   fluxes,
  const HaloView<Real const***>& f,
  const FieldType&               flow,
  HaloView<Real***> const&       nu_D,
  Real const                     D)
{
  using boundary_functor_t =
    FDBoundaryFunctor<Real, Kokkos::DefaultExecutionSpace>;

  auto const region =
    Kokkos::Profiling::ScopedRegion("SGS diffusion fluxes " + f.label());

  using MT                = std::decay_t<decltype(flow.mesh)>;
  const auto& mesh        = flow.mesh;
  const auto& grid        = mesh.grid;
  const auto [nx, ny, nz] = local_extents(f);
  const auto is_top       = grid.comm().is_last(2);
  const auto is_bottom    = grid.comm().is_first(2);

  const auto& zz    = mesh.zz;
  const auto& zw    = mesh.zw;
  const auto& dz    = mesh.dz;
  const auto& dzw   = mesh.dzw;
  const auto& invJ  = mesh.invJ;
  const auto& J     = mesh.J;
  const auto& eta_x = mesh.ex;
  const auto& eta_y = mesh.ey;
  const auto& exr   = mesh.exr;
  const auto& eyr   = mesh.eyr;

  HaloView<Real****, default_memory_pool> const tmp(
    Kokkos::view_alloc("viscous_tmp_derivatives", Kokkos::WithoutInitializing),
    Kokkos::LayoutLeft(f.extent(0), f.extent(1), f.extent(2), 3),
    {begin(f, 0), begin(f, 1), begin(f, 2), 0});
  const auto& tmp_x       = subview(tmp, ALL, ALL, ALL, 0);
  const auto& tmp_y       = subview(tmp, ALL, ALL, ALL, 1);
  const auto& tmp_z       = subview(tmp, ALL, ALL, ALL, 2);
  auto        tmp_x_inner = create_inner_view(tmp_x).view();
  auto        tmp_y_inner = create_inner_view(tmp_y).view();

  const auto     stream1   = get_next_stream();
  const auto     stream2   = get_next_stream();
  constexpr auto tile_size = []() -> Kokkos::Array<range_index_t, 3> {
    using Device = Kokkos::DefaultExecutionSpace;
    if constexpr (is_cuda_execution_space_v<Device>) return {16, 1, 8};
    if constexpr (is_hip_execution_space_v<Device>) return {64, 1, 4};
    return {0, 0, 0};
  }();

  auto req = async_update_halo_upper_z(grid.partition(), nu_D, 1);

  auto inner = create_inner_view(f).view();

  /*#region Calculate df/dξ, df/d𝜓, df/dζ */
  spectral::ddx(tmp_x_inner, inner, grid, stream1); // df/dξ
  spectral::ddy(tmp_y_inner, inner, grid, stream2); // df/d𝜓

  // Ghost cells update must be submitted after the calculations complete
  stream1.fence();
  auto reqs1 = async_update_halo_upper_z(grid.x_pencil(), tmp_x, 2);
  stream2.fence();
  auto reqs2 = async_update_halo_upper_z(grid.x_pencil(), tmp_y, 3);

  auto event = get_device_event();
  if (is_bottom) {
    auto const               policy = LoopPolicy<2>(stream2, {0, 0}, {nx, ny});
    boundary_functor_t const functor(
      subview(tmp_z, ALL, ALL, 0).view(),
      subview(f, ALL, ALL, index_range(0, 3)).view(),
      subview(mesh.dz_h, index_range(0, 2)).view(),
      1,
      LeftBoundary());
    parallel_for("du/dzeta bottom", policy, functor);
  }
  if (is_top) {
    auto const               policy = LoopPolicy<2>(stream2, {0, 0}, {nx, ny});
    boundary_functor_t const functor(
      subview(tmp_z, ALL, ALL, nz - 2).view(),
      subview(f, ALL, ALL, index_range(nz - 3, nz)).view(),
      subview(mesh.dz_h, index_range(nz - 3, nz - 1)).view(),
      1,
      RightBoundary());
    parallel_for("du/dzeta top", policy, functor);
  }
  if (is_bottom || is_top) enqueue(event, stream2);

  // df/dζ
  auto policy = LoopPolicy<3>(stream1,
                              {0, 0, is_bottom ? 1 : -1},
                              {nx, ny, is_top ? nz - 2 : nz},
                              tile_size);
  parallel_for(
    "du/dzeta", policy, KOKKOS_LAMBDA(int i, int j, int k) {
      auto alpha     = 1 / dz(k);
      tmp_z(i, j, k) = alpha * (f(i, j, k + 1) - f(i, j, k));
    });

  // Progress MPI ops while the above calculation executed
  reqs1.waitall();
  reqs2.waitall();
  if (is_bottom || is_top) wait_for(event, stream1);
  /*#endregion*/

  /*#region Calculate F_j = J^{-1} g^{ij} df/dξ_i */
  const auto& fx = fluxes.x;
  const auto& fy = fluxes.y;
  const auto& fz = fluxes.z;

  if (!check_same_layout_and_offset(fx, fy, fz, tmp_x, tmp_y, tmp_z, nu_D)) {
    throw std::runtime_error("Mismatched extents in f_{xyz}, tmp_{xyz}, nu_D");
  }

  req.waitall();
  const auto functor = KOKKOS_LAMBDA(int i, int j, int k)
  {
    auto offset    = &tmp_x(i, j, k) - tmp_x.data();
    auto stride    = tmp_x.stride(2);
    auto offset_ij = &invJ(i, j) - invJ.data();

    auto u_zeta =
      itp2center(tmp_z.data()[offset - stride], tmp_z.data()[offset]);
    auto coeff = nu_D.data()[offset] + D;
    fx.data()[offset] +=
      (M_INVJ_G11 * tmp_x.data()[offset] + M_INVJ_G13(zz(k)) * u_zeta) * coeff;
    fy.data()[offset] +=
      (M_INVJ_G22 * tmp_y.data()[offset] + M_INVJ_G23(zz(k)) * u_zeta) * coeff;
    auto ratio = dzw(k - 1) / (dzw(k - 1) + dzw(k));
    auto u_xi =
      itp2node(tmp_x.data()[offset], tmp_x.data()[offset + stride], ratio);
    auto u_psi =
      itp2node(tmp_y.data()[offset], tmp_y.data()[offset + stride], ratio);
    auto nud_w =
      itp2node(nu_D.data()[offset], nu_D.data()[offset + stride], ratio);
    fz.data()[offset] += (M_INVJ_G13(zw(k)) * u_xi + M_INVJ_G23(zw(k)) * u_psi
                          + M_INVJ_G33(zw(k)) * tmp_z.data()[offset])
                       * (nud_w + D);
  };
  parallel_for("viscous flux u",
               LoopPolicy<3>(
                 stream1, {0, 0, 0}, {nx, ny, is_top ? nz - 1 : nz}, tile_size),
               functor);
  stream1.fence(); // Make sure tmp_x, tmp_y not conflict with the next loop
  /*#endregion*/

  stream1.fence();
  stream2.fence();
}

#undef M_INVJ_G11
#undef M_INVJ_G22
#undef M_INVJ_G13
#undef M_INVJ_G23
#undef M_INVJ_G33

template<typename FlowType>
void add_sgs_diffusion_fluxes(Vector3Field<Real***> const& fluxes,
                              Vector3Field<Real***> const& grad_f,
                              FlowType const&              flow,
                              HaloView<Real***> const&     nuD,
                              Real const                   D)
{
  using Kokkos::parallel_for;

  if (nuD.begin(2) > -1 || grad_f.x.begin(2) > -1 || grad_f.y.begin(2) > -1) {
    throw std::runtime_error(
      "nuD and grad(f) must have at least one ghost cell in z");
  }

  auto const region =
    Kokkos::Profiling::ScopedRegion("SGS diffusion fluxes " + grad_f.x.label());

  using MT          = std::decay_t<decltype(flow.mesh)>;
  const auto& mesh  = flow.mesh;
  const auto& invJ  = mesh.invJ;
  const auto& eta_x = mesh.ex;
  const auto& eta_y = mesh.ey;
  const auto& zw    = mesh.zw;
  const auto& dzw   = mesh.dzw;
  const auto& dz    = mesh.dz;

  const auto [nx, ny, nz] = mesh.extents();

  const auto     stream1   = get_next_stream();
  constexpr auto tile_size = []() -> Kokkos::Array<range_index_t, 3> {
    using Device = Kokkos::DefaultExecutionSpace;
    if constexpr (is_cuda_execution_space_v<Device>) return {16, 1, 4};
    if constexpr (is_hip_execution_space_v<Device>) return {64, 1, 4};
    return {0, 0, 0};
  }();

  auto reqs1 = async_update_halo_upper_z(mesh.partition(), grad_f.x, 1);
  reqs1.push(async_update_halo_upper_z(mesh.partition(), grad_f.y, 2));
  reqs1.push(async_update_halo_upper_z(mesh.partition(), nuD, 3));

  if (!check_same_layout_and_offset(
        fluxes.x, fluxes.y, fluxes.z, grad_f.x, grad_f.y, grad_f.z, nuD)) {
    throw std::runtime_error("Mismatched extents in fluxes, grad_f, nuD");
  }

  const LoopPolicy<3> policy(stream1, {0, 0, 0}, {nx, ny, nz}, tile_size);
  reqs1.waitall();
  parallel_for(
    "nuD*grad(f)", policy, KOKKOS_LAMBDA(int i, int j, int k) {
      auto offset    = &grad_f.x(i, j, k) - grad_f.x.data();
      auto stride    = grad_f.x.stride(2);
      auto offset_ij = &invJ(i, j) - invJ.data();

      const auto fx_k = grad_f.x.data()[offset];
      const auto fy_k = grad_f.y.data()[offset];
      fluxes.x.data()[offset] +=
        (nuD.data()[offset] + D) * fx_k * invJ.data()[offset_ij];
      fluxes.y.data()[offset] +=
        (nuD.data()[offset] + D) * fy_k * invJ.data()[offset_ij];

      auto       ratio = dzw(k - 1) / 2 / dz(k);
      const auto fx_node =
        itp2node(fx_k, grad_f.x.data()[offset + stride], ratio);
      const auto fy_node =
        itp2node(fy_k, grad_f.y.data()[offset + stride], ratio);
      const auto fz = grad_f.z.data()[offset];
      const auto nuD_node =
        itp2node(nuD.data()[offset], nuD.data()[offset + stride], ratio);
      const auto qzeta =
        fx_node * MT::invJ_zeta_x(zw(k), eta_x.data()[offset_ij])
        + fy_node * MT::invJ_zeta_y(zw(k), eta_y.data()[offset_ij])
        + fz * MT::invJ_zeta_z();
      fluxes.z.data()[offset] += (nuD_node + D) * qzeta;
    });

  stream1.fence();
}

} // namespace alps::solver::detail
