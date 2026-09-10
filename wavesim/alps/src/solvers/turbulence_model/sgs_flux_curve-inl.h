//
// Created by xuanx004 on 7/11/24.
//

#pragma once

#include <common/container/matrix_field.h>
#include <common/container/view_types.h>
#include <common/container/view_utils.h>
#include <common/kokkos_abstraction/exec_policy.h>
#include <common/real_type.h>
#include <common/runtime/async_utils.h>
#include <decomp/ghost_cell_exchange.h>
#include <solvers/mesh/curvilinear_mesh.h>
#include <spectral/spectral.h>

#include <Kokkos_Profiling_ScopedRegion.hpp>

namespace alps::solver::detail {

template<typename FlowType>
void add_sgs_momentum_fluxes_from_Sij(Tensor33Field<Real***> const&     fluxes,
                                      FlowType const&                   flow,
                                      HaloView<Real***> const&          nu_t,
                                      SymmTensor33Field<Real***> const& Sij)
{
  using Kokkos::parallel_for;
  using std::tie;

  auto const region = Kokkos::Profiling::ScopedRegion("SGS momentum fluxes");

  using MT          = std::decay_t<decltype(flow.mesh)>;
  const auto& mesh  = flow.mesh;
  const auto& grid  = mesh.grid;
  const auto& invJ  = mesh.invJ;
  const auto& eta_x = mesh.ex;
  const auto& eta_y = mesh.ey;
  const auto& zz    = mesh.zz;
  const auto& zw    = mesh.zw;
  const auto& dzw   = mesh.dzw;
  const auto& dz    = mesh.dz;

  const auto [nx, ny, nz] = grid.extents();

  const auto     stream1   = get_next_stream();
  constexpr auto tile_size = []() -> Kokkos::Array<range_index_t, 3> {
    using Device = Kokkos::DefaultExecutionSpace;
    if constexpr (is_cuda_execution_space_v<Device>) return {16, 1, 8};
    if constexpr (is_hip_execution_space_v<Device>) return {64, 1, 2};
    return {0, 0, 0};
  }();

  update_halo_upper_z(mesh.grid.partition(), Sij.xx, 1);
  update_halo_upper_z(mesh.grid.partition(), Sij.xy, 2);
  update_halo_upper_z(mesh.grid.partition(), Sij.yy, 3);
  update_halo_upper_z(mesh.grid.partition(), nu_t, 4);
  auto req1 = async_update_halo_lower_z(mesh.grid.partition(), Sij.xz, 5);
  auto req2 = async_update_halo_lower_z(mesh.grid.partition(), Sij.yz, 6);

  // Add fluxes related to S11, S12, S22, S33
  auto const& flux_xx = create_inner_view(fluxes.xx).view();
  auto const& flux_xy = create_inner_view(fluxes.xy).view();
  auto const& flux_yx = create_inner_view(fluxes.yx).view();
  auto const& flux_yy = create_inner_view(fluxes.yy).view();
  auto const& flux_zz = create_inner_view(fluxes.zz).view();
  auto const& flux_xz = create_inner_view(fluxes.xz).view();
  auto const& flux_yz = create_inner_view(fluxes.yz).view();
  auto const& Sxx     = Sij.xx;
  auto const& Syy     = Sij.yy;
  auto const& Sxy     = Sij.xy;
  auto const& Szz     = Sij.zz;

  if (!check_same_layout_and_offset(Sxx, Syy, Sxy, Szz, nu_t)) {
    throw std::runtime_error("Mismatched extents in Sij and nu_t");
  }

  LoopPolicy<3> const policy(stream1, {0, 0, 0}, {nx, ny, nz}, tile_size);
  parallel_for(
    "nu*ux", policy, KOKKOS_LAMBDA(int i, int j, int k) {
      const auto offset    = &flux_xx(i, j, k) - flux_xx.data();
      const auto offset_hv = &Sxx(i, j, k) - Sxx.data();
      const auto offset_ij = &invJ(i, j) - invJ.data();
      const auto stride    = Sxx.stride(2);

      const auto Sxx_k    = Sxx.data()[offset_hv];
      const auto Syy_k    = Syy.data()[offset_hv];
      const auto Sxy_k    = Sxy.data()[offset_hv];
      const auto Szz_k    = Szz.data()[offset_hv];
      const auto twiceNuT = nu_t.data()[offset_hv] * 2;
      flux_xx.data()[offset] += twiceNuT * Sxx_k * invJ.data()[offset_ij];
      flux_xy.data()[offset] += twiceNuT * Sxy_k * invJ.data()[offset_ij];
      flux_yx.data()[offset] += twiceNuT * Sxy_k * invJ.data()[offset_ij];
      flux_yy.data()[offset] += twiceNuT * Syy_k * invJ.data()[offset_ij];
      flux_zz.data()[offset] += twiceNuT * Szz_k * MT::invJ_zeta_z();

      auto       ratio = dzw(k - 1) / 2 / dz(k);
      const auto Sxx_node =
        itp2node(Sxx_k, Sxx.data()[offset_hv + stride], ratio);
      const auto Syy_node =
        itp2node(Syy_k, Syy.data()[offset_hv + stride], ratio);
      const auto Sxy_node =
        itp2node(Sxy_k, Sxy.data()[offset_hv + stride], ratio);
      const auto twiceNuT_node =
        itp2node(twiceNuT, nu_t.data()[offset_hv + stride] * 2, ratio);
      flux_xz.data()[offset] +=
        twiceNuT_node
        * (Sxx_node * MT::invJ_zeta_x(zw(k), eta_x.data()[offset_ij])
           + Sxy_node * MT::invJ_zeta_y(zw(k), eta_y.data()[offset_ij]));
      flux_yz.data()[offset] +=
        twiceNuT_node
        * (Sxy_node * MT::invJ_zeta_x(zw(k), eta_x.data()[offset_ij])
           + Syy_node * MT::invJ_zeta_y(zw(k), eta_y.data()[offset_ij]));
    });

  // Add fluxes related to S13 and S23
  auto const& flux_zx = create_inner_view(fluxes.zx).view();
  auto const& flux_zy = create_inner_view(fluxes.zy).view();
  auto const& Sxz     = Sij.xz;
  auto const& Syz     = Sij.yz;
  req1.waitall();
  req2.waitall();
  parallel_for(
    "nu*uz", policy, KOKKOS_LAMBDA(int i, int j, int k) {
      auto offset    = &flux_xx(i, j, k) - flux_xx.data();
      auto offset_hv = &Sxx(i, j, k) - Sxx.data();
      auto stride    = Sxx.stride(2);
      auto offset_ij = &invJ(i, j) - invJ.data();

      const auto Sxz_k    = Sxz.data()[offset_hv];
      const auto Syz_k    = Syz.data()[offset_hv];
      auto       ratio    = dzw(k - 1) / 2 / dz(k);
      const auto nut_node = 2
                          * itp2node(nu_t.data()[offset_hv],
                                     nu_t.data()[offset_hv + stride],
                                     ratio);
      flux_zx.data()[offset] += nut_node * Sxz_k * invJ.data()[offset_ij];
      flux_zy.data()[offset] += nut_node * Syz_k * invJ.data()[offset_ij];
      flux_xz.data()[offset] += nut_node * Sxz_k * MT::invJ_zeta_z();
      flux_yz.data()[offset] += nut_node * Syz_k * MT::invJ_zeta_z();
      const auto Sxz_c = itp2center(Sxz_k, Sxz.data()[offset_hv - stride]);
      const auto Syz_c = itp2center(Syz_k, Syz.data()[offset_hv - stride]);
      flux_zz.data()[offset] +=
        nu_t.data()[offset_hv] * 2
        * (Sxz_c * MT::invJ_zeta_x(zz(k), eta_x.data()[offset_ij])
           + Syz_c * MT::invJ_zeta_y(zz(k), eta_y.data()[offset_ij]));
    });

  stream1.fence();
}

template<typename FlowType>
void add_sgs_momentum_fluxes_from_gradu(Tensor33Field<Real***> const& fluxes,
                                        FlowType const&               flow,
                                        HaloView<Real***> const&      nu_t,
                                        Vector3Field<Real***> const&  grad_u,
                                        Vector3Field<Real***> const&  grad_v,
                                        Vector3Field<Real***> const&  grad_w)
{
  using Kokkos::parallel_for;
  using std::tie;

  auto const region = Kokkos::Profiling::ScopedRegion("SGS momentum fluxes");

  using MT          = std::decay_t<decltype(flow.mesh)>;
  const auto& mesh  = flow.mesh;
  const auto& grid  = mesh.grid;
  const auto& invJ  = mesh.invJ;
  const auto& eta_x = mesh.ex;
  const auto& eta_y = mesh.ey;
  const auto& zz    = mesh.zz;
  const auto& zw    = mesh.zw;
  const auto& dzw   = mesh.dzw;
  const auto& dz    = mesh.dz;

  const auto [nx, ny, nz] = grid.extents();

  const auto     stream1   = get_next_stream();
  constexpr auto tile_size = []() -> Kokkos::Array<range_index_t, 3> {
    using Device = Kokkos::DefaultExecutionSpace;
    if constexpr (is_cuda_execution_space_v<Device>) return {16, 1, 4};
    if constexpr (is_hip_execution_space_v<Device>) return {64, 1, 2};
    return {0, 0, 0};
  }();

  auto reqs1 = async_update_halo_upper_z(mesh.grid.partition(), grad_u.x, 1);
  reqs1.push(async_update_halo_upper_z(mesh.grid.partition(), grad_u.y, 2));
  reqs1.push(async_update_halo_upper_z(mesh.grid.partition(), grad_v.x, 3));
  reqs1.push(async_update_halo_upper_z(mesh.grid.partition(), grad_v.y, 4));
  reqs1.push(async_update_halo_upper_z(mesh.grid.partition(), nu_t, 5));

  // Add fluxes related to S11, S12, S22, S33
  auto const& flux_xx = create_inner_view(fluxes.xx).view();
  auto const& flux_xy = create_inner_view(fluxes.xy).view();
  auto const& flux_yx = create_inner_view(fluxes.yx).view();
  auto const& flux_yy = create_inner_view(fluxes.yy).view();
  auto const& flux_zz = create_inner_view(fluxes.zz).view();
  auto const& flux_xz = create_inner_view(fluxes.xz).view();
  auto const& flux_yz = create_inner_view(fluxes.yz).view();
  auto const& ux      = grad_u.x;
  auto const& uy      = grad_u.y;
  auto const& vx      = grad_v.x;
  auto const& vy      = grad_v.y;
  auto const& wz      = grad_w.z;

  if (!check_same_layout_and_offset(ux, uy, vx, vy, wz, nu_t)) {
    throw std::runtime_error("Mismatched extents in grad_u and nu_t");
  }

  const LoopPolicy<3> policy(stream1, {0, 0, 0}, {nx, ny, nz}, tile_size);
  reqs1.waitall();
  parallel_for(
    "nu*ux", policy, KOKKOS_LAMBDA(int i, int j, int k) {
      auto offset     = &flux_xx(i, j, k) - flux_xx.data();
      auto offset_hv  = &ux(i, j, k) - ux.data();
      auto offset1_hv = &ux(i, j, k + 1) - ux.data();
      auto offset_ij  = &invJ(i, j) - invJ.data();

      const auto Sxx_k      = ux.data()[offset_hv];
      const auto Syy_k      = vy.data()[offset_hv];
      const auto twiceSxy_k = uy.data()[offset_hv] + vx.data()[offset_hv];
      const auto Szz_k      = wz.data()[offset_hv];
      flux_xx.data()[offset] +=
        nu_t.data()[offset_hv] * 2 * Sxx_k * invJ.data()[offset_ij];
      flux_xy.data()[offset] +=
        nu_t.data()[offset_hv] * twiceSxy_k * invJ.data()[offset_ij];
      flux_yx.data()[offset] +=
        nu_t.data()[offset_hv] * twiceSxy_k * invJ.data()[offset_ij];
      flux_yy.data()[offset] +=
        nu_t.data()[offset_hv] * 2 * Syy_k * invJ.data()[offset_ij];
      flux_zz.data()[offset] +=
        nu_t.data()[offset_hv] * 2 * Szz_k * MT::invJ_zeta_z();
      auto       ratio         = dzw(k - 1) / 2 / dz(k);
      const auto Sxx_node      = itp2node(Sxx_k, ux.data()[offset1_hv], ratio);
      const auto Syy_node      = itp2node(Syy_k, vy.data()[offset1_hv], ratio);
      const auto twiceSxy_node = itp2node(
        twiceSxy_k, uy.data()[offset1_hv] + vx.data()[offset1_hv], ratio);
      const auto nut_node =
        itp2node(nu_t.data()[offset_hv], nu_t.data()[offset1_hv], ratio);
      flux_xz.data()[offset] +=
        nut_node
        * (2 * Sxx_node * MT::invJ_zeta_x(zw(k), eta_x.data()[offset_ij])
           + twiceSxy_node * MT::invJ_zeta_y(zw(k), eta_y.data()[offset_ij]));
      flux_yz.data()[offset] +=
        nut_node
        * (twiceSxy_node * MT::invJ_zeta_x(zw(k), eta_x.data()[offset_ij])
           + 2 * Syy_node * MT::invJ_zeta_y(zw(k), eta_y.data()[offset_ij]));
    });

  // Add fluxes related to S13 and S23
  auto reqs2 = async_update_halo_lower_z(mesh.grid.partition(), grad_w.x, 1);
  reqs2.push(async_update_halo_lower_z(mesh.grid.partition(), grad_w.y, 2));
  reqs2.push(async_update_halo_lower_z(mesh.grid.partition(), grad_u.z, 3));
  reqs2.push(async_update_halo_lower_z(mesh.grid.partition(), grad_v.z, 4));

  auto const& flux_zx = create_inner_view(fluxes.zx).view();
  auto const& flux_zy = create_inner_view(fluxes.zy).view();
  auto const& uz      = grad_u.z;
  auto const& vz      = grad_v.z;
  auto const& wx      = grad_w.x;
  auto const& wy      = grad_w.y;
  reqs2.waitall();
  parallel_for(
    "nu*uz", policy, KOKKOS_LAMBDA(int i, int j, int k) {
      auto offset     = &flux_xx(i, j, k) - flux_xx.data();
      auto offset_hv  = &ux(i, j, k) - ux.data();
      auto offset1_hv = &ux(i, j, k - 1) - ux.data();
      auto offset_ij  = &invJ(i, j) - invJ.data();

      const auto twiceSxz_k = uz.data()[offset_hv] + wx.data()[offset_hv];
      const auto twiceSyz_k = vz.data()[offset_hv] + wy.data()[offset_hv];
      auto       ratio      = dzw(k - 1) / 2 / dz(k);
      const auto nut_node   = itp2node(
        nu_t.data()[offset_hv], nu_t.data()[offset_hv + nu_t.stride(2)], ratio);
      flux_zx.data()[offset] += nut_node * twiceSxz_k * invJ.data()[offset_ij];
      flux_zy.data()[offset] += nut_node * twiceSyz_k * invJ.data()[offset_ij];
      flux_xz.data()[offset] += nut_node * twiceSxz_k * MT::invJ_zeta_z();
      flux_yz.data()[offset] += nut_node * twiceSyz_k * MT::invJ_zeta_z();
      const auto twiceSxz_c =
        itp2center(twiceSxz_k, uz.data()[offset1_hv] + wx.data()[offset1_hv]);
      const auto twiceSyz_c =
        itp2center(twiceSyz_k, vz.data()[offset1_hv] + wy.data()[offset1_hv]);
      flux_zz.data()[offset] +=
        nu_t.data()[offset_hv]
        * (twiceSxz_c * MT::invJ_zeta_x(zz(k), eta_x.data()[offset_ij])
           + twiceSyz_c * MT::invJ_zeta_y(zz(k), eta_y.data()[offset_ij]));
    });

  stream1.fence();
}

} // namespace alps::solver::detail
