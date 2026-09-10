//
// Created by xuanx004 on 10/4/22.
//

#pragma once

#include <common/container/matrix_field.h>
#include <common/kokkos_abstraction/exec_policy.h>
#include <common/kokkos_abstraction/pool_space.h>
#include <common/real_type.h>
#include <common/runtime/async_utils.h>
#include <decomp/ghost_cell_exchange.h>
#include <solvers/operators/grad.h>
#include <solvers/operators/lap_curvilinear.h>
#include <spectral/spectral.h>

namespace alps::solver {

#define M_INVJ_G11     (invJ(i, j))
#define M_INVJ_G22     (invJ(i, j))
#define M_INVJ_G13(_z) (MT::invJ_zeta_x(_z, eta_x(i, j)))
#define M_INVJ_G23(_z) (MT::invJ_zeta_y(_z, eta_y(i, j)))
#define M_INVJ_G33(_z)                                            \
  (MT::invJ_zeta_x(_z, eta_x(i, j)) * MT::zeta_x(_z, exr(i, j))   \
   + MT::invJ_zeta_y(_z, eta_y(i, j)) * MT::zeta_y(_z, eyr(i, j)) \
   + MT::invJ_zeta_z() * J(i, j))

template<bool exclude_g33 = false, typename FieldType>
void add_viscous_fluxes_from_u(Tensor33Field<Real***> const& fluxes,
                               FieldType const&              flow,
                               Real const                    coeff)
{
  Kokkos::Profiling::pushRegion("viscous");

  using std::tie;
  using boundary_functor_t =
    FDBoundaryFunctor<Real, Kokkos::DefaultExecutionSpace>;

  using MT                = std::decay_t<decltype(flow.mesh)>;
  const auto& mesh        = flow.mesh;
  const auto& grid        = mesh.grid;
  const auto& [u, v, w]   = tie(flow.u.x, flow.u.y, flow.u.z);
  const auto [nx, ny, nz] = local_extents(u);
  const auto is_top       = grid.comm().is_last(2);
  const auto is_bottom    = grid.comm().is_first(2);

  const auto& zz    = mesh.zz;
  const auto& zw    = mesh.zw;
  const auto& dzw   = mesh.dzw;
  const auto& invJ  = mesh.invJ;
  const auto& J     = mesh.J;
  const auto& eta_x = mesh.ex;
  const auto& eta_y = mesh.ey;
  const auto& exr   = mesh.exr;
  const auto& eyr   = mesh.eyr;
  // for whatever reason, when using `coeff' in the following lambda kernels,
  // the kernels compiled by some compilers cannot execute at runtime;
  // workaround by creating a local variable
  const auto nu = coeff;

  if constexpr (!exclude_g33) {
    add_laplacian_fluxes(
      {fluxes.xx, fluxes.xy, fluxes.xz}, u, mesh, nu, "viscous flux u");
    add_laplacian_fluxes(
      {fluxes.yx, fluxes.yy, fluxes.yz}, v, mesh, nu, "viscous flux v");
  } else {
    add_laplacian_fluxes_no_g33(
      {fluxes.xx, fluxes.xy, fluxes.xz}, u, mesh, nu, "viscous flux u");
    add_laplacian_fluxes_no_g33(
      {fluxes.yx, fluxes.yy, fluxes.yz}, v, mesh, nu, "viscous flux v");
  }

  HaloView<Real****, default_memory_pool> const tmp(
    Kokkos::view_alloc("viscous_tmp_derivatives", Kokkos::WithoutInitializing),
    Kokkos::LayoutLeft(u.extent(0), u.extent(1), u.extent(2), 3),
    {begin(u, 0), begin(u, 1), begin(u, 2), 0});
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

  /* Calculate viscous fluxes of w component */
  auto inner = create_inner_view(w).view();

  /* Calculate dw/dξ, dw/d𝜓, dw/dζ */
  spectral::ddx(tmp_x_inner, inner, grid, stream1);
  spectral::ddy(tmp_y_inner, inner, grid, stream2);

  stream1.fence();
  auto reqs1 = async_update_halo_lower_z(grid.x_pencil(), tmp_x, 3);
  stream2.fence();
  auto reqs2 = async_update_halo_lower_z(grid.x_pencil(), tmp_y, 4);

  auto event = get_device_event();
  if (is_bottom) {
    LoopPolicy<2> const      policy(stream2, {0, 0}, {nx, ny});
    boundary_functor_t const functor(
      subview(tmp_z, ALL, ALL, 0).view(),
      subview(w, ALL, ALL, index_range(0, 3)).view(),
      subview(mesh.dzw_h, index_range(0, 2)).view(),
      1,
      LeftBoundary());
    parallel_for("dw/dzeta bottom", policy, functor);
  }
  if (is_top) {
    LoopPolicy<2> const      policy(stream2, {0, 0}, {nx, ny});
    boundary_functor_t const functor(
      subview(tmp_z, ALL, ALL, nz - 1).view(),
      subview(w, ALL, ALL, index_range(nz - 4, nz - 1)).view(),
      subview(mesh.dzw_h, index_range(nz - 4, nz - 2)).view(),
      1,
      RightBoundary());
    parallel_for("dw/dzeta top", policy, functor);
  }
  if (is_bottom || is_top) enqueue(event, stream2);

  auto policy = LoopPolicy<3>(stream1,
                              {0, 0, is_bottom ? 1 : 0},
                              {nx, ny, is_top ? nz - 1 : nz + 1},
                              tile_size);

  const auto& f = w;
  parallel_for(
    "dw/dzeta", policy, KOKKOS_LAMBDA(int i, int j, int k) {
      auto alpha     = 1 / dzw(k - 1);
      tmp_z(i, j, k) = alpha * (f(i, j, k) - f(i, j, k - 1));
    });

  reqs1.waitall();
  reqs2.waitall();
  if (is_bottom || is_top) wait_for(event, stream1);

  /* Calculate F_j = J^{-1} g^{ij} dw/dξ_i */
  const auto& fx = fluxes.zx;
  const auto& fy = fluxes.zy;
  const auto& fz = fluxes.zz;

  if (!check_same_layout_and_offset(fx, fy, fz, tmp_x, tmp_y, tmp_z)) {
    throw std::runtime_error("Mismatched extents in tmp_{xyz} and fluxes");
  }

  const auto functor = KOKKOS_LAMBDA(int i, int j, int k)
  {
    auto offset = &fx(i, j, k) - fx.data();
    auto stride = tmp_z.stride(2);

    auto ratio = dzw(k - 1) / (dzw(k - 1) + dzw(k));
    auto u_zeta =
      itp2node(tmp_z.data()[offset], tmp_z.data()[offset + stride], ratio);
    fx.data()[offset] +=
      (M_INVJ_G11 * tmp_x.data()[offset] + M_INVJ_G13(zw(k)) * u_zeta) * nu;
    fy.data()[offset] +=
      (M_INVJ_G22 * tmp_y.data()[offset] + M_INVJ_G23(zw(k)) * u_zeta) * nu;

    auto u_xi = itp2center(tmp_x.data()[offset], tmp_x.data()[offset - stride]);
    auto u_psi =
      itp2center(tmp_y.data()[offset], tmp_y.data()[offset - stride]);
    // g33_fz is computed outside constexpr-if because nvcc cannot capture them
    // in the constexpr-if below
    auto fz_contrib = M_INVJ_G13(zz(k)) * u_xi + M_INVJ_G23(zz(k)) * u_psi;
    [[maybe_unused]] auto g33_fz = M_INVJ_G33(zz(k)) * tmp_z.data()[offset];
    if constexpr (!exclude_g33) {
      fz_contrib += g33_fz;
    }
    fz.data()[offset] += fz_contrib * nu;
  };
  parallel_for("viscous flux w",
               LoopPolicy<3>(stream1,
                             {0, 0, is_bottom ? 1 : 0},
                             {nx, ny, is_top ? nz - 1 : nz},
                             tile_size),
               functor);
  if (is_bottom) {
    // the flux at the top boundary is probably not used for w, but we compute
    // it anyway for the sake of consistency and debugging
    LoopPolicy<2> const policy2(stream2, {0, 0}, {nx, ny});
    parallel_for(
      "viscous flux w bottom", policy2, KOKKOS_LAMBDA(int i, int j) {
        fx(i, j, 0) +=
          (M_INVJ_G11 * tmp_x(i, j, 0) + M_INVJ_G13(zw(0)) * tmp_z(i, j, 0))
          * nu;
        fy(i, j, 0) +=
          (M_INVJ_G22 * tmp_y(i, j, 0) + M_INVJ_G23(zw(0)) * tmp_z(i, j, 0))
          * nu;
        // g33_fz is computed outside constexpr-if because nvcc cannot capture
        // them in the constexpr-if below
        auto const     zz0 = zz(0);
        auto           fz_contrib =
          M_INVJ_G13(zz0) * tmp_x(i, j, 0) + M_INVJ_G23(zz0) * tmp_y(i, j, 0);
        [[maybe_unused]] auto g33_fz = M_INVJ_G33(zz0) * tmp_z(i, j, 0);
        if constexpr (!exclude_g33) {
          fz_contrib += g33_fz;
        }
        fz(i, j, 0) += fz_contrib * nu;
      });
  }
  if (is_top) {
    // the flux at the top boundary is probably not used for w, but we compute
    // it anyway for the sake of consistency and debugging
    LoopPolicy<2> const policy2(stream2, {0, 0}, {nx, ny});
    const auto          k = nz - 1;
    parallel_for(
      "viscous flux w top", policy2, KOKKOS_LAMBDA(int i, int j) {
        auto fz_contrib = M_INVJ_G13(zz(k)) * tmp_x(i, j, k - 1)
                        + M_INVJ_G23(zz(k)) * tmp_y(i, j, k - 1);
        [[maybe_unused]] auto g33_fz = M_INVJ_G33(zz(k)) * tmp_z(i, j, k);
        if constexpr (!exclude_g33) {
          fz_contrib += g33_fz;
        }
        fz(i, j, k) += fz_contrib * nu;
      });
  }

  stream1.fence();
  stream2.fence();

  Kokkos::Profiling::popRegion();
}

#undef M_INVJ_G11
#undef M_INVJ_G22
#undef M_INVJ_G13
#undef M_INVJ_G23
#undef M_INVJ_G33

} // namespace alps::solver
