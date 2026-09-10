#include "sgs_flux.h"

#include <common/container/view_utils.h>
#include <common/kokkos_abstraction/exec_policy.h>
#include <common/runtime/async_utils.h>
#include <decomp/ghost_cell_exchange.h>
#include <spectral/spectral.h>

#include <Kokkos_Profiling_ScopedRegion.hpp>

namespace alps::solver {
void add_sgs_momentum_fluxes_from_gradu(Tensor33Field<Real***> const& fluxes,
                                        FlowField const&              flow,
                                        HaloView<Real***> const&      nu_t,
                                        Real                          nu,
                                        Vector3Field<Real***> const&  grad_u,
                                        Vector3Field<Real***> const&  grad_v,
                                        Vector3Field<Real***> const&  grad_w,
                                        Real                          gamma)
{
  using Kokkos::parallel_for;
  using std::tie;

  auto const region = Kokkos::Profiling::ScopedRegion("SGS momentum fluxes");

  const auto  coeff = nu * gamma;
  const auto& mesh  = flow.mesh;
  const auto& grid  = mesh.grid;

  const auto stream1 = get_next_stream();
  const auto stream2 = get_next_stream();

  constexpr auto tile_size = []() -> Kokkos::Array<range_index_t, 3> {
    using Device = Kokkos::DefaultExecutionSpace;
    if constexpr (is_cuda_execution_space_v<Device>) return {32, 1, 4};
    if constexpr (is_hip_execution_space_v<Device>) return {64, 1, 4};
    return {0, 0, 0};
  }();

  auto req = async_update_halo_upper_z(grid.partition(), nu_t, 1);

  // Calculate the normal flux and shear flux 𝜈 S₁₂
  {
    const LoopPolicy<3> policy(
      stream1, local_begins(fluxes.xx), local_ends(fluxes.xx), tile_size);
    const auto ux = create_inner_view(grad_u.x).view();
    const auto vy = create_inner_view(grad_v.y).view();
    const auto wz = create_inner_view(grad_w.z).view();
    const auto xx = create_inner_view(fluxes.xx).view();
    const auto yy = create_inner_view(fluxes.yy).view();
    const auto zz = create_inner_view(fluxes.zz).view();
    const auto uy = create_inner_view(grad_u.y).view();
    const auto vx = create_inner_view(grad_v.x).view();
    const auto xy = create_inner_view(fluxes.xy).view();
    const auto yx = create_inner_view(fluxes.yx).view();
    parallel_for(
      "nu*Sii", policy, KOKKOS_LAMBDA(int i, int j, int k) {
        auto offset = &ux(i, j, k) - ux.data();
        auto nu_ijk = nu_t(i, j, k);

        xx.data()[offset] += (coeff + 2 * nu_ijk) * ux.data()[offset];
        yy.data()[offset] += (coeff + 2 * nu_ijk) * vy.data()[offset];
        zz.data()[offset] += (coeff + 2 * nu_ijk) * wz.data()[offset];

        auto twiceS12 = uy.data()[offset] + vx.data()[offset];
        xy.data()[offset] += coeff * uy.data()[offset] + nu_ijk * twiceS12;
        yx.data()[offset] += coeff * vx.data()[offset] + nu_ijk * twiceS12;
      });
  }

  // Calculate the shear flux 𝜈 S₁₃, 𝜈 S₂₃
  // The velocity gradients associated with these terms are available on node
  // locations, and the resulting fluxes are also evaluated on node locations.
  {
    req.waitall(); // nu_t needs interpolation, which needs the upper ghost
                   // cells to be updated

    const auto  nx      = local_extent(fluxes.xz, 0);
    const auto  ny      = local_extent(fluxes.xz, 1);
    const auto  nz      = local_extent(fluxes.xz, 2);
    const auto& dz      = mesh.dz;
    const auto& dzw     = mesh.dzw;
    const auto  uz      = create_inner_view(grad_u.z);
    const auto  wx      = create_inner_view(grad_w.x);
    const auto  vz      = create_inner_view(grad_v.z);
    const auto  wy      = create_inner_view(grad_w.y);
    const auto  f_13    = create_inner_view(fluxes.xz);
    const auto  f_31    = create_inner_view(fluxes.zx);
    const auto  f_23    = create_inner_view(fluxes.yz);
    const auto  f_32    = create_inner_view(fluxes.zy);
    const auto  is_top  = grid.comm().is_last(2);
    const auto  z_begin = 0;
    const auto  z_end   = is_top ? nz - 1 : nz;

    const LoopPolicy<3> policy(
      stream2, {0, 0, z_begin}, {nx, ny, z_end}, tile_size);
    parallel_for(
      "nu*S13", policy, KOKKOS_LAMBDA(int i, int j, int k) {
        auto inv_dzk = 1 / dz(k);
        auto ratio   = dzw(k - 1) / 2 * inv_dzk;
        auto nut_w   = itp2node(nu_t(i, j, k), nu_t(i, j, k + 1), ratio);

        auto offset   = &uz(i, j, k) - uz.data();
        auto twiceS13 = uz.data()[offset] + wx.data()[offset];
        f_13.data()[offset] += coeff * uz.data()[offset] + nut_w * twiceS13;
        f_31.data()[offset] += coeff * wx.data()[offset] + nut_w * twiceS13;
        auto twiceS23 = vz.data()[offset] + wy.data()[offset];
        f_23.data()[offset] += coeff * vz.data()[offset] + nut_w * twiceS23;
        f_32.data()[offset] += coeff * wy.data()[offset] + nut_w * twiceS23;
      });
  }

  stream1.fence();
  stream2.fence();
}

void add_sgs_momentum_fluxes_from_Sij(Tensor33Field<Real***> const&     fluxes,
                                      FlowField const&                  flow,
                                      HaloView<Real***> const&          nu_t,
                                      SymmTensor33Field<Real***> const& Sij)
{
  using Kokkos::parallel_for;
  using std::tie;

  auto const region = Kokkos::Profiling::ScopedRegion("SGS momentum fluxes");

  const auto stream1 = get_next_stream();
  const auto stream2 = get_next_stream();
  const auto stream3 = get_next_stream();

  const auto& mesh = flow.mesh;
  const auto& grid = mesh.grid;

  constexpr auto tile_size = []() -> Kokkos::Array<range_index_t, 3> {
    using Device = Kokkos::DefaultExecutionSpace;
    if constexpr (is_cuda_execution_space_v<Device>) return {32, 1, 4};
    if constexpr (is_hip_execution_space_v<Device>) return {64, 1, 4};
    return {0, 0, 0};
  }();

  auto req = async_update_halo_upper_z(grid.partition(), nu_t, 1);

  // Calculate the normal flux
  {
    const LoopPolicy<3> policy(
      stream1, local_begins(fluxes.xx), local_ends(fluxes.xx));
    const auto Sxx = create_inner_view(Sij.xx);
    const auto xx  = create_inner_view(fluxes.xx);
    const auto Syy = create_inner_view(Sij.yy);
    const auto yy  = create_inner_view(fluxes.yy);
    const auto Szz = create_inner_view(Sij.zz);
    const auto zz  = create_inner_view(fluxes.zz);
    parallel_for(
      "nu*Sii", policy, KOKKOS_LAMBDA(int i, int j, int k) {
        auto offset = &Sxx(i, j, k) - Sxx.data();
        auto nu_ijk = nu_t(i, j, k);

        xx.data()[offset] += nu_ijk * 2 * Sxx.data()[offset];
        yy.data()[offset] += nu_ijk * 2 * Syy.data()[offset];
        zz.data()[offset] += nu_ijk * 2 * Szz.data()[offset];
      });
  }

  // Calculate the shear flux 𝜈 S₁₂
  {
    const LoopPolicy<3> policy(
      stream2, local_begins(fluxes.xy), local_ends(fluxes.xy));
    const auto du_in    = create_inner_view(Sij.xy);
    const auto flux1_in = create_inner_view(fluxes.xy);
    const auto flux2_in = create_inner_view(fluxes.yx);
    parallel_for(
      "nu*ux", policy, KOKKOS_LAMBDA(int i, int j, int k) {
        flux1_in(i, j, k) += nu_t(i, j, k) * 2 * du_in(i, j, k);
        flux2_in(i, j, k) += nu_t(i, j, k) * 2 * du_in(i, j, k);
      });
  }

  // Calculate the shear flux 𝜈 S₁₃, 𝜈 S₂₃
  req.waitall();
  {
    const auto  Sxz    = create_inner_view(Sij.xz);
    const auto  Syz    = create_inner_view(Sij.yz);
    const auto  f_xz   = create_inner_view(fluxes.xz);
    const auto  f_zx   = create_inner_view(fluxes.zx);
    const auto  f_yz   = create_inner_view(fluxes.yz);
    const auto  f_zy   = create_inner_view(fluxes.zy);
    const auto  nx     = local_extent(f_xz, 0);
    const auto  ny     = local_extent(f_xz, 1);
    const auto  nz     = local_extent(f_xz, 2);
    const auto& dz     = mesh.dz;
    const auto& dzw    = mesh.dzw;
    const auto  is_top = grid.comm().is_last(2);

    const auto z_begin = 0;
    const auto z_end   = is_top ? nz - 1 : nz;

    const LoopPolicy<3> policy(
      stream3, {0, 0, z_begin}, {nx, ny, z_end}, tile_size);
    parallel_for(
      "nu*S13", policy, KOKKOS_LAMBDA(int i, int j, int k) {
        auto offset = &Sxz(i, j, k) - Sxz.data();
        auto ratio  = dzw(k - 1) / 2 / dz(k);
        auto nut_w  = itp2node(nu_t(i, j, k), nu_t(i, j, k + 1), ratio);
        f_xz.data()[offset] += nut_w * 2 * Sxz.data()[offset];
        f_zx.data()[offset] += nut_w * 2 * Sxz.data()[offset];
        f_yz.data()[offset] += nut_w * 2 * Syz.data()[offset];
        f_zy.data()[offset] += nut_w * 2 * Syz.data()[offset];
      });
  }

  stream1.fence();
  stream2.fence();
  stream3.fence();
}

} // namespace alps::solver
