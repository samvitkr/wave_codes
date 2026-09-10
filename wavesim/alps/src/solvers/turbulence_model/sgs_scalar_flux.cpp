//
// Created by xuanx004 on 7/14/24.
//

#include "sgs_scalar_flux.h"

#include <common/container/vector_field.h>
#include <common/container/view_utils.h>
#include <common/kokkos_abstraction/exec_policy.h>
#include <common/kokkos_abstraction/pool_space.h>
#include <common/runtime/async_utils.h>
#include <decomp/ghost_cell_exchange.h>
#include <solvers/field/flow_field.h>
#include <solvers/operators/grad.h>
#include <spectral/spectral.h>

#include <Kokkos_Profiling_ScopedRegion.hpp>

namespace alps::solver {

void add_SGS_and_molecular_diffusion_fluxes(Vector3Field<Real***> const& fluxes,
                                            HaloView<Real const***> const& f,
                                            FlowField const&               flow,
                                            HaloView<Real***> const&       nu_D,
                                            Real const                     D,
                                            Real const gamma)
{
  using Kokkos::parallel_for;

  auto const region =
    Kokkos::Profiling::ScopedRegion("SGS diffusion fluxes " + f.label());

  const Mesh& mesh      = flow.mesh;
  const auto& grid      = mesh.grid;
  const auto  nx        = mesh.extent(0);
  const auto  ny        = mesh.extent(1);
  const auto  nz        = mesh.extent(2);
  const auto  coeff     = D * gamma;
  const auto  is_top    = grid.comm().is_last(2);
  const auto  is_bottom = grid.comm().is_first(2);

  auto stream1 = get_next_stream();
  auto stream2 = get_next_stream();

  MDView<Real***, default_memory_pool> tmp3x(
    Kokkos::view_alloc("tmp_x", Kokkos::WithoutInitializing),
    create_local_layout(grid.pencil));
  MDView<Real***, default_memory_pool> tmp3y(
    Kokkos::view_alloc("tmp_y", Kokkos::WithoutInitializing),
    create_local_layout(grid.pencil));

  constexpr auto tile_size = []() -> Kokkos::Array<std::int64_t, 3> {
    using Device = Kokkos::DefaultExecutionSpace;
    if constexpr (is_cuda_execution_space_v<Device>) return {32, 1, 4};
    if constexpr (is_hip_execution_space_v<Device>) return {64, 1, 4};
    return {0, 0, 0};
  }();

  auto req = async_update_halo_upper_z(grid.partition(), nu_D, 1);

  {
    LoopPolicy<3> axpy_policy_x(
      stream1, local_begins(f), local_ends(f), tile_size);
    LoopPolicy<3> axpy_policy_y(
      stream2, local_begins(f), local_ends(f), tile_size);
    auto       inner        = create_inner_view(f).view();
    auto       flux_x_inner = create_inner_view(fluxes.x).view();
    auto       flux_y_inner = create_inner_view(fluxes.y).view();
    auto const nuD_inner    = create_inner_view(nu_D).view();

    spectral::ddx(tmp3x, inner, grid, stream1);
    Kokkos::parallel_for(
      "nuD*dc/dx", axpy_policy_x, KOKKOS_LAMBDA(int i, int j, int k) {
        flux_x_inner(i, j, k) = Kokkos::fma(
          nuD_inner(i, j, k) + coeff, tmp3x(i, j, k), flux_x_inner(i, j, k));
      });

    spectral::ddy(tmp3y, inner, grid, stream2);
    Kokkos::parallel_for(
      "nuD*dc/dy", axpy_policy_y, KOKKOS_LAMBDA(int i, int j, int k) {
        flux_y_inner(i, j, k) = Kokkos::fma(
          nuD_inner(i, j, k) + coeff, tmp3y(i, j, k), flux_y_inner(i, j, k));
      });
  }
  stream1.fence();

  auto const& flux = fluxes.z;
  // nuD*dc/dz at boundaries
  if (is_bottom) {
    using functor_t  = FDBoundaryFunctor<Real, Kokkos::DefaultExecutionSpace>;
    auto      policy = LoopPolicy<2>(stream1, {0, 0}, {nx, ny});
    functor_t functor(subview(tmp3x, ALL, ALL, 0),
                      subview(f, ALL, ALL, index_range(0, 3)).view(),
                      subview(mesh.dz_h, index_range(0, 2)).view(),
                      mesh.hbar,
                      LeftBoundary());
    parallel_for(policy, functor);

    parallel_for(
      policy, KOKKOS_LAMBDA(int i, int j) {
        flux(i, j, 0) =
          Kokkos::fma(nu_D(i, j, 0) + coeff, tmp3x(i, j, 0), flux(i, j, 0));
      });
  }
  if (is_top) {
    using functor_t  = FDBoundaryFunctor<Real, Kokkos::DefaultExecutionSpace>;
    auto      policy = LoopPolicy<2>(stream1, {0, 0}, {nx, ny});
    functor_t functor(
      subview(tmp3x, ALL, ALL, nz - 2),
      subview(f, ALL, ALL, index_range(nz - 3, nz)).view(),
      subview(mesh.dz_h, index_range(nz - 3, nz - 3 + 2)).view(),
      mesh.hbar,
      RightBoundary());
    parallel_for(policy, functor);

    parallel_for(
      policy, KOKKOS_LAMBDA(int i, int j) {
        flux(i, j, nz - 2) = Kokkos::fma(
          nu_D(i, j, nz - 1) + coeff, tmp3x(i, j, nz - 2), flux(i, j, nz - 2));
      });
  }

  const int   z_begin = is_bottom ? 1 : 0;
  const int   z_end   = is_top ? nz - 2 : nz;
  const auto& dz      = mesh.dz;
  const auto& dzw     = mesh.dzw;
  const auto  hbar    = mesh.hbar;
  const auto  policy =
    LoopPolicy<3>(stream2, {0, 0, z_begin}, {nx, ny, z_end}, tile_size);

  req.waitall();
  parallel_for(
    "diffusion " + f.label() + "_z",
    policy,
    KOKKOS_LAMBDA(int i, int j, int k) {
      auto ratio    = dzw(k - 1) / 2 / dz(k);
      auto nud_w    = itp2node(nu_D(i, j, k), nu_D(i, j, k + 1), ratio);
      auto cz       = (f(i, j, k + 1) - f(i, j, k)) / (dz(k) * hbar);
      flux(i, j, k) = Kokkos::fma(nud_w + coeff, cz, flux(i, j, k));
    });

  stream1.fence();
  stream2.fence();
}

void add_sgs_and_molecular_diffusion_fluxes_from_gradf(
  Vector3Field<Real***> const&       fluxes,
  Vector3Field<Real const***> const& grad_f,
  FlowField const&                   flow,
  HaloView<Real***> const&           nuD,
  Real const                         D,
  Real const                         gamma)
{
  using Kokkos::parallel_for;

  if (nuD.begin(2) > -1) {
    throw std::runtime_error("nuD must have at least one ghost cell in z");
  }

  auto const region =
    Kokkos::Profiling::ScopedRegion("SGS diffusion fluxes " + grad_f.x.label());

  const auto  coeff = D * gamma;
  const auto& mesh  = flow.mesh;

  const auto stream1 = get_next_stream();
  const auto stream2 = get_next_stream();

  constexpr auto tile_size = []() -> Kokkos::Array<std::int64_t, 3> {
    using Device = Kokkos::DefaultExecutionSpace;
    if constexpr (is_cuda_execution_space_v<Device>) return {32, 1, 4};
    if constexpr (is_hip_execution_space_v<Device>) return {64, 1, 4};
    return {0, 0, 0};
  }();

  auto req = async_update_halo_upper_z(mesh.partition(), nuD, 1);

  // Calculate q_x, q_y
  {
    const LoopPolicy<3> policy(
      stream1, local_begins(fluxes.x), local_ends(fluxes.x), tile_size);
    const auto fx = create_inner_view(grad_f.x).view();
    const auto fy = create_inner_view(grad_f.y).view();
    const auto qx = create_inner_view(fluxes.x).view();
    const auto qy = create_inner_view(fluxes.y).view();
    parallel_for(
      "nuD*grad(f)", policy, KOKKOS_LAMBDA(int i, int j, int k) {
        qx(i, j, k) += (coeff + nuD(i, j, k)) * fx(i, j, k);
        qy(i, j, k) += (coeff + nuD(i, j, k)) * fy(i, j, k);
      });
  }

  // Calculate q_z
  // The zeta-gradient is computed on node locations, and the resulting flux is
  // also evaluated on node locations.
  {
    req.waitall(); // nuD needs interpolation, which requires the upper ghost
                   // cells to be updated

    const auto [nx, ny, nz] = local_extents(fluxes.z);
    const auto& dz          = mesh.dz;
    const auto& dzw         = mesh.dzw;
    const auto  fz          = create_inner_view(grad_f.z);
    const auto  qz          = create_inner_view(fluxes.z);
    const auto  is_top      = mesh.comm().is_last(2);
    const auto  z_begin     = 0;
    const auto  z_end       = is_top ? nz - 1 : nz;

    const LoopPolicy<3> policy(
      stream2, {0, 0, z_begin}, {nx, ny, z_end}, tile_size);
    parallel_for(
      "nu*S13", policy, KOKKOS_LAMBDA(int i, int j, int k) {
        auto inv_dzk = 1 / dz(k);
        auto ratio   = dzw(k - 1) / 2 * inv_dzk;
        auto nuD_w   = itp2node(nuD(i, j, k), nuD(i, j, k + 1), ratio);
        qz(i, j, k) += (coeff + nuD_w) * fz(i, j, k);
      });
  }

  stream1.fence();
  stream2.fence();
}

} // namespace alps::solver
