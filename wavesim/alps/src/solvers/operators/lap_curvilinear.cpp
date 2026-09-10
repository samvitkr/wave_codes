#include "lap_curvilinear.h"

#include <common/container/view_utils.h>
#include <common/kokkos_abstraction/exec_policy.h>
#include <common/kokkos_abstraction/pool_space.h>
#include <common/runtime/async_utils.h>
#include <decomp/ghost_cell_exchange.h>
#include <solvers/operators/grad.h>
#include <spectral/spectral.h>

#include <Kokkos_Profiling_ScopedRegion.hpp>

#include <stdexcept>
#include <string>

namespace alps::solver {

namespace {

#define M_INVJ_G11     (invJ(i, j))
#define M_INVJ_G22     (invJ(i, j))
#define M_INVJ_G13(_z) (MT::invJ_zeta_x(_z, eta_x(i, j)))
#define M_INVJ_G23(_z) (MT::invJ_zeta_y(_z, eta_y(i, j)))
#define M_INVJ_G33(_z)                                            \
  (MT::invJ_zeta_x(_z, eta_x(i, j)) * MT::zeta_x(_z, exr(i, j))   \
   + MT::invJ_zeta_y(_z, eta_y(i, j)) * MT::zeta_y(_z, eyr(i, j)) \
   + MT::invJ_zeta_z() * J(i, j))

template<bool exclude_g33, typename MT>
void add_laplacian_fluxes_impl(const Vector3Field<Real***>&   fluxes,
                               const HaloView<Real const***>& f,
                               const MT&                      mesh,
                               Real const                     coeff,
                               std::string_view               label)
{
  auto const region = Kokkos::Profiling::ScopedRegion(std::string(label));

  using Kokkos::parallel_for;
  using boundary_functor_t =
    FDBoundaryFunctor<Real, Kokkos::DefaultExecutionSpace>;

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
  const auto  D     = coeff;

  HaloView<Real****, default_memory_pool> const tmp(
    Kokkos::view_alloc("lap_tmp_derivatives", Kokkos::WithoutInitializing),
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

  auto inner = create_inner_view(f).view();

  spectral::ddx(tmp_x_inner, inner, grid, stream1); // tmp_x = df/dξ
  spectral::ddy(tmp_y_inner, inner, grid, stream2); // tmp_y = df/d𝜓

  // Ghost cells update must be issued after the calculations complete
  stream1.fence();
  auto reqs1 = async_update_halo_upper_z(grid.x_pencil(), tmp_x, 1);
  stream2.fence();
  auto reqs2 = async_update_halo_upper_z(grid.x_pencil(), tmp_y, 2);

  auto event = get_device_event();
  if (is_bottom) {
    auto const               policy = LoopPolicy<2>(stream2, {0, 0}, {nx, ny});
    boundary_functor_t const functor(
      subview(tmp_z, ALL, ALL, 0).view(),
      subview(f, ALL, ALL, index_range(0, 3)).view(),
      subview(mesh.dz_h, index_range(0, 2)).view(),
      1,
      LeftBoundary());
    parallel_for("df/dzeta bottom", policy, functor);
  }
  if (is_top) {
    auto const               policy = LoopPolicy<2>(stream2, {0, 0}, {nx, ny});
    boundary_functor_t const functor(
      subview(tmp_z, ALL, ALL, nz - 2).view(),
      subview(f, ALL, ALL, index_range(nz - 3, nz)).view(),
      subview(mesh.dz_h, index_range(nz - 3, nz - 1)).view(),
      1,
      RightBoundary());
    parallel_for("df/dzeta top", policy, functor);
  }
  if (is_bottom || is_top) enqueue(event, stream2);

  // Compute df/dζ in the interior
  auto policy = LoopPolicy<3>(stream1,
                              {0, 0, is_bottom ? 1 : -1},
                              {nx, ny, is_top ? nz - 2 : nz},
                              tile_size);
  parallel_for(
    "df/dzeta", policy, KOKKOS_LAMBDA(int i, int j, int k) {
      auto alpha     = 1 / dz(k);
      tmp_z(i, j, k) = alpha * (f(i, j, k + 1) - f(i, j, k));
    });

  reqs1.waitall();
  reqs2.waitall();
  if (is_bottom || is_top) wait_for(event, stream1);

  /* Calculate F_j = J^{-1} g^{ij} df/dξ_i */
  const auto& fx = fluxes.x;
  const auto& fy = fluxes.y;
  const auto& fz = fluxes.z;

  if (!check_same_layout_and_offset(fx, fy, fz, tmp_x, tmp_y, tmp_z)) {
    throw std::runtime_error("Mismatched layout and offset in lap fluxes");
  }

  const auto functor = KOKKOS_LAMBDA(int i, int j, int k)
  {
    auto offset = &fx(i, j, k) - fx.data();
    auto stride = tmp_z.stride(2);

    // When is_bottom is true, this reads an uninitialized value at k = -1,
    // making f_zeta_c, fx, and fy invalid at k = 0.
    // This is harmless here because the fx, fy fluxes at k = 0 should not be
    // used in the solver
    auto f_zeta_c =
      itp2center(tmp_z.data()[offset - stride], tmp_z.data()[offset]);
    fx.data()[offset] +=
      (M_INVJ_G11 * tmp_x.data()[offset] + M_INVJ_G13(zz(k)) * f_zeta_c) * D;
    fy.data()[offset] +=
      (M_INVJ_G22 * tmp_y.data()[offset] + M_INVJ_G23(zz(k)) * f_zeta_c) * D;
    auto ratio = dzw(k - 1) / (dzw(k - 1) + dzw(k));
    auto f_xi =
      itp2node(tmp_x.data()[offset], tmp_x.data()[offset + stride], ratio);
    auto f_psi =
      itp2node(tmp_y.data()[offset], tmp_y.data()[offset + stride], ratio);
    // g33_fz is computed outside constexpr-if because nvcc cannot capture them
    // in the constexpr-if below
    auto fz_contrib = M_INVJ_G13(zw(k)) * f_xi + M_INVJ_G23(zw(k)) * f_psi;
    [[maybe_unused]] auto g33_fz = M_INVJ_G33(zw(k)) * tmp_z.data()[offset];
    if constexpr (!exclude_g33) {
      fz_contrib += g33_fz;
    }
    fz.data()[offset] += fz_contrib * D;
  };
  parallel_for("lap f",
               LoopPolicy<3>(
                 stream1, {0, 0, 0}, {nx, ny, is_top ? nz - 1 : nz}, tile_size),
               functor);

  stream1.fence();
  stream2.fence();
}

#undef M_INVJ_G11
#undef M_INVJ_G22
#undef M_INVJ_G13
#undef M_INVJ_G23
#undef M_INVJ_G33

} // anonymous namespace

void add_laplacian_fluxes(const Vector3Field<Real***>&   fluxes,
                          const HaloView<Real const***>& f,
                          const BottomWaveMesh&          mesh,
                          Real const                     coeff,
                          std::string_view               label)
{
  add_laplacian_fluxes_impl<false>(fluxes, f, mesh, coeff, label);
}

void add_laplacian_fluxes_no_g33(const Vector3Field<Real***>&   fluxes,
                                 const HaloView<Real const***>& f,
                                 const BottomWaveMesh&          mesh,
                                 Real const                     coeff,
                                 std::string_view               label)
{
  add_laplacian_fluxes_impl<true>(fluxes, f, mesh, coeff, label);
}

void add_laplacian_fluxes(const Vector3Field<Real***>&   fluxes,
                          const HaloView<Real const***>& f,
                          const TopWaveMesh&             mesh,
                          Real const                     coeff,
                          std::string_view               label)
{
  add_laplacian_fluxes_impl<false>(fluxes, f, mesh, coeff, label);
}

} // namespace alps::solver
