#include "grad_curvilinear.h"

#include "grad.h"
#include <common/container/view_utils.h>
#include <common/kokkos_abstraction/exec_policy.h>
#include <common/kokkos_abstraction/pool_space.h>
#include <common/runtime/async_utils.h>
#include <spectral/spectral.h>

#include <Kokkos_Profiling_ScopedRegion.hpp>

namespace alps::solver {

namespace detail {
template<typename MT>
void grad_impl(const Vector3Field<Real***>& grad_u,
               const HaloView<Real***>&     u,
               const MT&                    mesh,
               CenterPt /*tag*/)
{
  auto const region =
    Kokkos::Profiling::ScopedRegion("grad(" + u.label() + ")");

  using Kokkos::parallel_for;

  const auto& grid      = mesh.grid;
  const auto  is_top    = grid.comm().is_last(2);
  const auto  is_bottom = grid.comm().is_first(2);

  const auto nx = grid.extent(0);
  const auto ny = grid.extent(1);
  const auto nz = grid.extent(2);

  const auto stream1 = get_next_stream();
  const auto stream2 = get_next_stream();

  auto u_inner  = create_inner_view(u).view();
  auto ux_inner = create_inner_view(grad_u.x).view();
  auto uy_inner = create_inner_view(grad_u.y).view();
  spectral::ddx(ux_inner, u_inner, grid, stream1);
  spectral::ddy(uy_inner, u_inner, grid, stream2);

  stream1.fence();
  stream2.fence();

  const auto& uz      = grad_u.z;
  const auto  z_begin = is_bottom ? 1 : 0;
  const auto  z_end   = is_top ? nz - 2 : nz;
  const auto& invJ    = mesh.invJ;
  const auto& exr     = mesh.exr;
  const auto& eyr     = mesh.eyr;
  const auto& dz      = mesh.dz;
  const auto& dzw     = mesh.dzw;
  const auto& zz      = mesh.zz;
  if (is_bottom) {
    // Calculate u_x, u_y, u_z at k = 0 (bottom boundary)
    FDBoundaryFunctor<Real, std::decay_t<decltype(stream2)>> const functor(
      subview(uz, ALL, ALL, 0).view(),
      subview(u, ALL, ALL, index_range(0, 3)).view(),
      subview(mesh.dz_h, index_range(0, 2)).view(),
      1.0,
      LeftBoundary());
    parallel_for(LoopPolicy<2>(stream2, {0, 0}, {nx, ny}), functor);
    parallel_for(
      LoopPolicy<2>(stream2, {0, 0}, {nx, ny}), KOKKOS_LAMBDA(int i, int j) {
        ux_inner(i, j, 0) += uz(i, j, 0) * MT::zeta_x(zz(0), exr(i, j));
        uy_inner(i, j, 0) += uz(i, j, 0) * MT::zeta_y(zz(0), eyr(i, j));
        uz(i, j, 0) *= MT::zeta_z(Real(1), invJ(i, j));
      });
  }
  if (is_top) {
    // Calculate u_x, u_y at k = nz - 1 (top boundary), nz - 2 (cell center)
    // Calculate u_z at k = nz - 2 (top boundary)
    FDBoundaryFunctor<Real, std::decay_t<decltype(stream2)>> const functor(
      subview(uz, ALL, ALL, nz - 2).view(),
      subview(u, ALL, ALL, index_range(nz - 3, nz)).view(),
      subview(mesh.dz_h, index_range(nz - 3, nz - 1)).view(),
      1.0,
      RightBoundary());
    parallel_for(LoopPolicy<2>(stream2, {0, 0}, {nx, ny}), functor);
    parallel_for(
      LoopPolicy<2>(stream2, {0, 0}, {nx, ny}), KOKKOS_LAMBDA(int i, int j) {
        const auto inv_dzk1 = 1 / dz(nz - 2 - 1);
        const auto ratio1   = dzw(nz - 2 - 1) / 2 * inv_dzk1;
        const auto u_node_k1 =
          itp2node(u(i, j, nz - 2), u(i, j, nz - 3), ratio1); // u @ node nz-3
        const auto u_zeta_nz2 = (u(i, j, nz - 1) - u_node_k1) / dzw(nz - 3);
        ux_inner(i, j, nz - 2) +=
          u_zeta_nz2 * MT::zeta_x(zz(nz - 2), exr(i, j));
        uy_inner(i, j, nz - 2) +=
          u_zeta_nz2 * MT::zeta_y(zz(nz - 2), eyr(i, j));
        ux_inner(i, j, nz - 1) +=
          uz(i, j, nz - 2) * MT::zeta_x(zz(nz - 1), exr(i, j));
        uy_inner(i, j, nz - 1) +=
          uz(i, j, nz - 2) * MT::zeta_y(zz(nz - 1), eyr(i, j));
        uz(i, j, nz - 2) *= MT::zeta_z(Real(1), invJ(i, j));
      });
  }

  // Calculate u_x, u_y, u_z at all inner points
  GridPolicy<> const policy1(stream1, (z_end - z_begin) * ny, Kokkos::AUTO());
  parallel_for(
    "grad_ddz", policy1, KOKKOS_LAMBDA(GridPolicy<>::member_type const& team) {
      int  k        = team.league_rank() / ny + z_begin;
      int  j        = team.league_rank() % ny;
      auto inv_dzk  = 1 / dz(k);
      auto inv_dzk1 = 1 / dz(k - 1);
      auto alpha    = 1 / dzw(k - 1);
      auto ratio    = dzw(k - 1) / 2 * inv_dzk;
      auto ratio1   = dzw(k - 1) / 2 * inv_dzk1;
      parallel_for(Kokkos::TeamVectorRange(team, nx), [&](int i) {
        auto u_node_k  = itp2node(u(i, j, k), u(i, j, k + 1), ratio);
        auto u_node_k1 = itp2node(u(i, j, k), u(i, j, k - 1), ratio1);
        auto u_zeta    = (u_node_k - u_node_k1) * alpha; // u_zeta @ center k
        ux_inner(i, j, k) += u_zeta * MT::zeta_x(zz(k), exr(i, j));
        uy_inner(i, j, k) += u_zeta * MT::zeta_y(zz(k), eyr(i, j));
        uz(i, j, k) = (u(i, j, k + 1) - u(i, j, k)) * inv_dzk
                    * MT::zeta_z(Real(1), invJ(i, j));
      });
    });

  stream1.fence();
  stream2.fence();
}

template<typename MT>
void grad_impl(const Vector3Field<Real***>& grad_u,
               const HaloView<Real***>&     w,
               const MT&                    mesh,
               NodePt /*tag*/)
{
  auto const region =
    Kokkos::Profiling::ScopedRegion("grad(" + w.label() + ")");

  using Kokkos::parallel_for;

  const auto& grid      = mesh.grid;
  const auto  is_top    = grid.comm().is_last(2);
  const auto  is_bottom = grid.comm().is_first(2);

  const auto nx = grid.extent(0);
  const auto ny = grid.extent(1);
  const auto nz = grid.extent(2);

  const auto stream1 = get_next_stream();
  const auto stream2 = get_next_stream();

  auto w_inner  = create_inner_view(w).view();
  auto wx_inner = create_inner_view(grad_u.x).view();
  auto wy_inner = create_inner_view(grad_u.y).view();
  {
    const auto z_range = std::pair{0, (is_top ? nz - 1 : nz)};
    spectral::ddx(subview(wx_inner, ALL, ALL, z_range),
                  subview(w_inner, ALL, ALL, z_range),
                  grid,
                  stream1);
    spectral::ddy(subview(wy_inner, ALL, ALL, z_range),
                  subview(w_inner, ALL, ALL, z_range),
                  grid,
                  stream2);
  }

  stream1.fence();
  stream2.fence();

  const auto& wz      = grad_u.z;
  const auto  z_begin = is_bottom ? 1 : 0;
  const auto  z_end   = is_top ? nz - 2 : nz;
  const auto& invJ    = mesh.invJ;
  const auto& exr     = mesh.exr;
  const auto& eyr     = mesh.eyr;
  const auto& dz      = mesh.dz;
  const auto& dzw     = mesh.dzw;
  const auto& zw      = mesh.zw;
  if (is_bottom) {
    // Calculate u_x, u_y, u_z at k = 0 (bottom boundary)
    FDBoundaryFunctor<Real, std::decay_t<decltype(stream2)>> const functor(
      subview(wz, ALL, ALL, 0).view(),
      subview(w_inner, ALL, ALL, index_range(0, 3)),
      subview(mesh.dzw_h, index_range(0, 2)).view(),
      1.0,
      LeftBoundary());
    parallel_for(LoopPolicy<2>(stream2, {0, 0}, {nx, ny}), functor);
    parallel_for(LoopPolicy<2>(stream2, {0, 0}, {nx, ny}), functor);
    parallel_for(
      LoopPolicy<2>(stream2, {0, 0}, {nx, ny}), KOKKOS_LAMBDA(int i, int j) {
        wx_inner(i, j, 0) += wz(i, j, 0) * MT::zeta_x(zw(0), exr(i, j));
        wy_inner(i, j, 0) += wz(i, j, 0) * MT::zeta_y(zw(0), eyr(i, j));
        wz(i, j, 0) *= MT::zeta_z(Real(1), invJ(i, j));
      });
  }
  if (is_top) {
    // Calculate w_x, w_y at k = nz - 2 (top boundary)
    // Calculate w_z at k = nz - 1 (top boundary), nz - 2 (cell center)
    FDBoundaryFunctor<Real, std::decay_t<decltype(stream2)>> const functor(
      subview(wz, ALL, ALL, nz - 1).view(),
      subview(w_inner, ALL, ALL, index_range(nz - 4, nz - 1)),
      subview(mesh.dzw_h, index_range(nz - 4, nz - 2)).view(),
      1.0,
      RightBoundary());
    parallel_for(LoopPolicy<2>(stream2, {0, 0}, {nx, ny}), functor);
    parallel_for(
      LoopPolicy<2>(stream2, {0, 0}, {nx, ny}), KOKKOS_LAMBDA(int i, int j) {
        wx_inner(i, j, nz - 2) +=
          wz(i, j, nz - 1) * MT::zeta_x(zw(nz - 2), exr(i, j));
        wy_inner(i, j, nz - 2) +=
          wz(i, j, nz - 1) * MT::zeta_y(zw(nz - 2), eyr(i, j));
        wz(i, j, nz - 1) *= MT::zeta_z(Real(1), invJ(i, j));
        const auto alpha = 1 / dzw(nz - 3);
        wz(i, j, nz - 2) = (w(i, j, nz - 2) - w(i, j, nz - 3)) * alpha
                         * MT::zeta_z(Real(1), invJ(i, j));
      });
  }

  // Calculate w_x, w_y, w_z at all inner points
  GridPolicy<> const policy1(stream1, ny * (z_end - z_begin), Kokkos::AUTO());
  parallel_for(
    "grad_ddz", policy1, KOKKOS_LAMBDA(GridPolicy<>::member_type const& team) {
      int  k       = team.league_rank() / ny + z_begin;
      int  j       = team.league_rank() % ny;
      auto inv_dzk = 1 / dz(k);
      auto alpha   = 1 / dzw(k - 1);
      parallel_for(Kokkos::TeamVectorRange(team, nx), [&](int i) {
        auto w_c_k  = itp2center(w(i, j, k), w(i, j, k - 1)); // w @ center k
        auto w_c_k1 = itp2center(w(i, j, k), w(i, j, k + 1)); // w @ center k+1
        auto u_zeta = (w_c_k1 - w_c_k) * inv_dzk;
        wx_inner(i, j, k) += u_zeta * MT::zeta_x(zw(k), exr(i, j));
        wy_inner(i, j, k) += u_zeta * MT::zeta_y(zw(k), eyr(i, j));
        wz(i, j, k) = (w(i, j, k) - w(i, j, k - 1)) * alpha
                    * MT::zeta_z(Real(1), invJ(i, j));
      });
    });

  stream1.fence();
  stream2.fence();
}

template<typename MT, typename LocType>
Vector3Field<Real***>
grad_impl(const HaloView<Real***>& f, const MT& mesh, LocType /*tag*/)
{
  Vector3Field<Real***, default_memory_pool> grad_u(
    Kokkos::view_alloc("grad(" + f.label() + ")", Kokkos::WithoutInitializing),
    local_extents(f),
    {0, 0, 1});

  detail::grad_impl(grad_u, f, mesh, LocType{});

  return {grad_u.x, grad_u.y, grad_u.z};
}
} // namespace detail

Vector3Field<Real***>
grad(const HaloView<Real***>& f, const BottomWaveMesh& mesh, CenterPt tag)
{
  return detail::grad_impl(f, mesh, tag);
}

Vector3Field<Real***>
grad(const HaloView<Real***>& f, const TopWaveMesh& mesh, CenterPt tag)
{
  return detail::grad_impl(f, mesh, tag);
}

Vector3Field<Real***>
grad(const HaloView<Real***>& f, const BottomWaveMesh& mesh, NodePt tag)
{
  return detail::grad_impl(f, mesh, tag);
}

Vector3Field<Real***>
grad(const HaloView<Real***>& f, const TopWaveMesh& mesh, NodePt tag)
{
  return detail::grad_impl(f, mesh, tag);
}

void grad(const Vector3Field<Real***>& grad_f,
          const HaloView<Real***>&     f,
          const BottomWaveMesh&        mesh,
          CenterPt                     tag)
{
  detail::grad_impl(grad_f, f, mesh, tag);
}

void grad(const Vector3Field<Real***>& grad_f,
          const HaloView<Real***>&     f,
          const TopWaveMesh&           mesh,
          CenterPt                     tag)
{
  detail::grad_impl(grad_f, f, mesh, tag);
}

void grad(const Vector3Field<Real***>& grad_f,
          const HaloView<Real***>&     f,
          const BottomWaveMesh&        mesh,
          NodePt                       tag)
{
  detail::grad_impl(grad_f, f, mesh, tag);
}

void grad(const Vector3Field<Real***>& grad_f,
          const HaloView<Real***>&     f,
          const TopWaveMesh&           mesh,
          NodePt                       tag)
{
  detail::grad_impl(grad_f, f, mesh, tag);
}

} // namespace alps::solver
