//
// Created by xuanx004 on 10/4/22.
//

#include "solver.h"

#include <common/base/macros.h>
#include <common/device/devices.h>
#include <common/kokkos_abstraction/exec_policy.h>
#include <common/runtime/async_utils.h>
#include <decomp/ghost_cell_exchange.h>
#include <spectral/spectral.h>

#include <Kokkos_Profiling_ScopedRegion.hpp>

namespace alps::solver {

namespace {
struct scale_J_functor
{
  MDView<Real***> f;
  MDView<Real**>  J;

  scale_J_functor(MDView<Real***> f_, MDView<Real**> J_)
    : f{std::move(f_)}
    , J{std::move(J_)}
  {}

  KOKKOS_FUNCTION void operator()(int i, int j, int k) const
  {
    f(i, j, k) *= J(i, j);
  }
};

template<typename T>
KOKKOS_FORCEINLINE_FUNCTION constexpr auto zetax(T z, T ex)
{
  return TopWaveMesh::zeta_x(z, ex);
}

template<typename T>
KOKKOS_FORCEINLINE_FUNCTION constexpr auto zetay(T z, T ey)
{
  return TopWaveMesh::zeta_y(z, ey);
}
} // namespace

void FreeSurfaceSolver::correct(int rk_stage) const
{
  using Kokkos::parallel_for;

  Kokkos::Profiling::ScopedRegion region("Correct U");

  const auto  is_top  = flow_field.mesh.comm().is_last(2);
  const auto  delta_t = (Real)this->dt;
  const auto& u       = flow_field.u.x;
  const auto& v       = flow_field.u.y;
  const auto& w       = flow_field.u.z;
  const auto& p       = flow_field.pp;
  const auto& mesh    = flow_field.mesh;
  const auto& zw      = mesh.zw;
  const auto& dzw     = mesh.dzw;
  const auto& dz      = mesh.dz;
  const auto& exr     = mesh.exr;
  const auto& eyr     = mesh.eyr;
  const auto& J       = mesh.J;

  auto ends = local_extents(u);

  // update the ghost cells of pressure
  auto reqs = async_update_halo_z(mesh.grid.x_pencil(), p, 1);

  const auto p_inner = create_inner_view(p).view();

  MDView<Real***, default_memory_pool> p_x(
    Kokkos::view_alloc("dpx", Kokkos::WithoutInitializing), p_inner.layout());
  MDView<Real***, default_memory_pool> p_y(
    Kokkos::view_alloc("dpy", Kokkos::WithoutInitializing), p_inner.layout());

  auto stream1 = get_next_stream();
  auto stream2 = get_next_stream();

  spectral::ddx(p_x, p_inner, mesh.grid, stream1);
  reqs.waitall(); // overlap ghost cell exchange with ddx

  spectral::ddy(p_y, p_inner, mesh.grid, stream2);

  stream1.fence();
  stream2.fence();

  // auto tile_size = is_cuda_execution_space_v<Kokkos::DefaultExecutionSpace>
  //                  ? Kokkos::Array<int, 3>{16, 4, 2}
  //                  : Kokkos::Array<int, 3>{0, 0, 0};
  parallel_for(
    "correction_u",
    // LoopPolicy<3>(stream1,
    //               {0, 0, 0},
    //               {ends[0], ends[1], is_top ? ends[2] - 2 : ends[2]},
    //               tile_size),
    // KOKKOS_LAMBDA(int i, int j, int k) {
    GridPolicy<>(
      stream1, ends[1] * (is_top ? ends[2] - 2 : ends[2]), Kokkos::AUTO()),
    KOKKOS_LAMBDA(GridPolicy<>::member_type const& team) {
      int  k        = team.league_rank() / ends[1];
      int  j        = team.league_rank() % ends[1];
      auto inv_dzk  = 1 / dz(k);
      auto inv_dzk1 = 1 / dz(k - 1);
      auto alpha    = 1 / dzw(k - 1);
      auto ratio    = dzw(k - 1) / 2 * inv_dzk;
      auto ratio1   = dzw(k - 1) / 2 * inv_dzk1;
      parallel_for(
        Kokkos::TeamThreadRange(team, ends[0]), KOKKOS_TR_LAMBDA(int& i) {
          auto pw_k     = itp2node(p(i, j, k), p(i, j, k + 1), ratio);
          auto pw_k1    = itp2node(p(i, j, k), p(i, j, k - 1), ratio1);
          auto zetax_k  = zetax(zw(k), exr(i, j));
          auto zetax_k1 = zetax(zw(k - 1), exr(i, j));
          u(i, j, k) -=
            (p_x(i, j, k) + (zetax_k * pw_k - zetax_k1 * pw_k1) * alpha)
            * delta_t;
          auto zetay_k  = zetay(zw(k), eyr(i, j));
          auto zetay_k1 = zetay(zw(k - 1), eyr(i, j));
          v(i, j, k) -=
            (p_y(i, j, k) + (zetay_k * pw_k - zetay_k1 * pw_k1) * alpha)
            * delta_t;
          w(i, j, k) -=
            (p(i, j, k + 1) - p(i, j, k)) * inv_dzk * J(i, j) * delta_t;
        });
    });
  if (is_top) {
    int const   nz = ends[2];
    auto const& zz = mesh.zz;
    parallel_for(
      "correction_u top",
      LoopPolicy<2>(stream2, {0, 0}, {ends[0], ends[1]}),
      KOKKOS_LAMBDA(int i, int j) {
        auto dz0   = dz(nz - 2);
        auto dz1   = dz(nz - 3);
        auto exrij = exr(i, j);
        auto eyrij = eyr(i, j);

        auto alpha       = dz0 / (dz0 + dz1);
        auto coeff2      = (1 + alpha) / dz0;
        auto coeff1      = -1 / alpha / dz1;
        auto coeff0      = alpha / dz1;
        auto zetax_pzeta = zetax(zz(nz - 3), exrij) * p(i, j, nz - 3) * coeff0
                         + zetax(zz(nz - 2), exrij) * p(i, j, nz - 2) * coeff1
                         + zetax(zz(nz - 1), exrij) * p(i, j, nz - 1) * coeff2;
        auto zetay_pzeta = zetay(zz(nz - 3), eyrij) * p(i, j, nz - 3) * coeff0
                         + zetay(zz(nz - 2), eyrij) * p(i, j, nz - 2) * coeff1
                         + zetay(zz(nz - 1), eyrij) * p(i, j, nz - 1) * coeff2;
        auto p_zeta = p(i, j, nz - 3) * coeff0 + p(i, j, nz - 2) * coeff1
                    + p(i, j, nz - 1) * coeff2;

        u(i, j, nz - 1) -= (p_x(i, j, nz - 1) + zetax_pzeta) * delta_t;
        v(i, j, nz - 1) -= (p_y(i, j, nz - 1) + zetay_pzeta) * delta_t;
        w(i, j, nz - 2) -= p_zeta * J(i, j) * delta_t;

        auto pw_k     = p(i, j, nz - 1);
        auto ratio1   = dzw(nz - 3) / 2 / dz1;
        auto pw_k1    = itp2node(p(i, j, nz - 2), p(i, j, nz - 3), ratio1);
        auto zetax_k  = zetax(zw(nz - 2), exrij);
        auto zetax_k1 = zetax(zw(nz - 3), exrij);
        u(i, j, nz - 2) -= (p_x(i, j, nz - 2)
                            + (zetax_k * pw_k - zetax_k1 * pw_k1) / dzw(nz - 3))
                         * delta_t;
        auto zetay_k  = zetay(zw(nz - 2), eyrij);
        auto zetay_k1 = zetay(zw(nz - 3), eyrij);
        v(i, j, nz - 2) -= (p_y(i, j, nz - 2)
                            + (zetay_k * pw_k - zetay_k1 * pw_k1) / dzw(nz - 3))
                         * delta_t;
      });
  }

  stream1.fence();
  stream2.fence();

  ALPS_CHECK_LAST_DEVICE_ERROR();

  // set boundary condition
  // flow_field.top_bc->get_surface_invJw(
  //   subview(w, ALL, ALL, index_range(ends[2] - 3, ends[2] - 1)).view(),
  //   subview(u, ALL, ALL, index_range(ends[2] - 3, ends[2])).view(),
  //   subview(v, ALL, ALL, index_range(ends[2] - 3, ends[2])).view(),
  //   mesh,
  //   stream1);
  apply_bottom_bc(stream1);

  // dealias and update ghost cells
  LoopPolicy<3> norm_policy(stream1, local_begins(u), local_ends(u));
  auto          event_u = get_device_event();
  auto          event_v = get_device_event();
  auto          event_w = get_device_event();

  auto u_inner = create_inner_view(u).view();
  spectral::dealias(u_inner, mesh.grid, stream1);
  parallel_for("normalize u", norm_policy, scale_J_functor(u_inner, J));
  enqueue(event_u, stream1);

  auto v_inner = create_inner_view(v).view();
  spectral::dealias(v_inner, mesh.grid, stream1);
  parallel_for("normalize v", norm_policy, scale_J_functor(v_inner, J));
  enqueue(event_v, stream1);

  wait_for(event_u);
  update_halo_z(mesh.grid.x_pencil(), u, 2);

  auto w_inner = create_inner_view(w).view();
  spectral::dealias(w_inner, mesh.grid, stream1);
  parallel_for("normalize w", norm_policy, scale_J_functor(w_inner, J));
  enqueue(event_w, stream1);

  wait_for(event_v);
  update_halo_z(mesh.grid.x_pencil(), v, 3);

  parallel_for(
    "normalize p",
    LoopPolicy<3>(stream1, begins(p), alps::ends(p)),
    KOKKOS_LAMBDA(int i, int j, int k) { p(i, j, k) *= J(i, j); });

  wait_for(event_w);
  update_halo_z(mesh.grid.x_pencil(), w, 4);

  stream1.fence();

  ALPS_CHECK_LAST_DEVICE_ERROR();

  if (rk_stage == 2) {
    scale_scalars_by_J();
  }
}

} // namespace alps::solver
