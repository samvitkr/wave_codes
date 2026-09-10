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
} // namespace

void FlowOverWaveSolver::correct() const
{
  using MT = std::decay_t<decltype(flow_field.mesh)>;
  using Kokkos::parallel_for;

  Kokkos::Profiling::ScopedRegion region("Correct U");

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

  using teampolicy_t = Kokkos::TeamPolicy<Kokkos::IndexType<int>>;
  using member_t     = teampolicy_t::member_type;
  parallel_for(
    "correction_u",
    teampolicy_t(stream1, ends[1] * ends[2], Kokkos::AUTO),
    KOKKOS_LAMBDA(member_t team) {
      int k = team.league_rank() / ends[1];
      int j = team.league_rank() % ends[1];
      parallel_for(
        Kokkos::TeamVectorRange(team, ends[0]), KOKKOS_TR_LAMBDA(int& i) {
          auto inv_dzk  = 1 / dz(k);
          auto inv_dzk1 = 1 / dz(k - 1);
          auto alpha    = 1 / dzw(k - 1);
          auto ratio    = dzw(k - 1) / 2 * inv_dzk;
          auto pw_k     = itp2node(p(i, j, k), p(i, j, k + 1), ratio);
          auto ratio1   = dzw(k - 1) / 2 * inv_dzk1;
          auto pw_k1    = itp2node(p(i, j, k), p(i, j, k - 1), ratio1);
          u(i, j, k) -= (p_x(i, j, k)
                         + (MT::zeta_x(zw(k), exr(i, j)) * pw_k
                            - MT::zeta_x(zw(k - 1), exr(i, j)) * pw_k1)
                             * alpha)
                      * delta_t;
          v(i, j, k) -= (p_y(i, j, k)
                         + (MT::zeta_y(zw(k), eyr(i, j)) * pw_k
                            - MT::zeta_y(zw(k - 1), eyr(i, j)) * pw_k1)
                             * alpha)
                      * delta_t;
          w(i, j, k) -=
            (p(i, j, k + 1) - p(i, j, k)) * inv_dzk * J(i, j) * delta_t;
        });
    });

  // set boundary condition
  apply_bc(stream1);

  stream1.fence();

  ALPS_CHECK_LAST_DEVICE_ERROR();

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

  scale_scalars_by_J();
}

} // namespace alps::solver
