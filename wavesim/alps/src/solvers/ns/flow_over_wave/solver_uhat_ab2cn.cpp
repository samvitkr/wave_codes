//
// Created by xuanx004 on 10/4/22.
//

#include "solver.h"

#include <common/base/macros.h>
#include <common/container/view_utils.h>
#include <common/device/devices.h>
#include <common/kokkos_abstraction/exec_policy.h>
#include <common/kokkos_abstraction/pool_space.h>
#include <common/runtime/async_utils.h>
#include <solvers/ns/curvilinear_common/diffusion_cn_zeta_eqn.h>
#include <solvers/operators/grad.h>
#include <spectral/spectral.h>

#include <stdexcept>

namespace alps::solver {

namespace {
void update_mesh_coefficients_from_bc(
  const FlowOverWaveField&             flow,
  const VelocityBC&                    bc,
  const Kokkos::DefaultExecutionSpace& stream)
{
  auto const* ptr_bc = dynamic_cast<const NoSlipWallVarying*>(&bc);

  if (ptr_bc != nullptr) {
    Kokkos::deep_copy(stream, flow.eta, ptr_bc->eta_);
    Kokkos::deep_copy(stream, flow.mesh.et, ptr_bc->eta_t_);

    flow.mesh.update_metric_coefficients(flow.eta, stream);

    stream.fence();
  }
}

void get_cn_rhs(Vector3Field<Real***>&               Ru_s,
                Vector3Field<Real***>                vec_u,
                Vector3Field<Real const***>          Ru,
                Real                                 nu,
                Real                                 dt,
                BottomWaveMesh const&                mesh,
                bool                                 steps_initialized,
                Kokkos::DefaultExecutionSpace const& stream)
{
  using MT = std::decay_t<decltype(mesh)>;
  using boundary_functor_t =
    FDBoundaryFunctor<Real, Kokkos::DefaultExecutionSpace>;

  using Kokkos::parallel_for;
  auto const& u         = vec_u.x;
  auto const& v         = vec_u.y;
  auto const& w         = vec_u.z;
  auto const& invJ      = mesh.invJ;
  auto const& dz        = mesh.dz;
  auto const& dzw       = mesh.dzw;
  auto const& zz        = mesh.zz;
  auto const& zw        = mesh.zw;
  auto const& exr       = mesh.exr;
  auto const& eyr       = mesh.eyr;
  auto const  is_bottom = mesh.grid.comm().is_first(2);
  auto const  is_top    = mesh.grid.comm().is_last(2);

  const auto nx = mesh.extent(0);
  const auto ny = mesh.extent(1);
  const auto nz = mesh.extent(2);

  MDView<Real** [4], default_memory_pool> const boundary_dzeta(
    Kokkos::view_alloc("cn_rhs_boundary_dzeta", Kokkos::WithoutInitializing),
    nx,
    ny);

  auto const du_dzeta_bottom = subview(boundary_dzeta, ALL, ALL, 0);
  auto const du_dzeta_top    = subview(boundary_dzeta, ALL, ALL, 1);
  auto const dv_dzeta_bottom = subview(boundary_dzeta, ALL, ALL, 2);
  auto const dv_dzeta_top    = subview(boundary_dzeta, ALL, ALL, 3);

  auto event = get_device_event();
  enqueue(event, stream);
  auto stream2 = get_next_stream(); // used for w rhs since it does not require
                                    // the du/dzeta, dv/dzeta at the boundaries
  wait_for(event, stream2);

  if (is_bottom) {
    LoopPolicy<2> const      policy(stream, {0, 0}, {nx, ny});
    boundary_functor_t const du_functor(
      du_dzeta_bottom,
      subview(u, ALL, ALL, index_range(0, 3)).view(),
      subview(mesh.dz_h, index_range(0, 2)).view(),
      1,
      LeftBoundary());
    boundary_functor_t const dv_functor(
      dv_dzeta_bottom,
      subview(v, ALL, ALL, index_range(0, 3)).view(),
      subview(mesh.dz_h, index_range(0, 2)).view(),
      1,
      LeftBoundary());
    parallel_for("du/dzeta cn rhs bottom", policy, du_functor);
    parallel_for("dv/dzeta cn rhs bottom", policy, dv_functor);
  }
  if (is_top) {
    LoopPolicy<2> const      policy(stream, {0, 0}, {nx, ny});
    boundary_functor_t const du_functor(
      du_dzeta_top,
      subview(u, ALL, ALL, index_range(nz - 3, nz)).view(),
      subview(mesh.dz_h, index_range(nz - 3, nz - 1)).view(),
      1,
      RightBoundary());
    boundary_functor_t const dv_functor(
      dv_dzeta_top,
      subview(v, ALL, ALL, index_range(nz - 3, nz)).view(),
      subview(mesh.dz_h, index_range(nz - 3, nz - 1)).view(),
      1,
      RightBoundary());
    parallel_for("du/dzeta cn rhs top", policy, du_functor);
    parallel_for("dv/dzeta cn rhs top", policy, dv_functor);
  }

  /*
   * Build CN right-hand side for phi = J^{-1}u:
   *
   *   [I - (dt/2) nu L_zeta] phi^{n+1} = rhs,
   *   rhs = phi^n + dt * R^{AB2} + (dt/2) nu L_zeta(phi^n),
   *   R^{AB2} = (3/2)R^n - (1/2)R^{n-1}.
   *
   * First step fallback (no history): R^{AB2} -> R^n.
   */
  LoopPolicy<3> const policy_center_interior( // excl. boundary points
    stream,
    {0, 0, is_bottom ? 1 : 0},
    {nx, ny, is_top ? nz - 1 : nz});
  LoopPolicy<3> const policy_node_interior( // excl. boundary points
    stream2,
    {0, 0, is_bottom ? 1 : 0},
    {nx, ny, is_top ? nz - 2 : nz});
  if (steps_initialized) {
    const auto alpha = dt / 2;
    parallel_for(
      "w_rhs", policy_node_interior, KOKKOS_LAMBDA(int i, int j, int k) {
        using Kokkos::fma;
        auto invj      = invJ(i, j);
        auto g_l       = MT::g33(zz(k), exr(i, j), eyr(i, j), invj) * invj;
        auto g_u       = MT::g33(zz(k + 1), exr(i, j), eyr(i, j), invj) * invj;
        auto dw_dzeta2 = (g_u * (w(i, j, k + 1) - w(i, j, k)) / dzw(k)
                          - g_l * (w(i, j, k) - w(i, j, k - 1)) / dzw(k - 1))
                       / dz(k);
        auto N          = fma((Real)3, Ru.z(i, j, k), -Ru_s.z(i, j, k));
        Ru_s.z(i, j, k) = fma(fma(nu, dw_dzeta2, N), alpha, w(i, j, k) * invj);
      });
    parallel_for(
      "uv_rhs", policy_center_interior, KOKKOS_LAMBDA(int i, int j, int k) {
        using Kokkos::fma;
        auto invj = invJ(i, j);
        auto g_l  = MT::g33(zw(k - 1), exr(i, j), eyr(i, j), invj) * invj;
        auto g_u  = MT::g33(zw(k), exr(i, j), eyr(i, j), invj) * invj;

        auto du_dzeta_l = (is_bottom && k == 1)
                          ? du_dzeta_bottom(i, j)
                          : (u(i, j, k) - u(i, j, k - 1)) / dz(k - 1);
        auto du_dzeta_u = (is_top && k == nz - 2)
                          ? du_dzeta_top(i, j)
                          : (u(i, j, k + 1) - u(i, j, k)) / dz(k);
        auto du_dzeta2  = (g_u * du_dzeta_u - g_l * du_dzeta_l) / dzw(k - 1);
        auto Nu         = fma((Real)3, Ru.x(i, j, k), -Ru_s.x(i, j, k));
        Ru_s.x(i, j, k) = fma(fma(nu, du_dzeta2, Nu), alpha, u(i, j, k) * invj);

        auto dv_dzeta_l = (is_bottom && k == 1)
                          ? dv_dzeta_bottom(i, j)
                          : (v(i, j, k) - v(i, j, k - 1)) / dz(k - 1);
        auto dv_dzeta_u = (is_top && k == nz - 2)
                          ? dv_dzeta_top(i, j)
                          : (v(i, j, k + 1) - v(i, j, k)) / dz(k);
        auto dv_dzeta2  = (g_u * dv_dzeta_u - g_l * dv_dzeta_l) / dzw(k - 1);
        auto Nv         = fma((Real)3, Ru.y(i, j, k), -Ru_s.y(i, j, k));
        Ru_s.y(i, j, k) = fma(fma(nu, dv_dzeta2, Nv), alpha, v(i, j, k) * invj);
      });
  } else {
    Ru_s.x = Kokkos::create_mirror(default_memory_pool(), Ru.x);
    Ru_s.y = Kokkos::create_mirror(default_memory_pool(), Ru.y);
    Ru_s.z = Kokkos::create_mirror(default_memory_pool(), Ru.z);
    // Forward Euler for the R term
    parallel_for(
      "w_rhs_0", policy_node_interior, KOKKOS_LAMBDA(int i, int j, int k) {
        using Kokkos::fma;
        auto invj      = invJ(i, j);
        auto g_l       = MT::g33(zz(k), exr(i, j), eyr(i, j), invj) * invj;
        auto g_u       = MT::g33(zz(k + 1), exr(i, j), eyr(i, j), invj) * invj;
        auto dw_dzeta2 = (g_u * (w(i, j, k + 1) - w(i, j, k)) / dzw(k)
                          - g_l * (w(i, j, k) - w(i, j, k - 1)) / dzw(k - 1))
                       / dz(k);
        Ru_s.z(i, j, k) =
          fma(fma(nu / 2, dw_dzeta2, Ru.z(i, j, k)), dt, w(i, j, k) * invj);
      });
    parallel_for(
      "uv_rhs_0", policy_center_interior, KOKKOS_LAMBDA(int i, int j, int k) {
        using Kokkos::fma;
        auto invj       = invJ(i, j);
        auto g_l        = MT::g33(zw(k - 1), exr(i, j), eyr(i, j), invj) * invj;
        auto g_u        = MT::g33(zw(k), exr(i, j), eyr(i, j), invj) * invj;
        auto du_dzeta_l = (is_bottom && k == 1)
                          ? du_dzeta_bottom(i, j)
                          : (u(i, j, k) - u(i, j, k - 1)) / dz(k - 1);
        auto du_dzeta_u = (is_top && k == nz - 2)
                          ? du_dzeta_top(i, j)
                          : (u(i, j, k + 1) - u(i, j, k)) / dz(k);
        auto du_dzeta2  = (g_u * du_dzeta_u - g_l * du_dzeta_l) / dzw(k - 1);
        Ru_s.x(i, j, k) =
          fma(fma(nu / 2, du_dzeta2, Ru.x(i, j, k)), dt, u(i, j, k) * invj);

        auto dv_dzeta_l = (is_bottom && k == 1)
                          ? dv_dzeta_bottom(i, j)
                          : (v(i, j, k) - v(i, j, k - 1)) / dz(k - 1);
        auto dv_dzeta_u = (is_top && k == nz - 2)
                          ? dv_dzeta_top(i, j)
                          : (v(i, j, k + 1) - v(i, j, k)) / dz(k);
        auto dv_dzeta2  = (g_u * dv_dzeta_u - g_l * dv_dzeta_l) / dzw(k - 1);
        Ru_s.y(i, j, k) =
          fma(fma(nu / 2, dv_dzeta2, Ru.y(i, j, k)), dt, v(i, j, k) * invj);
      });
  }

  stream2.fence();
  stream.fence(); // ensure work done before deallocating the boundary_dzeta

  ALPS_CHECK_LAST_DEVICE_ERROR();
}

void impose_cn_rhs_bc(Vector3Field<Real***>                rhs,
                      FlowOverWaveField const&             flow,
                      double                               Re,
                      Kokkos::DefaultExecutionSpace const& stream)
{
  auto const is_bottom = flow.mesh.comm().is_first(2);
  auto const is_top    = flow.mesh.comm().is_last(2);
  auto const ends      = local_ends(rhs.x);

  auto const& invJ = flow.mesh.invJ;

  auto const rhs_bottom_x = subview(rhs.x, ALL, ALL, 0).view();
  auto const rhs_bottom_y = subview(rhs.y, ALL, ALL, 0).view();
  auto const rhs_bottom_z = subview(rhs.z, ALL, ALL, 0).view();
  auto const rhs_top_x    = subview(rhs.x, ALL, ALL, ends[2] - 1).view();
  auto const rhs_top_y    = subview(rhs.y, ALL, ALL, ends[2] - 1).view();
  auto const rhs_top_z    = subview(rhs.z, ALL, ALL, ends[2] - 2).view();

  auto const set_bottom_dirichlet = [&](Real u_bc, Real v_bc, Real w_bc) {
    Kokkos::parallel_for(
      LoopPolicy<2>(stream, {0, 0}, {ends[0], ends[1]}),
      KOKKOS_LAMBDA(int i, int j) {
        rhs_bottom_x(i, j) = u_bc * invJ(i, j);
        rhs_bottom_y(i, j) = v_bc * invJ(i, j);
        rhs_bottom_z(i, j) = w_bc * invJ(i, j);
      });
  };

  auto const set_bottom_dirichlet_varying = [&](NoSlipWallVarying const& bc) {
    auto u_bc = bc.u_;
    auto v_bc = bc.v_;
    auto w_bc = bc.w_;
    Kokkos::parallel_for(
      LoopPolicy<2>(stream, {0, 0}, {ends[0], ends[1]}),
      KOKKOS_LAMBDA(int i, int j) {
        rhs_bottom_x(i, j) = u_bc(i, j) * invJ(i, j);
        rhs_bottom_y(i, j) = v_bc(i, j) * invJ(i, j);
        rhs_bottom_z(i, j) = w_bc(i, j) * invJ(i, j);
      });
  };

  auto const set_bottom_neumann = [&](GradientWall const& bc) {
    // du/dz = g1, dv/dz = g2, w = 0 at the bottom boundary
    auto g1 = bc.grad_1;
    auto g2 = bc.grad_2;
    Kokkos::parallel_for(
      LoopPolicy<2>(stream, {0, 0}, {ends[0], ends[1]}),
      KOKKOS_LAMBDA(int i, int j) {
        // du/dzeta = J^{-2} * du/dz
        auto invj          = invJ(i, j);
        rhs_bottom_x(i, j) = g1 * invj * invj;
        rhs_bottom_y(i, j) = g2 * invj * invj;
        rhs_bottom_z(i, j) = 0;
      });
  };

  auto const set_top_dirichlet = [&](Real u_bc, Real v_bc, Real w_bc) {
    Kokkos::parallel_for(
      LoopPolicy<2>(stream, {0, 0}, {ends[0], ends[1]}),
      KOKKOS_LAMBDA(int i, int j) {
        rhs_top_x(i, j) = u_bc * invJ(i, j);
        rhs_top_y(i, j) = v_bc * invJ(i, j);
        rhs_top_z(i, j) = w_bc * invJ(i, j);
      });
  };

  auto const set_top_neumann = [&](Real grad_1, Real grad_2) {
    // du/dz = g1, dv/dz = g2, w = 0 at the top boundary
    Kokkos::parallel_for(
      LoopPolicy<2>(stream, {0, 0}, {ends[0], ends[1]}),
      KOKKOS_LAMBDA(int i, int j) {
        // du/dzeta = J^{-2} * du/dz
        auto invj       = invJ(i, j);
        rhs_top_x(i, j) = grad_1 * invj * invj;
        rhs_top_y(i, j) = grad_2 * invj * invj;
        rhs_top_z(i, j) = 0;
      });
  };

  auto const set_top_dirichlet_varying = [&](NoSlipWallVarying const& bc) {
    auto u_bc = bc.u_;
    auto v_bc = bc.v_;
    auto w_bc = bc.w_;
    Kokkos::parallel_for(
      LoopPolicy<2>(stream, {0, 0}, {ends[0], ends[1]}),
      KOKKOS_LAMBDA(int i, int j) {
        rhs_top_x(i, j) = u_bc(i, j) * invJ(i, j);
        rhs_top_y(i, j) = v_bc(i, j) * invJ(i, j);
        rhs_top_z(i, j) = w_bc(i, j) * invJ(i, j);
      });
  };

  if (is_bottom) {
    bool        bc_set = false;
    auto const* bc_ptr = flow.bottom_bc.get();
    if (auto const* bc = dynamic_cast<NoSlipWallVarying const*>(bc_ptr);
        bc != nullptr) {
      set_bottom_dirichlet_varying(*bc);
      bc_set = true;
    }
    if (auto const* bc = dynamic_cast<NoSlipWall const*>(bc_ptr);
        bc != nullptr) {
      set_bottom_dirichlet(bc->u_, bc->v_, bc->w_);
      bc_set = true;
    }
    if (auto const* bc = dynamic_cast<GradientWall const*>(bc_ptr);
        bc != nullptr) {
      set_bottom_neumann(*bc);
      bc_set = true;
    }
    if (dynamic_cast<TangentialStressWall const*>(bc_ptr) != nullptr) {
      throw std::runtime_error(
        "TangentialStressWall not implemented for bottom boundary");
    }
    if (!bc_set) {
      throw std::runtime_error(
        "Unsupported bottom BC type for AB2CN implicit solve");
    }
  }

  if (is_top) {
    bool        bc_set = false;
    auto const* bc_ptr = flow.top_bc.get();
    if (auto const* bc = dynamic_cast<NoSlipWallVarying const*>(bc_ptr);
        bc != nullptr) {
      set_top_dirichlet_varying(*bc);
      bc_set = true;
    }
    if (auto const* bc = dynamic_cast<NoSlipWall const*>(bc_ptr);
        bc != nullptr) {
      set_top_dirichlet(bc->u_, bc->v_, bc->w_);
      bc_set = true;
    }
    if (auto const* bc = dynamic_cast<GradientWall const*>(bc_ptr);
        bc != nullptr) {
      set_top_neumann(bc->grad_1, bc->grad_2);
      bc_set = true;
    }
    if (auto const* bc = dynamic_cast<TangentialStressWall const*>(bc_ptr);
        bc != nullptr) {
      auto const grad_1 = Real((double)bc->tau_1 * Re);
      auto const grad_2 = Real((double)bc->tau_2 * Re);
      set_top_neumann(grad_1, grad_2);
      bc_set = true;
    }
    if (!bc_set) {
      throw std::runtime_error(
        "Unsupported top BC type for AB2CN implicit solve");
    }
  }

  stream.fence();
}
} // namespace

/*
 * AB2CN velocity update in curvilinear zeta direction.
 *
 * Unknown: phi = J^{-1} uhat.
 * Solve per component:
 *   [I - (dt/2) nu L_zeta] phi^{n+1} = phi^n + dt*R^{AB2} + (dt/2)nu
 * L_zeta(phi^n), where L_zeta(phi) = D_zeta( g^{33} D_zeta phi ), and R
 * contains all explicit terms (advection, pressure gradient, etc) and explicit
 * part of diffusion.
 */
void FlowOverWaveSolver::calc_uhat_ab2cn(
  Vector3Field<Real***>                Ru,
  Kokkos::DefaultExecutionSpace const& stream)
{
  const auto& field = flow_field;

  // Construct and store CN RHS in Ru_saved (Ru_saved currently holds R^{n-1}).
  get_cn_rhs(this->Ru_saved,
             flow_field.u,
             Ru,
             1 / Real(options.Re),
             (Real)dt,
             flow_field.mesh,
             this->steps_initialized,
             stream);

  if (field.bottom_bc->is_time_dependent) {
    auto new_bc     = field.bottom_bc->get_updated_bc(get_time());
    field.bottom_bc = std::move(new_bc);
    update_mesh_coefficients_from_bc(field, *field.bottom_bc, stream);
  }

  impose_cn_rhs_bc(Ru_saved, flow_field, options.Re, stream);

  // Solve the CN system for uhat
  auto const inner_view = [](auto&& v) { return create_inner_view(v).view(); };
  ueqn->solve(inner_view(flow_field.u.x),
              inner_view(Ru_saved.x),
              flow_field.mesh,
              stream);
  ueqn->solve(inner_view(flow_field.u.y),
              inner_view(Ru_saved.y),
              flow_field.mesh,
              stream);

  weqn->solve(inner_view(flow_field.u.z),
              inner_view(Ru_saved.z),
              flow_field.mesh,
              stream);

  stream.fence();

  apply_bc(stream);

  // replace stored history with current explicit RHS: Ru_saved <- R^n.
  std::swap(Ru_saved, Ru);

  ALPS_CHECK_LAST_DEVICE_ERROR();
}

} // namespace alps::solver
