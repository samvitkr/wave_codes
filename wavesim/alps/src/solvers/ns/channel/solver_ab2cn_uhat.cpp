#include "solver_ab2cn.h"

#include "convection_diffusion.h"
#include <common/base/logging.h>
#include <common/base/macros.h>
#include <common/container/view_utils.h>
#include <common/device/devices.h>
#include <common/kokkos_abstraction/exec_policy.h>
#include <common/kokkos_abstraction/pool_space.h>
#include <common/runtime/async_utils.h>
#include <solvers/field/flow_field.h>
#include <solvers/ns/channel/diffusion_cn.h>
#include <solvers/operators/div.h>
#include <solvers/poisson/tridiagonal_solver.h>
#include <spectral/spectral.h>

#include <Kokkos_Core.hpp>
#include <Kokkos_Profiling_ScopedRegion.hpp>

#include <stdexcept>
#include <utility>

namespace alps::solver {

void ChannelFlowSolverAB2CN::solve_ueqn() const
{
  using Kokkos::parallel_for;
  using Kokkos::TeamThreadRange;

  auto const& grid      = flow_field.mesh.grid;
  auto const& u         = flow_field.u;
  auto const  is_top    = grid.comm().is_last(2);
  auto const  is_bottom = grid.comm().is_first(2);

  auto stream = get_next_stream();
  impose_u_eqn_bottom_bc(stream);
  impose_u_eqn_top_bc(stream);

  if (is_bottom) stream.fence();
  if (is_top) stream.fence();

  MDView<Real***, default_memory_pool> sol(
    Kokkos::view_alloc("", Kokkos::WithoutInitializing),
    ueqn->coeff_d.layout());

  spectral::fft_r2c_xy(sol, create_inner_view(u.x).view(), grid, stream);
  ueqn->solve(sol, stream);
  spectral::fft_c2r_xy(create_inner_view(u.x).view(), sol, true, grid, stream);

  spectral::fft_r2c_xy(sol, create_inner_view(u.y).view(), grid, stream);
  ueqn->solve(sol, stream);
  spectral::fft_c2r_xy(create_inner_view(u.y).view(), sol, true, grid, stream);

  // w component has a different layout
  auto w_subview = [w_n = weqn->coeff_d.extent_int(2)](auto const& f) {
    return subview(f, Kokkos::ALL, Kokkos::ALL, std::make_pair(0, w_n));
  };
  spectral::fft_r2c_xy(w_subview(sol), w_subview(u.z).view(), grid, stream);
  weqn->solve(w_subview(sol), stream);
  spectral::fft_c2r_xy(
    w_subview(u.z).view(), w_subview(sol), true, grid, stream);
  stream.fence();
}

void ChannelFlowSolverAB2CN::impose_u_eqn_bottom_bc(
  const Kokkos::DefaultExecutionSpace& space) const
{
  auto const is_bottom = flow_field.mesh.grid.comm().is_first(2);

  if (!is_bottom) return;

  auto ends = local_extents(flow_field.u.x);
  auto u    = subview(flow_field.u.x, ALL, ALL, 0);
  auto v    = subview(flow_field.u.y, ALL, ALL, 0);
  auto w    = subview(flow_field.u.z, ALL, ALL, 0);
  auto u1   = subview(flow_field.u.x, ALL, ALL, 1);
  auto v1   = subview(flow_field.u.y, ALL, ALL, 1);
  auto Dz   = flow_field.mesh.dzw_h(0) * flow_field.mesh.hbar;

  bool  bc_set   = false;
  auto* base_ptr = flow_field.bottom_bc.get();
  if (auto const* bc = dynamic_cast<NoSlipWall*>(base_ptr); bc != nullptr) {
    auto bc_u = bc->u_;
    auto bc_v = bc->v_;
    auto bc_w = bc->w_;
    Kokkos::parallel_for(
      LoopPolicy<2>(space, {0, 0}, {ends[0], ends[1]}),
      KOKKOS_LAMBDA(int i, int j) {
        u(i, j) = bc_u;
        v(i, j) = bc_v;
        w(i, j) = bc_w;
      });
    bc_set = true;
  }
  if (auto const* bc = dynamic_cast<GradientWall*>(base_ptr); bc != nullptr) {
    auto bc_coeff = Real(ueqn->alpha / double(Dz));
    auto grad_1   = bc->grad_1;
    auto grad_2   = bc->grad_2;
    Kokkos::parallel_for(
      LoopPolicy<2>(space, {0, 0}, {ends[0], ends[1]}),
      KOKKOS_LAMBDA(int i, int j) {
        u(i, j) = grad_1;
        v(i, j) = grad_2;
        w(i, j) = 0;
        u1(i, j) -= bc_coeff * grad_1;
        v1(i, j) -= bc_coeff * grad_2;
      });
    bc_set = true;
  }
  if (auto const* bc = dynamic_cast<TangentialStressWall*>(base_ptr);
      bc != nullptr) {
    throw std::runtime_error(
      "TangentialStressWall not implemented for bottom boundary");
  }

  if (!bc_set) {
    logger->warn("No bottom boundary condition set");
  }
}

void ChannelFlowSolverAB2CN::impose_u_eqn_top_bc(
  const Kokkos::DefaultExecutionSpace& space) const
{
  auto const is_top = flow_field.mesh.grid.comm().is_last(2);

  if (!is_top) return;

  auto ends = local_extents(flow_field.u.x);
  auto u    = subview(flow_field.u.x, ALL, ALL, ends[2] - 1);
  auto v    = subview(flow_field.u.y, ALL, ALL, ends[2] - 1);
  auto w    = subview(flow_field.u.z, ALL, ALL, ends[2] - 2);
  auto u1   = subview(flow_field.u.x, ALL, ALL, ends[2] - 2);
  auto v1   = subview(flow_field.u.y, ALL, ALL, ends[2] - 2);
  auto Dz   = flow_field.mesh.dzw_h(ends[2] - 3) * flow_field.mesh.hbar;

  bool  bc_set   = false;
  auto* base_ptr = flow_field.top_bc.get();
  if (auto const* bc = dynamic_cast<NoSlipWall*>(base_ptr); bc != nullptr) {
    auto bc_u = bc->u_;
    auto bc_v = bc->v_;
    auto bc_w = bc->w_;
    Kokkos::parallel_for(
      LoopPolicy<2>(space, {0, 0}, {ends[0], ends[1]}),
      KOKKOS_LAMBDA(int i, int j) {
        u(i, j) = bc_u;
        v(i, j) = bc_v;
        w(i, j) = bc_w;
      });
    bc_set = true;
  }
  if (auto const* bc = dynamic_cast<GradientWall*>(base_ptr); bc != nullptr) {
    auto bc_coeff = Real(ueqn->alpha / double(Dz));
    auto grad_1   = bc->grad_1;
    auto grad_2   = bc->grad_2;
    Kokkos::parallel_for(
      LoopPolicy<2>(space, {0, 0}, {ends[0], ends[1]}),
      KOKKOS_LAMBDA(int i, int j) {
        u(i, j) = grad_1;
        v(i, j) = grad_2;
        w(i, j) = 0;
        u1(i, j) += bc_coeff * grad_1;
        v1(i, j) += bc_coeff * grad_2;
      });
    bc_set = true;
  }
  if (auto const* bc = dynamic_cast<TangentialStressWall*>(base_ptr);
      bc != nullptr) {
    auto bc_coeff = Real(ueqn->alpha / double(Dz));
    auto grad_1   = Real((double)bc->tau_1 * options.Re);
    auto grad_2   = Real((double)bc->tau_2 * options.Re);
    Kokkos::parallel_for(
      LoopPolicy<2>(space, {0, 0}, {ends[0], ends[1]}),
      KOKKOS_LAMBDA(int i, int j) {
        u(i, j) = grad_1;
        v(i, j) = grad_2;
        w(i, j) = 0;
        u1(i, j) += bc_coeff * grad_1;
        v1(i, j) += bc_coeff * grad_2;
      });
    bc_set = true;
  }

  if (!bc_set) {
    logger->warn("No top boundary condition set.");
  }
}

void ChannelFlowSolverAB2CN::calc_uhat(URhs rhs)
{
  using Kokkos::fma;
  using Kokkos::parallel_for;
  using Kokkos::TeamVectorRange;

  auto& R   = rhs.Ru;
  auto& Rnu = rhs.Ru_viscous;

  const auto& u = flow_field.u.x;
  const auto& v = flow_field.u.y;
  const auto& w = flow_field.u.z;

  auto ends = local_extents(u);

  Kokkos::Profiling::pushRegion("U stepping");
  // Advance the velocity field
  auto policy =
    GridPolicy<>(get_next_stream(), ends[1] * ends[2], Kokkos::AUTO);
  using member_t = GridPolicy<>::member_type;
  if (steps_initialized) {
    // Second-order Adam Bashforth
    const auto  alpha = Real(dt / 2);
    const auto& Rs    = this->Ru_saved;
    parallel_for(
      "u_stepping", policy, KOKKOS_LAMBDA(member_t team) {
        int k = team.league_rank() / ends[1];
        int j = team.league_rank() % ends[1];
        parallel_for(TeamVectorRange(team, ends[0]), [&](int i) {
          auto f = fma((Real)3, R.x(i, j, k), -Rs.x(i, j, k) + Rnu.x(i, j, k));
          u(i, j, k) = fma(alpha, f, u(i, j, k));
        });
        parallel_for(TeamVectorRange(team, ends[0]), [&](int i) {
          auto f = fma((Real)3, R.y(i, j, k), -Rs.y(i, j, k) + Rnu.y(i, j, k));
          v(i, j, k) = fma(alpha, f, v(i, j, k));
        });
        parallel_for(TeamVectorRange(team, ends[0]), [&](int i) {
          auto f = fma((Real)3, R.z(i, j, k), -Rs.z(i, j, k) + Rnu.z(i, j, k));
          w(i, j, k) = fma(alpha, f, w(i, j, k));
        });
      });
  } else {
    // Forward Euler for the first step
    auto alph = (Real)dt;
    parallel_for(
      "u_stepping_0", policy, KOKKOS_LAMBDA(member_t team) {
        int k = team.league_rank() / ends[1];
        int j = team.league_rank() % ends[1];
        parallel_for(TeamVectorRange(team, ends[0]), [&](int i) {
          u(i, j, k) = fma(R.x(i, j, k) + Rnu.x(i, j, k) / 2, alph, u(i, j, k));
          v(i, j, k) = fma(R.y(i, j, k) + Rnu.y(i, j, k) / 2, alph, v(i, j, k));
          w(i, j, k) = fma(R.z(i, j, k) + Rnu.z(i, j, k) / 2, alph, w(i, j, k));
        });
      });
  }
  policy.space().fence();

  ALPS_CHECK_LAST_DEVICE_ERROR();

  // Save the right-hand side
  std::swap(Ru_saved, rhs.Ru);

  solve_ueqn();
  // set boundary condition
  apply_bc(policy.space());

  ALPS_CHECK_LAST_DEVICE_ERROR();

  if (!steps_initialized) {
    steps_initialized = !steps_initialized;
  }

  Kokkos::Profiling::popRegion();
}

ChannelFlowSolverAB2CN::URhs ChannelFlowSolverAB2CN::calc_explicit_rhs() const
{
  auto const rhs_region = Kokkos::Profiling::ScopedRegion("U rhs");

  // Calculate viscous terms
  const Mesh& mesh = flow_field.mesh;

  Vector3Field<Real***, default_memory_pool> flux_u(
    Kokkos::view_alloc("U_flux", Kokkos::WithoutInitializing),
    mesh.extents(),
    {0, 0, 1});
  Vector3Field<Real***, default_memory_pool> flux_v(
    Kokkos::view_alloc("V_flux", Kokkos::WithoutInitializing),
    mesh.extents(),
    {0, 0, 1});
  Vector3Field<Real***, default_memory_pool> flux_w(
    Kokkos::view_alloc("W_flux", Kokkos::WithoutInitializing),
    mesh.extents(),
    {0, 0, 1});
  Kokkos::deep_copy(Kokkos::DefaultExecutionSpace(), flux_u, 0);
  Kokkos::deep_copy(Kokkos::DefaultExecutionSpace(), flux_v, 0);
  Kokkos::deep_copy(Kokkos::DefaultExecutionSpace(), flux_w, 0);
  Kokkos::DefaultExecutionSpace().fence();

  add_viscous_fluxes_from_u(Tensor33Field<Real***>(flux_u, flux_v, flux_w),
                            flow_field,
                            1 / Real(options.Re));
  ALPS_CHECK_LAST_DEVICE_ERROR();

  Vector3Field<Real***> fu_viscous;
  fu_viscous.z = div(flux_w, flow_field.mesh, NodePt(), "div(w_flux)");
  fu_viscous.y = div(flux_v, flow_field.mesh, CenterPt(), "div(v_flux)");
  fu_viscous.x = div(flux_u, flow_field.mesh, CenterPt(), "div(u_flux)");

  Kokkos::DefaultExecutionSpace().fence();
  ALPS_CHECK_LAST_DEVICE_ERROR();

  // Calculate other explicit terms
  Vector3Field<Real***> fu;

  calculate_convection_fluxes(Tensor33Field<Real***>(flux_u, flux_v, flux_w),
                              flow_field);
  fu.z = div(flux_w, flow_field.mesh, NodePt(), "div(w_flux)");
  fu.y = div(flux_v, flow_field.mesh, CenterPt(), "div(v_flux)");
  fu.x = div(flux_u, flow_field.mesh, CenterPt(), "div(u_flux)");
  for (const auto& body_force : body_forces) {
    body_force->add_forces(fu, Kokkos::DefaultExecutionSpace());
  }
  Kokkos::DefaultExecutionSpace().fence();
  ALPS_CHECK_LAST_DEVICE_ERROR();

  return {fu, fu_viscous};
}

} // namespace alps::solver
