#include "solver.h"

#include <common/base/logging.h>
#include <common/container/view_utils.h>
#include <common/kokkos_abstraction/exec_policy.h>
#include <common/runtime/async_utils.h>
#include <decomp/mdcomm.h>
#include <solvers/ns/curvilinear_common/apply_bc.h>

namespace alps {
namespace solver {

void apply_top_bc(const FlowOverWaveField&             flow,
                  const GradientWall&                  bc,
                  const Kokkos::DefaultExecutionSpace& stream)
{
  if (!flow.mesh.comm().is_last(2)) return;

  const auto& u    = flow.u.x;
  const auto& v    = flow.u.y;
  const auto& w    = flow.u.z;
  const auto& dz   = flow.mesh.dz;
  const auto& invJ = flow.mesh.invJ;
  auto        ends = local_ends(u);

  /* See the similar function in src/solvers/ns/bc.cpp for the derivation of
   * the coefficients. */
  auto nz     = ends[2];
  auto grad_1 = bc.grad_1;
  auto grad_2 = bc.grad_2;
  parallel_for(
    "set GradientWall top bc",
    LoopPolicy<2>(stream, {0, 0}, {ends[0], ends[1]}),
    KOKKOS_LAMBDA(int i, int j) {
      auto alpha  = dz(nz - 2) + dz(nz - 3);
      auto beta   = alpha + dz(nz - 2);
      auto coeff1 = alpha * alpha / dz(nz - 3) / beta;
      auto coeff2 = -dz(nz - 2) * dz(nz - 2) / dz(nz - 3) / beta;
      auto coeff0 = dz(nz - 2) * alpha / beta;
      // Because the function is invoked after the velocity is multiplied by
      // J^{-1}, the d/dz gradient value is multiplied by J^{-2}
      u(i, j, nz - 1) = coeff1 * u(i, j, nz - 2) + coeff2 * u(i, j, nz - 3)
                      + coeff0 * grad_1 * invJ(i, j) * invJ(i, j);
      v(i, j, nz - 1) = coeff1 * v(i, j, nz - 2) + coeff2 * v(i, j, nz - 3)
                      + coeff0 * grad_2 * invJ(i, j) * invJ(i, j);
      w(i, j, nz - 2) = 0;
    });
}

void apply_top_bc(const FlowOverWaveSolver&            solver,
                  const TangentialStressWall&          bc,
                  const Kokkos::DefaultExecutionSpace& stream)
{
  auto const& flow = solver.flow_field;

  if (!flow.mesh.comm().is_last(2)) return;

  using real_t = decltype(GradientWall::grad_1);
  // For resolved boundary layer, compute the velocity gradient at the wall
  // using stress and viscosity
  GradientWall grad_bc;
  grad_bc.grad_1 = static_cast<real_t>((double)bc.tau_1 * solver.options.Re);
  grad_bc.grad_2 = static_cast<real_t>((double)bc.tau_2 * solver.options.Re);
  apply_top_bc(flow, grad_bc, stream);
  stream.fence(); // Wait for kernel before grad_bc is out of scope
}

void apply_top_bc(const FlowOverWaveField&             flow,
                  const NoSlipWall&                    bc,
                  const Kokkos::DefaultExecutionSpace& stream)
{
  if (!flow.mesh.comm().is_last(2)) return;

  auto        ends = local_ends(flow.u.x);
  const auto  u    = subview(flow.u.x, ALL, ALL, ends[2] - 1);
  const auto  v    = subview(flow.u.y, ALL, ALL, ends[2] - 1);
  const auto  w    = subview(flow.u.z, ALL, ALL, ends[2] - 2);
  const auto& invJ = flow.mesh.invJ;

  auto u_ = bc.u_;
  auto v_ = bc.v_;
  auto w_ = bc.w_;
  parallel_for(
    "set noslipwall top bc",
    LoopPolicy<2>(stream, {0, 0}, {ends[0], ends[1]}),
    KOKKOS_LAMBDA(int i, int j) {
      u(i, j) = u_ * invJ(i, j);
      v(i, j) = v_ * invJ(i, j);
      w(i, j) = w_ * invJ(i, j);
    });
}

void apply_bottom_bc(const FlowOverWaveField&             flow,
                     const NoSlipWallVarying&             bc,
                     const Kokkos::DefaultExecutionSpace& stream)
{
  if (!flow.mesh.comm().is_first(2)) return;

  auto        ends = local_ends(flow.u.x);
  const auto  u    = subview(flow.u.x, ALL, ALL, 0);
  const auto  v    = subview(flow.u.y, ALL, ALL, 0);
  const auto  w    = subview(flow.u.z, ALL, ALL, 0);
  const auto& invJ = flow.mesh.invJ;

  auto u_ = bc.u_;
  auto v_ = bc.v_;
  auto w_ = bc.w_;
  parallel_for(
    "set noslipwall bottom bc",
    LoopPolicy<2>(stream, {0, 0}, {ends[0], ends[1]}),
    KOKKOS_LAMBDA(int i, int j) {
      u(i, j) = u_(i, j) * invJ(i, j);
      v(i, j) = v_(i, j) * invJ(i, j);
      w(i, j) = w_(i, j) * invJ(i, j);
    });
}

void FlowOverWaveSolver::apply_bc(const Kokkos::DefaultExecutionSpace& space,
                                  const WhichBoundary boundary) const
{
  auto stream1 = get_next_stream();

  auto const& flow = flow_field;
  if (boundary == WhichBoundary::TopBC || boundary == WhichBoundary::Both) {
    auto event = get_device_event();
    enqueue(event, space);
    wait_for(event, stream1);

    bool bc_set = false;
    if (auto const* bc = dynamic_cast<NoSlipWall const*>(flow.top_bc.get());
        bc != nullptr) {
      apply_top_bc(flow, *bc, stream1);
      bc_set = true;
    }
    if (auto const* bc = dynamic_cast<GradientWall const*>(flow.top_bc.get());
        bc != nullptr) {
      apply_top_bc(flow, *bc, stream1);
      bc_set = true;
    }
    if (auto const* bc =
          dynamic_cast<TangentialStressWall const*>(flow.top_bc.get());
        bc != nullptr) {
      apply_top_bc(*this, *bc, stream1);
      bc_set = true;
    }
    if (!bc_set) {
      logger->warn("No top boundary condition set.");
    }
  }

  if (boundary == WhichBoundary::BottomBC || boundary == WhichBoundary::Both) {
    bool bc_set = false;
    if (auto const* bc = dynamic_cast<NoSlipWall const*>(flow.bottom_bc.get());
        bc != nullptr) {
      apply_bottom_bc(flow, *bc, space);
      bc_set = true;
    }
    if (auto const* bc =
          dynamic_cast<GradientWall const*>(flow.bottom_bc.get());
        bc != nullptr) {
      apply_bottom_bc(flow, *bc, space);
      bc_set = true;
    }
    if (auto const* bc =
          dynamic_cast<NoSlipWallVarying const*>(flow.bottom_bc.get());
        bc != nullptr) {
      apply_bottom_bc(flow, *bc, space);
      bc_set = true;
    }

    if (!bc_set) {
      logger->warn("No bottom boundary condition set.");
    }
  }

  space.fence();
  stream1.fence();
}

void FlowOverWaveSolver::set_boundary_stress_flux(
  Tensor33Field<Real***> const& fluxes,
  WhichBoundary                 boundary) const
{
  auto const& stream1 = get_next_stream();

  auto const& flow = flow_field;
  auto const& mesh = flow.mesh;
  using MT         = std::decay_t<decltype(mesh)>;

  auto const is_top    = mesh.comm().is_last(2);
  auto const is_bottom = mesh.comm().is_first(2);

  if (is_top
      && (boundary == WhichBoundary::TopBC
          || boundary == WhichBoundary::Both)) {
    auto const* bc =
      dynamic_cast<TangentialStressWall const*>(flow.top_bc.get());
    if (bc == nullptr) return;

    auto const  ends       = local_ends(fluxes.xz);
    auto const& tau_13_top = subview(fluxes.xz, ALL, ALL, ends[2] - 2).view();
    auto const& tau_23_top = subview(fluxes.yz, ALL, ALL, ends[2] - 2).view();
    auto const  tau_1      = bc->tau_1;
    auto const  tau_2      = bc->tau_2;

    Kokkos::parallel_for(
      LoopPolicy<2>(stream1, {0, 0}, {ends[0], ends[1]}),
      KOKKOS_LAMBDA(int i, int j) {
        tau_13_top(i, j) = tau_1 * (Real)MT::invJ_zeta_z();
        tau_23_top(i, j) = tau_2 * (Real)MT::invJ_zeta_z();
      });
  }

  if (is_bottom
      && (boundary == WhichBoundary::BottomBC
          || boundary == WhichBoundary::Both)) {
    (void)boundary; // suppress unused variable warning
  }

  stream1.fence();
}

} // namespace solver
} // namespace alps
