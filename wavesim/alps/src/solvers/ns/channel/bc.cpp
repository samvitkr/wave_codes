#include "bc.h"

#include <common/container/view_utils.h>
#include <common/kokkos_abstraction/exec_policy.h>
#include <decomp/mdcomm.h>

namespace alps {
namespace solver {

void apply_top_bc(const FlowField&                     flow,
                  const GradientWall&                  bc,
                  const Kokkos::DefaultExecutionSpace& stream)
{
  if (!flow.mesh.comm().is_last(2)) return;

  const auto& u    = flow.u.x;
  const auto& v    = flow.u.y;
  const auto& w    = flow.u.z;
  const auto& dz   = flow.mesh.dz;
  const auto  hbar = flow.mesh.hbar;
  auto        ends = local_ends(u);

  /* Compute the boundary velocity with prescribed gradient using a
   * three-point stencil. The coefficients are the same as those used at the
   * bottom except for coeff0, which is the negative of the bottom coeff0. */
  auto nz     = ends[2];
  auto grad_1 = bc.grad_1;
  auto grad_2 = bc.grad_2;
  parallel_for(
    LoopPolicy<2>(stream, {0, 0}, {ends[0], ends[1]}),
    KOKKOS_LAMBDA(int i, int j) {
      auto alpha  = dz(nz - 2) + dz(nz - 3);
      auto beta   = alpha + dz(nz - 2);
      auto coeff1 = alpha * alpha / dz(nz - 3) / beta;
      auto coeff2 = -dz(nz - 2) * dz(nz - 2) / dz(nz - 3) / beta;
      auto coeff0 = dz(nz - 2) * alpha / beta * hbar;
      u(i, j, nz - 1) =
        coeff1 * u(i, j, nz - 2) + coeff2 * u(i, j, nz - 3) + coeff0 * grad_1;
      v(i, j, nz - 1) =
        coeff1 * v(i, j, nz - 2) + coeff2 * v(i, j, nz - 3) + coeff0 * grad_2;
      w(i, j, nz - 2) = 0;
    });
}

void apply_bottom_bc(const FlowField&                     flow,
                     const GradientWall&                  bc,
                     const Kokkos::DefaultExecutionSpace& stream)
{
  if (!flow.mesh.comm().is_first(2)) return;

  const auto& u    = flow.u.x;
  const auto& v    = flow.u.y;
  const auto& w    = flow.u.z;
  const auto& dz   = flow.mesh.dz;
  const auto  hbar = flow.mesh.hbar;
  auto        ends = local_ends(u);

  /* Compute the boundary velocity with prescribed gradient using a
   * three-point stencil
   * u0 = -(dz0*(dz0 + dz1)/(2*dz0 + dz1)) * H * grad +
   *       (dz0 + dz1)^2/(dz1*(2*dz0 + dz1)) * u1-
   *       (dz0^2)/(dz1*(2*dz0 + dz1)) * u2
   */
  auto grad_1 = bc.grad_1;
  auto grad_2 = bc.grad_2;
  parallel_for(
    LoopPolicy<2>(stream, {0, 0}, {ends[0], ends[1]}),
    KOKKOS_LAMBDA(int i, int j) {
      auto alpha  = dz(0) + dz(1);
      auto beta   = alpha + dz(0);
      auto coeff1 = alpha * alpha / dz(1) / beta;
      auto coeff2 = -dz(0) * dz(0) / dz(1) / beta;
      auto coeff0 = -dz(0) * alpha / beta * hbar;
      u(i, j, 0)  = coeff1 * u(i, j, 1) + coeff2 * u(i, j, 2) + coeff0 * grad_1;
      v(i, j, 0)  = coeff1 * v(i, j, 1) + coeff2 * v(i, j, 2) + coeff0 * grad_2;
      w(i, j, 0)  = 0;
    });
}

void apply_top_bc(const FlowField&                     flow,
                  const NoSlipWall&                    bc,
                  const Kokkos::DefaultExecutionSpace& stream)
{
  if (!flow.mesh.comm().is_last(2)) return;

  auto       ends = local_ends(flow.u.x);
  const auto u    = subview(flow.u.x, ALL, ALL, ends[2] - 1);
  const auto v    = subview(flow.u.y, ALL, ALL, ends[2] - 1);
  const auto w    = subview(flow.u.z, ALL, ALL, ends[2] - 2);

  auto u_ = bc.u_;
  auto v_ = bc.v_;
  auto w_ = bc.w_;
  parallel_for(
    LoopPolicy<2>(stream, {0, 0}, {ends[0], ends[1]}),
    KOKKOS_LAMBDA(int i, int j) {
      u(i, j) = u_;
      v(i, j) = v_;
      w(i, j) = w_;
    });
}

void apply_bottom_bc(const FlowField&                     flow,
                     const NoSlipWall&                    bc,
                     const Kokkos::DefaultExecutionSpace& stream)
{
  if (!flow.mesh.comm().is_first(2)) return;

  auto       ends = local_ends(flow.u.x);
  const auto u    = subview(flow.u.x, ALL, ALL, 0);
  const auto v    = subview(flow.u.y, ALL, ALL, 0);
  const auto w    = subview(flow.u.z, ALL, ALL, 0);

  auto u_ = bc.u_;
  auto v_ = bc.v_;
  auto w_ = bc.w_;
  parallel_for(
    LoopPolicy<2>(stream, {0, 0}, {ends[0], ends[1]}),
    KOKKOS_LAMBDA(int i, int j) {
      u(i, j) = u_;
      v(i, j) = v_;
      w(i, j) = w_;
    });
}

} // namespace solver
} // namespace alps
