//
// Created by xuananqing on 6/14/24.
//

#pragma once

#include <common/container/view_utils.h>
#include <common/kokkos_abstraction/exec_policy.h>
#include <solvers/field/bc_types.h>
#include <solvers/field/traits.h>

#include <Kokkos_Core.hpp>

#include <type_traits>

namespace alps::solver {

template<class FieldType>
void apply_bottom_bc(
  FieldType const&                     flow,
  NoSlipWall const&                    bc,
  Kokkos::DefaultExecutionSpace const& stream,
  std::enable_if_t<is_curvilinear_field<FieldType>::value>* = nullptr)
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
      u(i, j) = u_ * invJ(i, j);
      v(i, j) = v_ * invJ(i, j);
      w(i, j) = w_ * invJ(i, j);
    });
}

template<class FieldType>
void apply_bottom_bc(
  const FieldType&                     flow,
  const GradientWall&                  bc,
  const Kokkos::DefaultExecutionSpace& stream,
  std::enable_if_t<is_curvilinear_field<FieldType>::value>* = nullptr)
{
  if (!flow.mesh.comm().is_first(2)) return;

  const auto& u    = flow.u.x;
  const auto& v    = flow.u.y;
  const auto& w    = flow.u.z;
  const auto& dz   = flow.mesh.dz;
  const auto& invJ = flow.mesh.invJ;
  auto        ends = local_ends(u);

  /* See the similar function in src/solvers/ns/bc.cpp for the derivation of
   * the coefficients. */
  auto grad_1 = bc.grad_1;
  auto grad_2 = bc.grad_2;
  parallel_for(
    "set GradientWall bottom bc",
    LoopPolicy<2>(stream, {0, 0}, {ends[0], ends[1]}),
    KOKKOS_LAMBDA(int i, int j) {
      auto alpha  = dz(0) + dz(1);
      auto beta   = alpha + dz(0);
      auto coeff1 = alpha * alpha / dz(1) / beta;
      auto coeff2 = -dz(0) * dz(0) / dz(1) / beta;
      auto coeff0 = -dz(0) * alpha / beta;
      // Because the function is invoked after the velocity is multiplied by
      // J^{-1}, the d/dz gradient value is multiplied by J^{-2}
      u(i, j, 0) = coeff1 * u(i, j, 1) + coeff2 * u(i, j, 2)
                 + coeff0 * grad_1 * invJ(i, j) * invJ(i, j);
      v(i, j, 0) = coeff1 * v(i, j, 1) + coeff2 * v(i, j, 2)
                 + coeff0 * grad_2 * invJ(i, j) * invJ(i, j);
      w(i, j, 0) = 0;
    });
}

} // namespace alps::solver
