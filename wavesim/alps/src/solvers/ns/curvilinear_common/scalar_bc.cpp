//
// Created by xuanx004 on 7/13/24.
//

#include "scalar_bc.h"

#include <common/kokkos_abstraction/exec_policy.h>
#include <decomp/mdcomm.h>
#include <solvers/mesh/curvilinear_mesh.h>

namespace alps::solver {

void apply_top_bc(HaloView<Real***> const&             f,
                  ConstantDirichletBC const&           bc,
                  CurvilinearMesh const&               mesh,
                  Kokkos::DefaultExecutionSpace const& stream)
{
  if (!mesh.comm().is_last(2)) return;

  auto const  nz   = mesh.extent(2);
  auto const  val  = bc.value_;
  auto const& invJ = mesh.invJ;
  parallel_for(
    LoopPolicy<2>(stream, {0, 0}, {mesh.extent(0), mesh.extent(1)}),
    KOKKOS_LAMBDA(int i, int j) { f(i, j, nz - 1) = val * invJ(i, j); });
}

void apply_bottom_bc(HaloView<Real***> const&             f,
                     ConstantDirichletBC const&           bc,
                     CurvilinearMesh const&               mesh,
                     Kokkos::DefaultExecutionSpace const& stream)
{
  if (!mesh.comm().is_first(2)) return;

  auto const  val  = bc.value_;
  auto const& invJ = mesh.invJ;
  parallel_for(
    LoopPolicy<2>(stream, {0, 0}, {mesh.extent(0), mesh.extent(1)}),
    KOKKOS_LAMBDA(int i, int j) { f(i, j, 0) = val * invJ(i, j); });
}

void apply_top_bc(HaloView<Real***> const&             f,
                  ConstantGradientBC const&            bc,
                  CurvilinearMesh const&               mesh,
                  const Kokkos::DefaultExecutionSpace& stream)
{
  if (!mesh.comm().is_last(2)) return;

  auto const  nz   = mesh.extent(2);
  auto const  grad = bc.grad_;
  auto const& dz   = mesh.dz;
  auto const& invJ = mesh.invJ;
  parallel_for(
    LoopPolicy<2>(stream, {0, 0}, {mesh.extent(0), mesh.extent(1)}),
    KOKKOS_LAMBDA(int i, int j) {
      auto alpha      = dz(nz - 2) + dz(nz - 3);
      auto beta       = alpha + dz(nz - 2);
      auto coeff1     = alpha * alpha / dz(nz - 3) / beta;
      auto coeff2     = -dz(nz - 2) * dz(nz - 2) / dz(nz - 3) / beta;
      auto coeff0     = dz(nz - 2) * alpha / beta;
      f(i, j, nz - 1) = coeff1 * f(i, j, nz - 2) + coeff2 * f(i, j, nz - 3)
                      + coeff0 * grad * invJ(i, j) * invJ(i, j);
    });
}

void apply_bottom_bc(HaloView<Real***> const&             f,
                     ConstantGradientBC const&            bc,
                     CurvilinearMesh const&               mesh,
                     const Kokkos::DefaultExecutionSpace& stream)
{
  if (!mesh.comm().is_first(2)) return;

  const auto& dz   = mesh.dz;
  const auto& invJ = mesh.invJ;

  /* See the similar function in src/solvers/ns/bc.cpp for the derivation of
   * the coefficients. */
  auto const grad = bc.grad_;
  parallel_for(
    "set GradientWall bottom bc",
    LoopPolicy<2>(stream, {0, 0}, {mesh.extent(0), mesh.extent(1)}),
    KOKKOS_LAMBDA(int i, int j) {
      auto alpha  = dz(0) + dz(1);
      auto beta   = alpha + dz(0);
      auto coeff1 = alpha * alpha / dz(1) / beta;
      auto coeff2 = -dz(0) * dz(0) / dz(1) / beta;
      auto coeff0 = -dz(0) * alpha / beta;
      // Because the function is invoked after the velocity is multiplied by
      // J^{-1}, the d/dz gradient value is multiplied by J^{-2}
      f(i, j, 0) = coeff1 * f(i, j, 1) + coeff2 * f(i, j, 2)
                 + coeff0 * grad * invJ(i, j) * invJ(i, j);
    });
}

void apply_top_bc(HaloView<Real***> const&             f,
                  ConstantFluxBC const&                bc,
                  Real                                 D,
                  CurvilinearMesh const&               mesh,
                  Kokkos::DefaultExecutionSpace const& stream)
{
  if (!mesh.comm().is_last(2)) return;

  ConstantGradientBC grad_bc;
  grad_bc.grad_ = -bc.flux_ / D;
  apply_top_bc(f, grad_bc, mesh, stream);
}

void apply_bottom_bc(HaloView<Real***> const&             f,
                     ConstantFluxBC const&                bc,
                     Real                                 D,
                     CurvilinearMesh const&               mesh,
                     Kokkos::DefaultExecutionSpace const& stream)
{
  if (!mesh.comm().is_first(2)) return;

  ConstantGradientBC grad_bc;
  grad_bc.grad_ = -bc.flux_ / D;
  apply_bottom_bc(f, grad_bc, mesh, stream);
}

} // namespace alps::solver
