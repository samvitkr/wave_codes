#pragma once

#include <common/real_type.h>
#include <solvers/mesh/curvilinear_mesh.h>
#include <solvers/ns/channel/diffusion_cn.h>

#include <type_traits>

namespace alps::solver {

template<typename VarLoc, typename MT>
struct DiffusionCNZetaCoeffProvider
{
  HaloView<Real const*> dz;
  HaloView<Real const*> dzw;
  HaloView<Real const*> zeta;
  MDView<Real const**>  exr;
  MDView<Real const**>  eyr;
  MDView<Real const**>  invJ;
  bool                  is_bottom;
  bool                  is_top;
  int                   nz;
  DiffusionCNBCType     bottom_bc;
  DiffusionCNBCType     top_bc;
  Real                  alpha;

  struct Coefficients
  {
    Real lower;
    Real diag;
    Real upper;
  };

  DiffusionCNZetaCoeffProvider(MT const&         mesh,
                               DiffusionCNBCType bottom_bc_,
                               DiffusionCNBCType top_bc_,
                               Real              alpha_)
    : dz{mesh.dz}
    , dzw{mesh.dzw}
    , zeta{std::is_same_v<VarLoc, CenterPt> ? mesh.zw : mesh.zz}
    , exr{mesh.exr}
    , eyr{mesh.eyr}
    , invJ{mesh.invJ}
    , is_bottom{mesh.grid.comm().is_first(2)}
    , is_top{mesh.grid.comm().is_last(2)}
    , nz{mesh.extent(2)}
    , bottom_bc{bottom_bc_}
    , top_bc{top_bc_}
    , alpha{alpha_}
  {}

  KOKKOS_FUNCTION Coefficients coeff(int i, int j, int k) const
  {
    return coeff(i, j, k, VarLoc{});
  }

  KOKKOS_INLINE_FUNCTION
  Coefficients coeff(int i, int j, int k, CenterPt /*VarLoc*/) const
  {
    if (is_bottom && k == 0 && bottom_bc == DiffusionCNBCType::Dirichlet) {
      return {0, 1, 0};
    }
    if (is_bottom && k == 1 && bottom_bc == DiffusionCNBCType::Dirichlet) {
      // uses a three-point stencil to compute the boundary gradient
      auto dz0  = dz(0);
      auto dz1  = dz(1);
      auto beta = dz0 / (dz0 + dz1);
      auto c0   = -(1 + beta) / dz0; // u(0) coeff for du/dz at the boundary
      auto c1   = 1 / beta / dz1;    // u(1) coeff for du/dz at the boundary
      auto c2   = -beta / dz1;       // u(2) coeff for du/dz at the boundary

      auto const g_l   = MT::g33(zeta(0), exr(i, j), eyr(i, j), invJ(i, j));
      auto const g_u   = MT::g33(zeta(1), exr(i, j), eyr(i, j), invJ(i, j));
      auto       lower = alpha * g_l * c0 / dzw(0);
      auto       diag  = Real(1) + alpha * (g_l * c1 + g_u / dz(1)) / dzw(0);
      auto       upper = alpha * (g_l * c2 - g_u / dz(1)) / dzw(0);
      return {lower, diag, upper};
    }
    if (is_bottom && k == 0 && bottom_bc == DiffusionCNBCType::Neumann) {
      // The apparent one-sided stencil is only a way to express the boundary
      // value implied by the Neumann condition. The first interior equation
      // uses the same boundary-flux stencil, so substituting this boundary
      // value collapses exactly to imposing the prescribed Neumann flux at the
      // boundary. Thus the low-order boundary-value reconstruction does not
      // modify the tridiagonal operator or the interior unknowns. The boundary
      // value is reset after the implicit step.
      return {0, -Real(1) / dz(0), Real(1) / dz(0)};
    }
    if (is_top && k == nz - 2 && top_bc == DiffusionCNBCType::Dirichlet) {
      // uses a three-point stencil to compute the boundary gradient
      auto dz0  = dz(nz - 2);
      auto dz1  = dz(nz - 3);
      auto beta = dz0 / (dz0 + dz1);
      auto c2   = (1 + beta) / dz0; // u(nz-1) coeff for du/dz at the boundary
      auto c1   = -1 / beta / dz1;  // u(nz-2) coeff for du/dz at the boundary
      auto c0   = beta / dz1;       // u(nz-3) coeff for du/dz at the boundary

      auto const g_l = MT::g33(zeta(nz - 3), exr(i, j), eyr(i, j), invJ(i, j));
      auto const g_u = MT::g33(zeta(nz - 2), exr(i, j), eyr(i, j), invJ(i, j));

      auto lower = -alpha * (g_u * c0 + g_l / dz1) / dzw(nz - 3);
      auto diag  = Real(1) + alpha * (g_l / dz1 - g_u * c1) / dzw(nz - 3);
      auto upper = -alpha * g_u * c2 / dzw(nz - 3);
      return {lower, diag, upper};
    }
    if (is_top && k == nz - 1 && top_bc == DiffusionCNBCType::Dirichlet) {
      return {0, 1, 0};
    }
    if (is_top && k == nz - 1 && top_bc == DiffusionCNBCType::Neumann) {
      return {-Real(1) / dz(k - 1), Real(1) / dz(k - 1), 0};
    }

    auto const g_l = MT::g33(zeta(k - 1), exr(i, j), eyr(i, j), invJ(i, j));
    auto const g_u = MT::g33(zeta(k), exr(i, j), eyr(i, j), invJ(i, j));

    auto lower = -alpha * g_l / (dz(k - 1) * dzw(k - 1));
    auto diag  = Real(1) + alpha * (g_l / dz(k - 1) + g_u / dz(k)) / dzw(k - 1);
    auto upper = -alpha * g_u / (dz(k) * dzw(k - 1));
    return {lower, diag, upper};
  }

  KOKKOS_INLINE_FUNCTION Coefficients coeff(int i,
                                            int j,
                                            int k,
                                            NodePt /*VarLoc*/) const
  {
    if (is_bottom && k == 0) {
      return {0, 1, 0};
    }
    if (is_top && k == nz - 2) {
      return {0, 1, 0};
    }

    auto const g_l = MT::g33(zeta(k), exr(i, j), eyr(i, j), invJ(i, j));
    auto const g_u = MT::g33(zeta(k + 1), exr(i, j), eyr(i, j), invJ(i, j));

    auto lower = -alpha * g_l / (dzw(k - 1) * dz(k));
    auto diag  = Real(1) + alpha * (g_l / dzw(k - 1) + g_u / dzw(k)) / dz(k);
    auto upper = -alpha * g_u / (dzw(k) * dz(k));
    return {lower, diag, upper};
  }
};

} // namespace alps::solver
