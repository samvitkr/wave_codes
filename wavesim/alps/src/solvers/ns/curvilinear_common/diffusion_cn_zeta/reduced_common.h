#pragma once

#include <common/base/macros.h>
#include <common/container/view_types.h>
#include <common/kokkos_abstraction/exec_policy.h>
#include <common/real_type.h>

#include <Kokkos_Core.hpp>
#include <mpipp/comm.h>

namespace alps::solver {

struct DistTridiagMpiTag
{
  static int constexpr Initial  = 20;
  static int constexpr Reduced  = 21;
  static int constexpr Backward = 22;
};

template<bool shift_c0_on_rank0, typename CoeffFunctor>
struct DistTridiagLocalElimination
{
  struct Forward
  {};

  struct Backward
  {};

  MDView<Real***>       dst;
  MDView<Real const***> rhs;
  CoeffFunctor          coeff;

  MDView<Real***> aa;
  MDView<Real***> cc;
  MDView<Real***> first_line;
  MDView<Real***> last_line;
  MDView<Real**>  first_sol;
  MDView<Real**>  next_sol;

  int  nx;
  int  n_eqns;
  bool is_bottom;
  bool is_top;

  DistTridiagLocalElimination(MDView<Real***> const&       dst_,
                              MDView<Real const***> const& rhs_,
                              CoeffFunctor const&          coeff_,
                              MDView<Real***> const&       aa_,
                              MDView<Real***> const&       cc_,
                              MDView<Real***> const&       first_line_,
                              MDView<Real***> const&       last_line_,
                              mpipp::communicator const&   comm)
    : dst{dst_}
    , rhs{rhs_}
    , coeff{coeff_}
    , aa{aa_}
    , cc{cc_}
    , first_line{first_line_}
    , last_line{last_line_}
    , first_sol{subview(first_line, Kokkos::ALL, Kokkos::ALL, 2)}
    // Reuse last_line's rhs slot (index 2) as the buffer for the next
    // rank's boundary solution; last_line is no longer needed after the
    // InitialCollapse exchange
    , next_sol{subview(last_line, Kokkos::ALL, Kokkos::ALL, 2)}
    , nx{aa_.extent_int(0)}
    , n_eqns{aa_.extent_int(2)}
    , is_bottom{comm.rank() == 0}
    , is_top{comm.rank() == comm.size() - 1}
  {}

  /**
   * Eliminates local interior unknowns along each tridiagonal line and stores
   * normalized reduced equations in first_line and last_line.
   *
   * For an interior rank r > 0, with F and L denoting the first and last local
   * unknowns, the stored equations are
   *
   *   a_r x_{r-1,L} + x_{r,F} + c_r x_{r+1,F} = d_r,
   *   A_r x_{r,F}   + x_{r,L} + C_r x_{r+1,F} = D_r.
   *
   * The coefficients a_r, c_r, d_r, A_r, C_r, D_r are reduced coefficients,
   * not the original tridiagonal coefficients.
   *
   * On rank 0, if shift_c0_on_rank0=true, the rank-0 interior unknowns are
   * also eliminated from the first equation, giving
   *
   *   x_{0,F} + c_0 x_{1,F} = d_0.
   *
   * If shift_c0_on_rank0=false, the first rank's first equation is only
   * normalized, not fully reduced across the local block:
   *
   *   x_{0,F} + c_0 x_{0,F+1} = d_0.
   */
  KOKKOS_FUNCTION void operator()(Forward /*tag*/,
                                  GridPolicy<>::member_type team) const
  {
    using Kokkos::TeamThreadRange;

    int const j = team.league_rank();
    Kokkos::parallel_for(
      TeamThreadRange(team, nx), KOKKOS_TR_LAMBDA(int& i) {
        auto const c0_coeff = coeff.coeff(i, j, 0);
        Real       b0       = c0_coeff.diag;
        Real       c0       = c0_coeff.upper;
        Real       d0       = rhs(i, j, 0);

        auto c1_coeff = coeff.coeff(i, j, 1);
        if (is_bottom) {
          Real factor  = c1_coeff.lower / b0;
          Real bb      = Real(1) / (c1_coeff.diag - factor * c0);
          dst(i, j, 1) = bb * (rhs(i, j, 1) - factor * d0);
          cc(i, j, 1)  = bb * c1_coeff.upper;
          aa(i, j, 1)  = 0;
          if constexpr (shift_c0_on_rank0) {
            d0 = d0 - c0 * dst(i, j, 1);
            c0 = -c0 * cc(i, j, 1);
          }
        } else {
          dst(i, j, 1) = rhs(i, j, 1) / c1_coeff.diag;
          cc(i, j, 1)  = c1_coeff.upper / c1_coeff.diag;
          aa(i, j, 1)  = c1_coeff.lower / c1_coeff.diag;
          b0           = b0 - c0 * aa(i, j, 1);
          d0           = d0 - c0 * dst(i, j, 1);
          if (is_top && n_eqns == 2) {
            c0 = 0;
          } else {
            c0 = -c0 * cc(i, j, 1);
          }
        }

        for (int k = 2; k < n_eqns; ++k) {
          auto coeffs = coeff.coeff(i, j, k);
          Real bb = Real(1) / (coeffs.diag - coeffs.lower * cc(i, j, k - 1));
          dst(i, j, k) = (rhs(i, j, k) - coeffs.lower * dst(i, j, k - 1)) * bb;
          if (k == n_eqns - 1 && is_top) {
            cc(i, j, k) = 0;
          } else {
            cc(i, j, k) = coeffs.upper * bb;
          }
          if (is_bottom) {
            aa(i, j, k) = 0;
            if constexpr (shift_c0_on_rank0) {
              d0 = d0 - c0 * dst(i, j, k);
              c0 = -c0 * cc(i, j, k);
            }
          } else {
            aa(i, j, k) = (-coeffs.lower * aa(i, j, k - 1)) * bb;
            b0          = b0 - c0 * aa(i, j, k);
            d0          = d0 - c0 * dst(i, j, k);
            c0          = -c0 * cc(i, j, k);
          }
        }

        aa(i, j, 0)  = is_bottom ? Real(0) : c0_coeff.lower / b0;
        cc(i, j, 0)  = c0 / b0;
        dst(i, j, 0) = d0 / b0;

        first_line(i, j, 0) = aa(i, j, 0);
        first_line(i, j, 1) = cc(i, j, 0);
        first_line(i, j, 2) = dst(i, j, 0);
        last_line(i, j, 0)  = aa(i, j, n_eqns - 1);
        last_line(i, j, 1)  = cc(i, j, n_eqns - 1);
        last_line(i, j, 2)  = dst(i, j, n_eqns - 1);
      });
  }

  /** Reconstructs local solutions from the reduced boundary solution. */
  KOKKOS_FUNCTION void operator()(Backward /*tag*/,
                                  GridPolicy<>::member_type team) const
  {
    int const j = team.league_rank();
    Kokkos::parallel_for(
      Kokkos::TeamThreadRange(team, nx), KOKKOS_TR_LAMBDA(int& i) {
        Real dd0   = first_sol(i, j);
        Real dd_p1 = is_top ? dst(i, j, n_eqns - 1)
                            : dst(i, j, n_eqns - 1)
                                - next_sol(i, j) * cc(i, j, n_eqns - 1);
        if (!is_bottom) {
          dd_p1 += -aa(i, j, n_eqns - 1) * dd0;
        }
        dst(i, j, n_eqns - 1) = dd_p1;
        for (int k = n_eqns - 2; k > 0; --k) {
          dd_p1 = dst(i, j, k) - cc(i, j, k) * dd_p1;
          if (!is_bottom) {
            dd_p1 += -aa(i, j, k) * dd0;
          }
          dst(i, j, k) = dd_p1;
        }
        if (is_bottom && !shift_c0_on_rank0) {
          dst(i, j, 0) = dd0 - cc(i, j, 0) * dd_p1;
        } else {
          dst(i, j, 0) = dd0;
        }
      });
  }
};

} // namespace alps::solver
