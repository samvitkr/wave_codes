//
// Created by xuanx004 on 5/27/24.
//

#pragma once

#include <common/container/view_types.h>
#include <common/real_type.h>

#include <Kokkos_Core_fwd.hpp>

#include <variant>

namespace alps::solver::hos {

struct RK4
{
  static constexpr char const* name = "Runge-Kutta 4 stage, 4th order";
}; // classic Runge-Kutta 4th order

struct LowStorageRK2N
{
  LowStorageRK2N(int                        stages,
                 std::vector<double> const& alpha_,
                 std::vector<double> const& beta_,
                 std::vector<double> const& c_)
    : S{stages}
    , alpha(alpha_)
    , beta(beta_)
    , c(c_)
  {}

  int S{};

  std::vector<double> alpha;
  std::vector<double> beta;
  std::vector<double> c;
}; // Williamson's 2N-storage Runge-Kutta method

struct LowStorageRK2C
{
  LowStorageRK2C(int                        stages,
                 std::vector<double> const& gamma_,
                 std::vector<double> const& b_,
                 std::vector<double> const& c_)
    : S{stages}
    , gamma(gamma_)
    , b(b_)
    , c(c_)
  {}

  int S{};

  std::vector<double> gamma;
  std::vector<double> b;
  std::vector<double> c;
};

struct ORK256 : LowStorageRK2N
{
  using base_t = LowStorageRK2N;

  static constexpr char const* name =
    "Runge-Kutta 5-stage, 2nd order optimized for wave propagation";

  ORK256();
}; // 5-stage, 2nd order low-storage optimized Runge-Kutta method

struct RK46NL : LowStorageRK2N
{
  using base_t = LowStorageRK2N;

  static constexpr char const* name =
    "Runge-Kutta 6-stage, 4th order, low-dissipation, low-dispersion";

  RK46NL();
}; // 6-stage, 4th order low-storage optimized Runge-Kutta method

struct TSLDDRK74 : LowStorageRK2C
{
  using base_t = LowStorageRK2C;

  static constexpr char const* name =
    "Runge-Kutta 7-stage, 4th order, low-dissipation, low-dispersion";

  TSLDDRK74();
}; // 7-stage, 4th order low-storage optimized Runge-Kutta method

using RKIntegrator = std::variant<RK4, ORK256, RK46NL, TSLDDRK74>;

RKIntegrator make_integrator(std::string const& name);

void update_intermediate_solution_rk4_async(
  MDView<Real** [2]>                   sol,
  MDView<Real const** [2]>             sol0,
  MDView<Real const** [2]>             dFdt,
  Real                                 dt,
  Kokkos::DefaultExecutionSpace const& space);

void update_final_solution_rk4_async(
  MDView<Real** [2]>                   sol,
  MDView<Real const** [2][4]>          dFdt,
  Real                                 dt,
  Kokkos::DefaultExecutionSpace const& space);

void update_solution_williamson_async(
  MDView<Real** [2]>                   U,
  MDView<Real** [2]>                   W,
  MDView<Real const** [2]>             dFdt,
  LowStorageRK2N                       integrator,
  int                                  stage,
  Real                                 dt,
  Kokkos::DefaultExecutionSpace const& space);

void update_solution_RK2C_async(MDView<Real** [2]>                   U,
                                MDView<Real** [2]>                   W,
                                MDView<Real const** [2]>             dFdt,
                                LowStorageRK2C                       integrator,
                                int                                  stage,
                                Real                                 dt,
                                Kokkos::DefaultExecutionSpace const& space);

} // namespace alps::solver::hos
