//
// Created by xuanx004 on 5/27/24.
//

#include "solver_rk_stage.h"

#include <common/device/device_traits.h>
#include <common/kokkos_abstraction/exec_policy.h>
#include <common/utils/to_lower_case.h>

namespace alps::solver::hos {

namespace {
constexpr auto default_tile_size = []() -> Kokkos::Array<size_t, 3> {
  using Device = Kokkos::DefaultExecutionSpace;
  if constexpr (is_cuda_execution_space_v<Device>) return {32, 4, 2};
  if constexpr (is_hip_execution_space_v<Device>) return {64, 4, 1};
  return {};
}();
} // namespace

/// Calculate a forward Euler: sol = sol0 + dt * dFdt
void update_intermediate_solution_rk4_async(
  MDView<Real** [2]>                   sol,
  MDView<Real const** [2]>             sol0,
  MDView<Real const** [2]>             dFdt,
  Real                                 dt,
  Kokkos::DefaultExecutionSpace const& space)
{
  LoopPolicy<3> policy(space,
                       {0, 0, 0},
                       {sol.extent(0), sol.extent(1), sol.extent(2)},
                       default_tile_size);
  Kokkos::parallel_for(
    "update", policy, KOKKOS_LAMBDA(int i, int j, int k) {
      sol(i, j, k) = sol0(i, j, k) + dFdt(i, j, k) * dt;
    });
}

/// Calculate final step of RK4: F(t+dt) = F(t) + (k1 + 2k2 + 2k3 + k4)/6*dt
void update_final_solution_rk4_async(MDView<Real** [2]>                   sol,
                                     MDView<Real const** [2][4]>          dFdt,
                                     Real                                 dt,
                                     Kokkos::DefaultExecutionSpace const& space)
{
  LoopPolicy<3> policy(space,
                       {0, 0, 0},
                       {sol.extent(0), sol.extent(1), sol.extent(2)},
                       default_tile_size);

  Kokkos::parallel_for(
    "update", policy, KOKKOS_LAMBDA(int i, int j, int k) {
      sol(i, j, k) += ((dFdt(i, j, k, 0) + dFdt(i, j, k, 3))
                       + (dFdt(i, j, k, 1) + dFdt(i, j, k, 2)) * 2)
                    * (dt / 6);
    });
}

ORK256::ORK256()
  : base_t(5,
           {0.0, -1.0, -1.55798, -1.0, -0.45031},
           {0.2, 0.83204, 0.6, 0.35394, 0.2},
           {0.0, 0.2, 0.2, 0.8, 0.8})
{}

RK46NL::RK46NL()
  : base_t(6,
           {0.0,
            -0.737101392796,
            -1.634740794343,
            -0.744739003780,
            -1.469897351522,
            -2.813971388035},
           {0.032918605146,
            0.823256998200,
            0.381530948900,
            0.200092213184,
            1.718581042715,
            0.27},
           {0.0,
            0.032918605146,
            0.249351723343,
            0.466911705055,
            0.582030414044,
            0.847252983783})
{}

TSLDDRK74::TSLDDRK74()
  : base_t(7,
           {0.241566650129646868,
            0.0423866513027719953,
            0.215602732678803776,
            0.232328007537583987,
            0.256223412574146438,
            0.0978694102142697230,
            0.0},
           {0.0941840925477795334,
            0.149683694803496998,
            0.285204742060440058,
            -0.122201846148053668,
            0.0605151571191401122,
            0.345986987898399296,
            0.186627171718797670},
           {0.0,
            0.335750742677426401,
            0.286254438654048527,
            0.744675262090520366,
            0.639198690801246909,
            0.723609252956949472,
            0.91124223849547205})
{}

/// Update the solution following the Williamson's 2N-storage RK
/// Wi = alpha_i * W_(i-1) + delta_t * F(U_(i-1), t_i)
/// U_i = U_(i-1) + beta_i * W_i
void update_solution_williamson_async(
  MDView<Real** [2]>                   U,
  MDView<Real** [2]>                   W,
  MDView<Real const** [2]>             dFdt,
  LowStorageRK2N                       integrator,
  int                                  stage,
  Real                                 dt,
  Kokkos::DefaultExecutionSpace const& space)
{
  LoopPolicy<3> policy(space,
                       {0, 0, 0},
                       {U.extent(0), U.extent(1), U.extent(2)},
                       default_tile_size);

  auto alpha = static_cast<Real>(integrator.alpha.at(stage));
  auto beta  = static_cast<Real>(integrator.beta.at(stage));

  if (alpha == 0) {
    Kokkos::parallel_for(
      "update", policy, KOKKOS_LAMBDA(int i, int j, int k) {
        W(i, j, k) = dt * dFdt(i, j, k);
        U(i, j, k) += beta * W(i, j, k);
      });
  } else {
    Kokkos::parallel_for(
      "update", policy, KOKKOS_LAMBDA(int i, int j, int k) {
        W(i, j, k) = alpha * W(i, j, k) + dt * dFdt(i, j, k);
        U(i, j, k) += beta * W(i, j, k);
      });
  }
}

void update_solution_RK2C_async(MDView<Real** [2]>                   U,
                                MDView<Real** [2]>                   W,
                                MDView<Real const** [2]>             dFdt,
                                LowStorageRK2C                       integrator,
                                int                                  stage,
                                Real                                 dt,
                                Kokkos::DefaultExecutionSpace const& space)
{
  LoopPolicy<3> policy(space,
                       {0, 0, 0},
                       {U.extent(0), U.extent(1), U.extent(2)},
                       default_tile_size);

  auto gamma = static_cast<Real>(integrator.gamma.at(stage));
  auto b     = static_cast<Real>(integrator.b.at(stage));

  if (stage == 0) {
    Kokkos::parallel_for(
      "update", policy, KOKKOS_LAMBDA(int i, int j, int k) {
        auto u0    = U(i, j, k);
        auto w     = u0 + dt * b * dFdt(i, j, k);
        U(i, j, k) = w + dt * gamma * dFdt(i, j, k);
        W(i, j, k) = w;
      });
  } else if (stage != integrator.S - 1) {
    Kokkos::parallel_for(
      "update", policy, KOKKOS_LAMBDA(int i, int j, int k) {
        W(i, j, k) += dt * b * dFdt(i, j, k);
        U(i, j, k) = W(i, j, k) + dt * gamma * dFdt(i, j, k);
      });
  } else {
    // last stage
    Kokkos::parallel_for(
      "update", policy, KOKKOS_LAMBDA(int i, int j, int k) {
        U(i, j, k) = W(i, j, k) + dt * b * dFdt(i, j, k);
      });
  }
}

RKIntegrator make_integrator(std::string const& name)
{
  auto const name_lower = to_lower_case(name);

  if (name_lower == "rk4") {
    return RK4{};
  }
  if (name_lower == "ork256") {
    return ORK256{};
  }
  if (name_lower == "rk46nl") {
    return RK46NL{};
  }
  if (name_lower == "tslddrk74") {
    return TSLDDRK74{};
  }
  throw std::runtime_error("Unknown integrator: " + name);
}

} // namespace alps::solver::hos
