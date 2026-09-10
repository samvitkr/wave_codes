//
// Created by xuanx004 on 10/4/22.
//

#include "solver.h"

#include <common/base/macros.h>
#include <common/container/view_utils.h>
#include <common/device/devices.h>
#include <common/kokkos_abstraction/exec_policy.h>
#include <common/runtime/async_utils.h>
#include <spectral/spectral.h>

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
} // namespace

void FlowOverWaveSolver::calc_uhat_ab2(
  Vector3Field<Real***>                Ru,
  Kokkos::DefaultExecutionSpace const& stream)
{
  using Kokkos::parallel_for;
  const auto& field = flow_field;
  const auto& u     = field.u.x;
  const auto& v     = field.u.y;
  const auto& w     = field.u.z;
  const auto& invJ  = field.mesh.invJ;

  const auto ends = field.mesh.extents();

  /* Advance the velocity field
   * U = J(n+1)^{-1}\hat{u} = J(n)^{-1}u(n) + Δt*R(n+1/2)
   */
  GridPolicy<> const policy(stream, ends[1] * ends[2], Kokkos::AUTO);
  using member_t = decltype(policy)::member_type;
  if (steps_initialized) {
    // Second-order Adam Bashforth
    const auto  alpha = Real(-dt / 2);
    const auto  beta  = Real(1.5 * dt);
    const auto& Ru_s  = this->Ru_saved;
    parallel_for(
      "u_stepping", policy, KOKKOS_LAMBDA(member_t team) {
        using Kokkos::fma;
        int k = team.league_rank() / ends[1];
        int j = team.league_rank() % ends[1];
        parallel_for(
          Kokkos::TeamVectorRange(team, ends[0]), KOKKOS_TR_LAMBDA(int& i) {
            u(i, j, k) =
              fma(beta,
                  Ru.x(i, j, k),
                  fma(alpha, Ru_s.x(i, j, k), u(i, j, k) * invJ(i, j)));
          });
        parallel_for(
          Kokkos::TeamVectorRange(team, ends[0]), KOKKOS_TR_LAMBDA(int& i) {
            v(i, j, k) =
              fma(beta,
                  Ru.y(i, j, k),
                  fma(alpha, Ru_s.y(i, j, k), v(i, j, k) * invJ(i, j)));
          });
        parallel_for(
          Kokkos::TeamVectorRange(team, ends[0]), KOKKOS_TR_LAMBDA(int& i) {
            w(i, j, k) =
              fma(beta,
                  Ru.z(i, j, k),
                  fma(alpha, Ru_s.z(i, j, k), w(i, j, k) * invJ(i, j)));
          });
      });
  } else {
    // Forward Euler for the first step
    const auto delta_t = (Real)dt;
    parallel_for(
      "u_stepping_0", policy, KOKKOS_LAMBDA(member_t team) {
        using Kokkos::fma;
        int k = team.league_rank() / ends[1];
        int j = team.league_rank() % ends[1];
        parallel_for(
          Kokkos::TeamVectorRange(team, ends[0]), KOKKOS_TR_LAMBDA(int& i) {
            u(i, j, k) = fma(Ru.x(i, j, k), delta_t, u(i, j, k) * invJ(i, j));
            v(i, j, k) = fma(Ru.y(i, j, k), delta_t, v(i, j, k) * invJ(i, j));
            w(i, j, k) = fma(Ru.z(i, j, k), delta_t, w(i, j, k) * invJ(i, j));
          });
      });
  }
  policy.space().fence();

  ALPS_CHECK_LAST_DEVICE_ERROR();

  // Save the right-hand side
  std::swap(Ru_saved, Ru);

  if (field.bottom_bc->is_time_dependent) {
    auto new_bc     = field.bottom_bc->get_updated_bc(get_time());
    field.bottom_bc = std::move(new_bc);
    update_mesh_coefficients_from_bc(field, *field.bottom_bc, policy.space());
  }

  // set boundary condition
  apply_bc(policy.space());

  ALPS_CHECK_LAST_DEVICE_ERROR();
}

void FlowOverWaveSolver::calc_uhat(Vector3Field<Real***> Ru)
{
  // Update the scalar fields J^{-1}c
  advance_scalars_ab2();

  // cleanup SGS model if needed
  if (auto* model = std::get_if<DynamicSmagorinsky>(&turbulence_model);
      model != nullptr) {
    model->cleanup_Smag();
  }
  if (auto* model =
        std::get_if<AnisotropicMinimumDissipation>(&turbulence_model);
      model != nullptr) {
    model->cleanup_grad();
  }

  using Kokkos::parallel_for;
  const auto& field = flow_field;
  const auto& u     = field.u.x;
  const auto& v     = field.u.y;
  const auto& w     = field.u.z;

  Kokkos::Profiling::pushRegion("U stepping");

  auto stream = get_next_stream();
  if (options.integrator == "ab2cn") {
    calc_uhat_ab2cn(Ru, stream);
  } else if (options.integrator == "ab2") {
    // AB2 without Crank-Nicolson treatment of viscous term
    calc_uhat_ab2(Ru, stream);
  } else {
    throw std::runtime_error("Unsupported integrator type: "
                             + options.integrator);
  }

  spectral::dealias(create_inner_view(u).view(), field.mesh.grid, stream);
  spectral::dealias(create_inner_view(v).view(), field.mesh.grid, stream);
  spectral::dealias(create_inner_view(w).view(), field.mesh.grid, stream);

  if (!steps_initialized) {
    steps_initialized = !steps_initialized;
  }
  stream.fence();

  Kokkos::Profiling::popRegion();
}

} // namespace alps::solver
