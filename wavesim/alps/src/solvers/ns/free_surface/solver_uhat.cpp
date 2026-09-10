//
// Created by xuanx004 on 10/4/22.
//

#include "solver.h"

#include <common/base/macros.h>
#include <common/container/view_utils.h>
#include <common/device/devices.h>
#include <common/kokkos_abstraction/exec_policy.h>
#include <common/kokkos_abstraction/pool_space.h>
#include <common/runtime/async_utils.h>
#include <solvers/ns/free_surface/bc.h>
#include <spectral/spectral.h>

namespace alps::solver {

namespace {

void advance_velocity_AB2(MDView<Real***> const&                       u,
                          MDView<Real***> const&                       v,
                          MDView<Real***> const&                       w,
                          decltype(FreeSurfaceSolver::Ru_saved) const& Ru,
                          decltype(FreeSurfaceSolver::Ru_saved) const& Ru_s,
                          MDView<Real const**> const&                  invJ0,
                          Real                                         dt,
                          bool                                 is_first_step,
                          Kokkos::DefaultExecutionSpace const& stream)
{
  auto const         ends = local_extents(Ru.x);
  GridPolicy<> const policy(stream, ends[1] * ends[2], Kokkos::AUTO);
  using member_t = decltype(policy)::member_type;
  if (is_first_step) {
    // Forward Euler for the first step
    parallel_for(
      "u_stepping_0", policy, KOKKOS_LAMBDA(member_t team) {
        using Kokkos::fma;
        int k = team.league_rank() / ends[1];
        int j = team.league_rank() % ends[1];
        Kokkos::parallel_for(
          Kokkos::TeamVectorRange(team, ends[0]), KOKKOS_TR_LAMBDA(int& i) {
            u(i, j, k)      = fma(Ru.x(i, j, k), dt, u(i, j, k) * invJ0(i, j));
            Ru_s.x(i, j, k) = u(i, j, k); // store J^{-1}u for 2nd stage RK
            v(i, j, k)      = fma(Ru.y(i, j, k), dt, v(i, j, k) * invJ0(i, j));
            Ru_s.y(i, j, k) = v(i, j, k); // store J^{-1}u for 2nd stage RK
            w(i, j, k)      = fma(Ru.z(i, j, k), dt, w(i, j, k) * invJ0(i, j));
            Ru_s.z(i, j, k) = w(i, j, k); // store J^{-1}u for 2nd stage RK
          });
      });
    return;
  }
  parallel_for(
    "u_stepping", policy, KOKKOS_LAMBDA(member_t team) {
      using Kokkos::fma;
      int k = team.league_rank() / ends[1];
      int j = team.league_rank() % ends[1];
      // Second-order Adam Bashforth
      Kokkos::parallel_for(
        Kokkos::TeamVectorRange(team, ends[0]), KOKKOS_TR_LAMBDA(int& i) {
          u(i, j, k)      = fma(fma((Real)3, Ru.x(i, j, k), -Ru_s.x(i, j, k)),
                           dt / 2,
                           u(i, j, k) * invJ0(i, j));
          Ru_s.x(i, j, k) = u(i, j, k); // store J^{-1}u for 2nd stage RK
        });
      Kokkos::parallel_for(
        Kokkos::TeamVectorRange(team, ends[0]), KOKKOS_TR_LAMBDA(int& i) {
          v(i, j, k)      = fma(fma((Real)3, Ru.y(i, j, k), -Ru_s.y(i, j, k)),
                           dt / 2,
                           v(i, j, k) * invJ0(i, j));
          Ru_s.y(i, j, k) = v(i, j, k); // store J^{-1}v for 2nd stage RK
        });
      Kokkos::parallel_for(
        Kokkos::TeamVectorRange(team, ends[0]), KOKKOS_TR_LAMBDA(int& i) {
          w(i, j, k)      = fma(fma((Real)3, Ru.z(i, j, k), -Ru_s.z(i, j, k)),
                           dt / 2,
                           w(i, j, k) * invJ0(i, j));
          Ru_s.z(i, j, k) = w(i, j, k); // store J^{-1}w for 2nd stage RK
        });
    });
}

} // anonymous namespace

void FreeSurfaceSolver::calc_uhat_rk2(Vector3Field<Real***> Ru, int const stage)
{
  if (stage != 1 && stage != 2) {
    throw std::runtime_error("Invalid RK stage");
  }

  if (stage == 1) {
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
  }

  const auto& field = flow_field;
  const auto& u_in  = create_inner_view(field.u.x).view();
  const auto& v_in  = create_inner_view(field.u.y).view();
  const auto& w_in  = create_inner_view(field.u.z).view();

  const auto is_top = field.mesh.comm().is_last(2);
  const auto ends   = field.mesh.extents();

  auto stream = get_next_stream();

  decltype(TopWaveMesh::invJ) invJ0;
  if (stage == 1) {
    // Save the inverse Jacobian before updating the surface elevation
    invJ0 =
      MDView<Real**, default_memory_pool>("invJ0", field.mesh.invJ.layout());
    Kokkos::deep_copy(stream, invJ0, field.mesh.invJ);
  }

  // Update the surface elevation
  field.top_bc->get_updated_eta(field.eta,
                                subview(u_in, ALL, ALL, ends[2] - 1),
                                subview(v_in, ALL, ALL, ends[2] - 1),
                                subview(w_in, ALL, ALL, ends[2] - 2),
                                1 / options.Re,
                                (Real)dt,
                                stage,
                                field.mesh);
  field.mesh.update_metric_coefficients(field.eta, stream);
  stream.fence();

  ALPS_CHECK_LAST_DEVICE_ERROR();

  MDView<Real** [3]> vec_uz;
  if (is_top) {
    vec_uz =
      MDView<Real** [3], default_memory_pool>("uz surface", ends[0], ends[1]);
  }
  field.top_bc->get_surface_uz(
    vec_uz, field.u, 1 / Real(options.Re), field.mesh);

  if (stage == 1) {
    if (steps_initialized) {
      advance_velocity_AB2(
        u_in, v_in, w_in, Ru, Ru_saved, invJ0, (Real)dt, false, stream);
    } else {
      if (Ru_saved.x.span() < 1) { // Initialize the saved right-hand side
        Ru_saved = Vector3Field<Real***, default_memory_pool>(
          Kokkos::view_alloc("Ru", Kokkos::WithoutInitializing),
          field.mesh.extents());
      }
      advance_velocity_AB2(
        u_in, v_in, w_in, Ru, Ru_saved, invJ0, (Real)dt, true, stream);
    }
  } else {
    // \hat{u} is unchanged in the 2nd stage and was stored in Ru_saved
    deep_copy(stream, u_in, create_inner_view(Ru_saved.x).view());
    deep_copy(stream, v_in, create_inner_view(Ru_saved.y).view());
    deep_copy(stream, w_in, create_inner_view(Ru_saved.z).view());
  }
  stream.fence();

  ALPS_CHECK_LAST_DEVICE_ERROR();

  if (stage == 2) {
    // Save the right-hand side
    std::swap(Ru_saved, Ru);
  }

  // set boundary condition
  apply_bottom_bc(stream);
  stream.fence();
  field.top_bc->get_surface_uhat(
    field.u, vec_uz, field.pp, (Real)dt, field.mesh);

  ALPS_CHECK_LAST_DEVICE_ERROR();

  dealias(u_in, field.mesh.grid, stream);
  dealias(v_in, field.mesh.grid, stream);
  dealias(w_in, field.mesh.grid, stream);
  stream.fence();

  if (stage == 2 && !steps_initialized) {
    steps_initialized = !steps_initialized;
  }
}

void FreeSurfaceSolver::calc_uhat(Vector3Field<Real***> Ru, int rk_stage)
{
  if (rk_stage == 1) {
    advance_scalars_ab2();

    Kokkos::Profiling::pushRegion("U stepping");
    calc_uhat_rk2(std::move(Ru), 1);
    Kokkos::Profiling::popRegion();
  } else if (rk_stage == 2) {
    Kokkos::Profiling::pushRegion("U stepping");
    calc_uhat_rk2(std::move(Ru), 2);
    Kokkos::Profiling::popRegion();
  } else {
    throw std::runtime_error("Invalid RK stage");
  }
}
} // namespace alps::solver
