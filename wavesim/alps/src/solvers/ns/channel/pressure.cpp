#include "pressure.h"

#include <common/container/view_utils.h>
#include <common/kokkos_abstraction/exec_policy.h>
#include <common/kokkos_abstraction/pool_space.h>
#include <common/runtime/async_utils.h>
#include <solvers/field/flow_field.h>
#include <solvers/mesh/mesh.h>
#include <solvers/ns/pressure_eqn_coeffs.h>
#include <solvers/poisson/tridiagonal_solver.h>
#include <spectral/spectral.h>

namespace alps::solver {

PressureEqn::PressureEqn(const FlowField& flow)
  : flow_field{flow}
  , d(MDView<Real***, default_memory_pool>(
      Kokkos::view_alloc("d", Kokkos::WithoutInitializing),
      flow.mesh.grid
        .get_r2c_xy_output_layout<Real, Kokkos::DefaultExecutionSpace>()))
  , dl(MDView<Real***, default_memory_pool>(
      Kokkos::view_alloc("dl", Kokkos::WithoutInitializing),
      d.layout()))
  , du(MDView<Real***, default_memory_pool>(
      Kokkos::view_alloc("du", Kokkos::WithoutInitializing),
      d.layout()))
  , solver{create_tridiagonal_solver<Real>(d.extent_int(0),
                                           d.extent_int(1),
                                           d.extent_int(2),
                                           flow.mesh.grid.comm().axis_comm[2])}
{}

void PressureEqn::initialize() const
{
  set_pressure_eqn_coefficients(d,
                                dl,
                                du,
                                PressureBCType::NEUMANN,
                                flow_field.mesh,
                                Kokkos::DefaultExecutionSpace());
  Kokkos::DefaultExecutionSpace().fence();
  solver->setup(d, dl, du);
}

PressureEqn::~PressureEqn() = default;

void solve(const MDView<Real***>& div_u, const PressureEqn& peqn, Real dt)
{
  using Kokkos::parallel_for;
  using Kokkos::TeamThreadRange;

  const auto&                  flow      = peqn.flow_field;
  const Mesh&                  mesh      = flow.mesh;
  const HaloView<const Real*>& dz        = mesh.dz;
  const HaloView<const Real*>& dzw       = mesh.dzw;
  const auto                   hbar      = mesh.hbar;
  const auto&                  grid      = mesh.grid;
  const auto                   is_top    = grid.comm().is_last(2);
  const auto                   is_bottom = grid.comm().is_first(2);

  const auto& d  = peqn.d;
  const auto& dl = peqn.dl;
  const auto& du = peqn.du;

  auto ends = local_ends(flow.pp);

  auto rhs = create_inner_view(flow.pp);

  MDView<Real***, default_memory_pool> sigma_y(
    Kokkos::view_alloc("", Kokkos::WithoutInitializing), d.layout());

  auto stream0 = get_next_stream();

  // if the solver overwrites coefficients, the coefficients need to be reset
  if (peqn.solver->is_coeff_overwritten) {
    set_pressure_eqn_coefficients(
      peqn.d, peqn.dl, peqn.du, PressureBCType::NEUMANN, mesh, stream0);
  }

  auto           stream1   = get_next_stream();
  auto           stream2   = get_next_stream();
  constexpr auto tile_size = [&]() -> Kokkos::Array<std::int64_t, 3> {
    using Device = decltype(stream0);
    if constexpr (is_cuda_execution_space_v<Device>) return {16, 1, 8};
    if constexpr (is_hip_execution_space_v<Device>) return {64, 1, 4};
    return {0, 0, 0};
  }();
  auto policy =
    LoopPolicy<3>(stream0,
                  {0, 0, is_bottom ? 1 : 0},
                  {ends[0], ends[1], is_top ? ends[2] - 1 : ends[2]},
                  tile_size);
  parallel_for(
    policy, KOKKOS_LAMBDA(int i, int j, int k) {
      auto delta   = dzw(k - 1) * hbar;
      rhs(i, j, k) = (div_u(i, j, k) / dt) * delta * delta;
    });
  if (is_bottom) {
    parallel_for(
      LoopPolicy<2>(stream1, {0, 0}, {ends[0], ends[1]}),
      KOKKOS_LAMBDA(int i, int j) {
        auto delta   = dzw(0) * hbar;
        auto alpha   = 2 + dz(1) / dz(0);
        rhs(i, j, 0) = ((div_u(i, j, 1) / dt) * delta * delta) / alpha;
      });
  }
  if (is_top) {
    parallel_for(
      LoopPolicy<2>(stream2, {0, 0}, {ends[0], ends[1]}),
      KOKKOS_LAMBDA(int i, int j) {
        auto delta = dzw(ends[2] - 3) * hbar;
        auto alpha = -(2 + dz(ends[2] - 3) / dz(ends[2] - 2));
        rhs(i, j, ends[2] - 1) =
          ((div_u(i, j, ends[2] - 2) / dt) * delta * delta) / alpha;
      });
  }

  if (is_bottom) stream1.fence();
  if (is_top) stream2.fence();

  spectral::fft_r2c_xy(sigma_y, rhs.view(), grid, stream0);

  if (grid.comm().is_first(1) && is_top) {
    auto sigma_y_0 = Kokkos::subview(sigma_y, 0, 0, ends[2] - 1);
    Kokkos::deep_copy(stream0, sigma_y_0, 0);
  }

  stream0.fence();
  peqn.solver->comm_.barrier();

  peqn.solver->solve(sigma_y, d, dl, du);

  spectral::fft_c2r_xy(rhs.view(), sigma_y, true, grid, stream0);

  stream0.fence();
}

} // namespace alps::solver
