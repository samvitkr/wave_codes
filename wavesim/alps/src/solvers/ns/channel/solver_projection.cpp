#include "pressure.h"
#include "solver_ab2.h"
#include "solver_ab2cn.h"
#include <common/base/macros.h>
#include <common/container/view_types.h>
#include <common/container/view_utils.h>
#include <common/device/devices.h>
#include <common/runtime/async_utils.h>
#include <decomp/ghost_cell_exchange.h>
#include <solvers/field/flow_field.h>
#include <solvers/operators/div.h>
#include <spectral/spectral.h>

#include <Kokkos_Core.hpp>
#include <Kokkos_Profiling_ScopedRegion.hpp>

namespace alps::solver {

namespace {
template<typename SolverType>
void projection_correct(SolverType const& solver)
{
  using Kokkos::parallel_for;

  Kokkos::Profiling::ScopedRegion region("Correct U");

  const auto& u       = solver.flow_field.u.x;
  const auto& v       = solver.flow_field.u.y;
  const auto& w       = solver.flow_field.u.z;
  const auto& pp      = solver.flow_field.pp;
  const Mesh& mesh    = solver.flow_field.mesh;
  const auto& dz      = mesh.dz;
  const auto& grid    = mesh.grid;
  const auto  hbar    = mesh.hbar;
  const auto  delta_t = (Real)solver.dt;

  const auto ends = local_extents(u);

  const auto stream1 = get_next_stream();
  const auto stream2 = get_next_stream();

  // update the ghost cells of pressure
  auto reqs = async_update_halo_z(mesh.grid.x_pencil(), pp, 1);

  const auto p_inner = create_inner_view(pp);

  const MDView<Real***, default_memory_pool> p_x(
    Kokkos::view_alloc("dpx", Kokkos::WithoutInitializing), p_inner.layout());
  const MDView<Real***, default_memory_pool> p_y(
    Kokkos::view_alloc("dpy", Kokkos::WithoutInitializing), p_inner.layout());

  spectral::ddx(p_x, p_inner.view(), grid, stream1);
  // u = u - delta_t * dp/dx
  using teampolicy_t = Kokkos::TeamPolicy<Kokkos::IndexType<int>>;
  using member_t     = teampolicy_t::member_type;
  parallel_for(
    "correction_u",
    teampolicy_t(stream1, ends[1] * ends[2], Kokkos::AUTO),
    KOKKOS_LAMBDA(member_t team) {
      int k = team.league_rank() / ends[1];
      int j = team.league_rank() % ends[1];
      parallel_for(
        Kokkos::TeamVectorRange(team, ends[0]), KOKKOS_TR_LAMBDA(int& i) {
          u(i, j, k) = Kokkos::fma(p_x(i, j, k), -delta_t, u(i, j, k));
        });
    });
  reqs.waitall(); // overlap ghost cell exchange of p with ddx

  // v = v - delta_t * dp/dy
  // w = w - delta_t * dp/dz
  spectral::ddy(p_y, p_inner.view(), grid, stream2);
  parallel_for(
    "correction_vw",
    teampolicy_t(stream2, ends[1] * ends[2], Kokkos::AUTO),
    KOKKOS_LAMBDA(member_t team) {
      const int  k     = team.league_rank() / ends[1];
      const int  j     = team.league_rank() % ends[1];
      const auto alpha = -delta_t / (dz(k) * hbar);
      parallel_for(
        Kokkos::TeamVectorRange(team, ends[0]), KOKKOS_TR_LAMBDA(const int& i) {
          v(i, j, k) = Kokkos::fma(p_y(i, j, k), -delta_t, v(i, j, k));
          w(i, j, k) =
            Kokkos::fma(pp(i, j, k + 1) - pp(i, j, k), alpha, w(i, j, k));
        });
    });

  stream1.fence();
  stream2.fence();

  ALPS_CHECK_LAST_DEVICE_ERROR();

  // set boundary condition
  solver.apply_bc(stream1);

  // dealias and update ghost cells
  // operations are interleaved to overlap communication and dealiasing
  auto event = get_device_event();

  const auto u_inner = create_inner_view(u).view();
  spectral::dealias(u_inner, grid, stream1);
  enqueue(event, stream1);

  const auto v_inner = create_inner_view(v).view();
  spectral::dealias(v_inner, grid, stream1);
  wait_for(event); // wait for the dealiasing of u before updating ghost cells
  enqueue(event, stream1);
  update_halo_z(mesh.partition(), u, 2);

  const auto w_inner = create_inner_view(w).view();
  spectral::dealias(w_inner, grid, stream1);
  wait_for(event); // wait for the dealiasing of v before updating ghost cells
  enqueue(event, stream1);
  update_halo_z(mesh.partition(), v, 3);
  wait_for(event);
  update_halo_z(mesh.partition(), w, 4);

  ALPS_CHECK_LAST_DEVICE_ERROR();
}
} // anonymous namespace

void ChannelFlowSolverAB2::project() const
{
  Kokkos::Profiling::pushRegion("Project U: div");
  auto div_u = div(flow_field.u, flow_field.mesh, CenterPt());
  Kokkos::Profiling::popRegion();

  Kokkos::Profiling::pushRegion("Project U: poisson");
  solve(div_u, *peqn, (Real)dt);

  ALPS_CHECK_LAST_DEVICE_ERROR();
  Kokkos::Profiling::popRegion();
}

void ChannelFlowSolverAB2::correct() const
{
  projection_correct(*this);
}

void ChannelFlowSolverAB2CN::project() const
{
  Kokkos::Profiling::pushRegion("Project U: div");
  auto div_u = div(flow_field.u, flow_field.mesh, CenterPt());
  Kokkos::Profiling::popRegion();

  Kokkos::Profiling::pushRegion("Project U: poisson");
  solve(div_u, *peqn, (Real)dt);

  ALPS_CHECK_LAST_DEVICE_ERROR();
  Kokkos::Profiling::popRegion();
}

void ChannelFlowSolverAB2CN::correct() const
{
  projection_correct(*this);
}

} // namespace alps::solver
