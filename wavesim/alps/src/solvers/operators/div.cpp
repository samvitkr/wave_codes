#include "div.h"

#include <common/base/macros.h>
#include <common/container/view_types.h>
#include <common/container/view_utils.h>
#include <common/kokkos_abstraction/exec_policy.h>
#include <common/kokkos_abstraction/reducer.h>
#include <common/runtime/async_utils.h>
#include <decomp/ghost_cell_exchange.h>
#include <solvers/mesh/mesh.h>
#include <spectral/spectral.h>

#include <Kokkos_Profiling_ScopedRegion.hpp>

using Kokkos::parallel_for;
using teampolicy_t = alps::GridPolicy<>;
using member_t     = teampolicy_t::member_type;

namespace alps {
namespace solver {

MDView<Real***> div(const Vector3Field<Real***>& vec_u,
                    const Mesh&                  mesh,
                    CenterPt /*unused*/,
                    std::string_view label)
{
  auto const region = Kokkos::Profiling::ScopedRegion(std::string(label));

  const auto& u    = vec_u.x;
  const auto& v    = vec_u.y;
  const auto& w    = vec_u.z;
  auto        ends = local_extents(u);

  /* update_halo_lower_z(mesh.grid, w); */
  auto reqs = async_update_halo_lower_z(mesh.grid.x_pencil(), w, 1);

  // calculate du/dx
  auto u_inner = create_inner_view(u);

  MDView<Real***, default_memory_pool> div_u(
    Kokkos::view_alloc(std::string(label), Kokkos::WithoutInitializing),
    u_inner.layout());
  auto stream = get_next_stream();
  spectral::ddx(div_u, u_inner.view(), mesh.grid, stream);

  // calculate dv/dy
  spectral::ddy_and_add(div_u, create_inner_view(v).view(), mesh.grid, stream);
  auto event = get_device_event();
  enqueue(event, stream);

  // calculate dw/dz
  const auto  hbar = mesh.hbar;
  const auto& dzw  = mesh.dzw;
  auto policy1 = teampolicy_t(stream, ends[1] * (ends[2] - 1), Kokkos::AUTO());
  auto policy2 = teampolicy_t(get_next_stream(), ends[1], Kokkos::AUTO());
  parallel_for(
    "div_uz_CenterPt", policy1, KOKKOS_LAMBDA(member_t team) {
      int  k     = team.league_rank() / ends[1] + 1;
      int  j     = team.league_rank() % ends[1];
      auto alpha = 1 / (dzw(k - 1) * hbar);
      parallel_for(
        Kokkos::TeamVectorRange(team, ends[0]), KOKKOS_TR_LAMBDA(int& i) {
          div_u(i, j, k) += (w(i, j, k) - w(i, j, k - 1)) * alpha;
        });
    });
  reqs.waitall(); // wait for ghost cells
  if (!mesh.grid.comm().is_first(2)) {
    wait_for(event, policy2.space()); // avoid race condition on div_u
    parallel_for(
      "div_uz_CenterPt_bottom", policy2, KOKKOS_LAMBDA(member_t team) {
        int  j     = team.league_rank();
        auto alpha = 1 / (dzw(0 - 1) * hbar);
        parallel_for(
          Kokkos::TeamVectorRange(team, ends[0]), KOKKOS_TR_LAMBDA(int& i) {
            div_u(i, j, 0) += (w(i, j, 0) - w(i, j, -1)) * alpha;
          });
      });
    policy2.space().fence();
  }

  policy1.space().fence();
  return div_u;
}

/// Calculate the maximum divergence of u (halo of u must be updated beforehand)
std::pair<Real, std::array<int, 3>> max_div(const Vector3Field<Real***>& vec_u,
                                            const Mesh&                  mesh,
                                            CenterPt /*unused*/)
{
  Kokkos::Profiling::ScopedRegion region("Max div");

  const auto& u    = vec_u.x;
  const auto& v    = vec_u.y;
  const auto& w    = vec_u.z;
  auto        ends = local_extents(u);

  auto stream = get_next_stream();

  // calculate du/dx
  auto u_inner = create_inner_view(u).view();

  MDView<Real***, default_memory_pool> div_u(
    Kokkos::view_alloc("div_u", Kokkos::WithoutInitializing), u_inner.layout());
  spectral::ddx(div_u, u_inner, mesh.grid, stream);

  // calculate dv/dy
  spectral::ddy_and_add(div_u, create_inner_view(v).view(), mesh.grid, stream);

  // calculate dw/dz
  const auto  hbar   = mesh.hbar;
  const auto& dzw    = mesh.dzw;
  auto        begins = Kokkos::Array<int, 3>{0, 0, 0};
  if (mesh.grid.comm().is_first(2)) begins[2] = 1;
  if (mesh.grid.comm().is_last(2)) ends[2] = ends[2] - 1;
  constexpr auto tile_size = []() -> decltype(begins) {
    using Device = Kokkos::DefaultExecutionSpace;
    if constexpr (alps::is_cuda_execution_space_v<Device>) return {16, 1, 8};
    if constexpr (alps::is_hip_execution_space_v<Device>) return {64, 1, 4};
    return {0, 0, 0};
  }();
  auto policy1 = LoopPolicy<3>(stream, begins, ends, tile_size);

  using max_div_reducer_t = Kokkos::MaxLoc<Real, Kokkos::Array<int, 3>>;
  max_div_reducer_t::value_type maxDivAndLoc{};
  Kokkos::parallel_reduce(
    "max_div reduction",
    policy1,
    KOKKOS_LAMBDA(int i, int j, int k, max_div_reducer_t::value_type& max) {
      const auto alpha = 1 / (dzw(k - 1) * hbar);
      const auto div =
        Kokkos::abs(div_u(i, j, k) + (w(i, j, k) - w(i, j, k - 1)) * alpha);
      if (div > max.val) {
        max.val = div;
        max.loc = {i, j, k};
      }
    },
    max_div_reducer_t(maxDivAndLoc));

  policy1.space().fence();
  const std::array<int, 3> loc{
    maxDivAndLoc.loc[0], maxDivAndLoc.loc[1], maxDivAndLoc.loc[2]};
  return {maxDivAndLoc.val, loc};
}

MDView<Real***> div(const Vector3Field<Real***>& vec_u,
                    const Mesh&                  mesh,
                    NodePt /*unused*/,
                    std::string_view label)
{
  auto const region = Kokkos::Profiling::ScopedRegion(std::string(label));

  const auto& u    = vec_u.x;
  const auto& v    = vec_u.y;
  const auto& w    = vec_u.z;
  auto        ends = local_extents(u);

  /* update_halo_upper_z(mesh.grid, w); */
  auto reqs = async_update_halo_upper_z(mesh.grid.x_pencil(), w, 1);

  // calculate du/dx
  auto u_inner = create_inner_view(u);

  MDView<Real***, default_memory_pool> div_u(
    Kokkos::view_alloc(std::string(label), Kokkos::WithoutInitializing),
    u_inner.layout());
  auto stream = get_next_stream();
  spectral::ddx(div_u, u_inner.view(), mesh.grid, stream);

  // calculate dv/dy
  spectral::ddy_and_add(div_u, create_inner_view(v).view(), mesh.grid, stream);
  auto event = get_device_event();
  enqueue(event, stream);

  // calculate dw/dz
  const auto  hbar = mesh.hbar;
  const auto& dz   = mesh.dz;
  if (mesh.grid.comm().is_last(2)) ends[2] -= 1; // the node is
  auto policy1 = teampolicy_t(stream, ends[1] * (ends[2] - 1), Kokkos::AUTO());
  auto policy2 = teampolicy_t(get_next_stream(), ends[1], Kokkos::AUTO());
  parallel_for(
    "div_uz_NodePt", policy1, KOKKOS_LAMBDA(member_t team) {
      int  k     = team.league_rank() / ends[1];
      int  j     = team.league_rank() % ends[1];
      auto alpha = 1 / (dz(k) * hbar);
      parallel_for(
        Kokkos::TeamVectorRange(team, ends[0]), KOKKOS_TR_LAMBDA(int& i) {
          div_u(i, j, k) += (w(i, j, k + 1) - w(i, j, k)) * alpha;
        });
    });
  reqs.waitall(); // wait for ghost cells
  if (!mesh.grid.comm().is_last(2)) {
    wait_for(event, policy2.space()); // avoid race condition on div_u
    parallel_for(
      "div_uz_NodePt_top", policy2, KOKKOS_LAMBDA(member_t team) {
        int  j     = team.league_rank();
        auto alpha = 1 / (dz(ends[2] - 1) * hbar);
        parallel_for(
          Kokkos::TeamVectorRange(team, ends[0]), KOKKOS_TR_LAMBDA(int& i) {
            div_u(i, j, ends[2] - 1) +=
              (w(i, j, ends[2]) - w(i, j, ends[2] - 1)) * alpha;
          });
      });
    policy2.space().fence();
  }

  policy1.space().fence();
  return div_u;
}

} // namespace solver
} // namespace alps
