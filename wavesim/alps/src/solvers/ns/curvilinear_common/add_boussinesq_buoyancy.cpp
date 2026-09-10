//
// Created by xuanx004 on 7/24/24.
//

#include "add_boussinesq_buoyancy.h"

#include <common/base/macros.h>
#include <common/container/view_utils.h>
#include <common/device/device_traits.h>
#include <common/kokkos_abstraction/exec_policy.h>
#include <common/kokkos_abstraction/pool_space.h>
#include <common/runtime/async_utils.h>
#include <decomp/mdcomm.h>

#include <mpipp/collectives.h>

namespace alps::solver {

void add_buoyancy_force_with_ref_scalar_curve(
  HaloView<Real***> const&             fz,
  HaloView<Real const***> const&       theta,
  Real const                           theta0,
  Real const                           beta,
  CurvilinearMesh const&               mesh,
  Kokkos::DefaultExecutionSpace const& stream)
{
  auto const& dzw  = mesh.dzw;
  auto const& invJ = mesh.invJ;

  Kokkos::parallel_for(
    "add BoussinesqForce",
    LoopPolicy<3>(stream, local_begins(fz), local_ends(fz)),
    KOKKOS_LAMBDA(int i, int j, int k) {
      auto ratio   = dzw(k - 1) / (dzw(k - 1) + dzw(k));
      auto theta_w = itp2node(theta(i, j, k), theta(i, j, k + 1), ratio);
      fz(i, j, k) -= beta * (theta_w - theta0) * invJ(i, j);
    });
}

void add_buoyancy_force_without_ref_scalar_curve(
  HaloView<Real***> const&             fz,
  HaloView<Real const***> const&       theta,
  Real const                           beta,
  CurvilinearMesh const&               mesh,
  Kokkos::DefaultExecutionSpace const& stream)
{
  // fz = fz - beta * (J^{-1} * theta - \bar{J^{-1} * theta}) where \bar{} is
  // the average for a constant zeta
  auto const& invJ = mesh.invJ;
  auto const  nx   = mesh.extent(0);
  auto const  ny   = mesh.extent(1);
  auto const  nz   = mesh.extent(2);

  // compute the average of J^{-1}*theta
  // length includes the upper ghost cell
  MDView<Real*, default_memory_pool> const local_sum(
    "local_sum " + theta.label(), nz + 1);
  auto const local_sum_h = Kokkos::create_mirror_view(
    PoolSpace<Kokkos::SharedHostPinnedSpace>(), local_sum);

  auto constexpr y_block_size = 16;
  auto const ny_global        = mesh.global_extent(1);
  auto const n_y_blocks       = (ny + y_block_size - 1) / y_block_size;
  auto const policy           = [stream, nz, n_y_blocks] {
    using Device = Kokkos::DefaultExecutionSpace;
    if constexpr (is_cuda_execution_space_v<Device>) {
      return GridPolicy<>(stream, (nz + 1) * n_y_blocks, Kokkos::AUTO(), 16);
    }
    if constexpr (is_hip_execution_space_v<Device>) {
      return GridPolicy<>(stream, (nz + 1) * n_y_blocks, Kokkos::AUTO(), 64);
    }
    return GridPolicy<>(stream, (nz + 1) * n_y_blocks, Kokkos::AUTO());
  }();
  Kokkos::parallel_for(
    "local sum", policy, KOKKOS_LAMBDA(GridPolicy<>::member_type const& team) {
      const auto k       = team.league_rank() / n_y_blocks;
      const auto y_begin = (team.league_rank() % n_y_blocks) * y_block_size;
      const auto y_end   = Kokkos::min(y_begin + y_block_size, ny);

      Real team_sum{0};
      Kokkos::parallel_reduce(
        Kokkos::TeamThreadRange(team, y_begin, y_end),
        KOKKOS_TR_LAMBDA(int j, Real& j_sum) {
          Real i_sum{0};
          Kokkos::parallel_reduce(
            Kokkos::ThreadVectorRange(team, nx),
            KOKKOS_TR_LAMBDA(int i, Real& threadSum) {
              threadSum += (invJ(i, j) * theta(i, j, k)) / nx;
            },
            i_sum);
          Kokkos::single(Kokkos::PerThread(team),
                         [&] { j_sum += i_sum / ny_global; });
        },
        team_sum);

      Kokkos::single(Kokkos::PerTeam(team), [k, team_sum, local_sum] {
        Kokkos::atomic_add(&local_sum(k), team_sum);
      });
    });

  Kokkos::deep_copy(stream, local_sum_h, local_sum);
  stream.fence();

  std::vector<Real> reduce_buffer(local_sum_h.data(),
                                  local_sum_h.data() + local_sum_h.size());
  mpipp::allreduce(nonstd::span(reduce_buffer),
                   nonstd::span(local_sum_h.data(), local_sum_h.size()),
                   mpipp::plus<Real>(),
                   mesh.comm().axis_comm[1]);

  Kokkos::deep_copy(stream, local_sum, local_sum_h);

  auto const& dzw = mesh.dzw;
  Kokkos::parallel_for(
    "add BoussinesqForce",
    LoopPolicy<3>(stream, local_begins(fz), local_ends(fz)),
    KOKKOS_LAMBDA(int i, int j, int k) {
      auto ratio    = dzw(k - 1) / (dzw(k - 1) + dzw(k));
      auto theta_w  = itp2node(theta(i, j, k), theta(i, j, k + 1), ratio);
      auto theta0_w = itp2node(local_sum(k), local_sum(k + 1), ratio);
      fz(i, j, k) -= beta * (theta_w * invJ(i, j) - theta0_w);
    });
}

} // namespace alps::solver
