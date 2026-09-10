//
// Created by xuanx004 on 12/29/23.
//

#include "translate_velocity.h"

#include <common/container/view_utils.h>
#include <common/kokkos_abstraction/exec_policy.h>
#include <decomp/pencil_plan.h>

#include <mpipp/collectives.h>

void apply_translation_to_fix_bottom_mean_u(
  alps::solver::FlowField const&       flow,
  Kokkos::DefaultExecutionSpace const& space)
{
  const auto& u = flow.u.x;
  using RealT   = typename std::decay_t<decltype(u)>::non_const_value_type;

  RealT u_sum_local{0};
  if (flow.mesh.comm().is_first(2)) {
    // sum bottom u along x and y locally
    Kokkos::parallel_reduce(
      alps::LoopPolicy<2>(
        space, {0, 0}, {alps::local_end(u, 0), alps::local_end(u, 1)}),
      KOKKOS_LAMBDA(int i, int j, RealT& sum) { sum += u(i, j, 0); },
      u_sum_local);
  }

  // allreduce to get the sum of u at the bottom
  RealT u_sum{};
  mpipp::allreduce(u_sum_local, u_sum, mpipp::plus<RealT>{}, flow.mesh.comm());
  auto u_mean = u_sum / flow.mesh.global_extent(0);
  u_mean /= flow.mesh.global_extent(1);

  // subtract the mean from u (including ghost cells)
  Kokkos::parallel_for(
    alps::LoopPolicy<3>(space, alps::begins(u), alps::ends(u)),
    KOKKOS_LAMBDA(int i, int j, int k) { u(i, j, k) -= u_mean; });
}
