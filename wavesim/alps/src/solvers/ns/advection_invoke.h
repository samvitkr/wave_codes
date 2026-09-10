#pragma once

#include <common/container/vector_field.h>
#include <common/container/view_utils.h>
#include <common/device/device_traits.h>
#include <common/kokkos_abstraction/exec_policy.h>
#include <common/real_type.h>
#include <common/runtime/async_utils.h>
#include <decomp/mdcomm.h>

#include <Kokkos_Core.hpp>

namespace alps::solver {

template<class Functor, class FieldType>
void invoke_advection_fluxes_functor(const Vector3Field<Real***>&   fluxes,
                                     const HaloView<Real const***>& f,
                                     const FieldType&               flow)
{
  Kokkos::Profiling::pushRegion("advection " + f.label());
  std::vector<Kokkos::DefaultExecutionSpace> streams;

  auto       begins  = local_begins(flow.u.x);
  auto       ends    = local_ends(flow.u.x);
  const auto functor = Functor(fluxes, flow.mesh, flow.u, f);
  if (flow.mesh.comm().is_first(2)) {
    begins[2] = 1;

    LoopPolicy<2, typename Functor::Bottom_0> policy_bottom_0(
      streams.emplace_back(get_next_stream()),
      {begins[0], begins[1]},
      {ends[0], ends[1]},
      {16, 8});
    Kokkos::parallel_for(
      "advection " + f.label() + " bottom_0", policy_bottom_0, functor);
  }
  if (flow.mesh.comm().is_last(2)) {
    ends[2] -= 2;

    LoopPolicy<2, typename Functor::Top_NzMinus1> policy_top_1(
      streams.emplace_back(get_next_stream()),
      {begins[0], begins[1]},
      {ends[0], ends[1]},
      {16, 8});
    LoopPolicy<2, typename Functor::Top_NzMinus2> policy_top_2(
      policy_top_1.space(),
      {begins[0], begins[1]},
      {ends[0], ends[1]},
      {16, 8});
    Kokkos::parallel_for(
      "advection " + f.label() + " top_2", policy_top_2, functor);
    Kokkos::parallel_for(
      "advection " + f.label() + " top_1", policy_top_1, functor);
  }

  constexpr auto tile = []() -> Kokkos::Array<std::int64_t, 3> {
    using Device = Kokkos::DefaultExecutionSpace;
    if constexpr (is_cuda_execution_space_v<Device>) return {16, 1, 8};
    if constexpr (is_hip_execution_space_v<Device>) return {64, 1, 4};
    return {0, 0, 0};
  }();
  LoopPolicy<3, typename Functor::General> policy(
    streams.emplace_back(get_next_stream()), begins, ends, tile);
  Kokkos::parallel_for("advection " + f.label(), policy, functor);

  fence(std::move(streams));
  Kokkos::Profiling::popRegion();
}

} // namespace alps::solver
