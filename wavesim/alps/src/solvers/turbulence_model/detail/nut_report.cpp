//
// Created by xuananqing on 3/29/23.
//

#include "nut_report.h"

#include <common/base/logging.h>
#include <common/device/device_traits.h>
#include <common/kokkos_abstraction/exec_policy.h>
#include <common/kokkos_abstraction/pool_space.h>
#include <common/runtime/async_utils.h>
#include <decomp/block_partition.h>
#include <linear_algebra/sum.h>
#include <solvers/mesh/mesh.h>

#include <mpipp/collectives.h>

namespace alps::solver::detail {

void report_nu_t_stats(MDView<Real const***> const& nu_t,
                       double const                 time,
                       Mesh const&                  mesh,
                       std::string                  logger_name,
                       RotatingFileSinkConfig       file_config)
{
  assert(mesh.extent(2) == nu_t.extent_int(2));

  MDView<Real*, default_memory_pool> const mean_nu_local("mean nu",
                                                         nu_t.extent(2));
  SumXY<Real, std::decay_t<decltype(nu_t)>::memory_space> const sum_functor(
    nu_t, mean_nu_local);

  const auto stream = get_next_stream();
  auto const policy = [&stream, teams = nu_t.extent_int(2)] {
    using Device = Kokkos::DefaultExecutionSpace;
    if constexpr (is_cuda_execution_space_v<Device>) {
      return GridPolicy<>(stream, teams, 16, 32);
    }
    if constexpr (is_hip_execution_space_v<Device>) {
      return GridPolicy<>(stream, teams, 8, 64);
    }
    return GridPolicy<>(stream, teams, Kokkos::AUTO());
  }();
  sum_functor.execute(policy);

  // sum reduction in y direction
  MDView<Real*, default_host_memory_pool> const mean_nu_1("mean nu 1",
                                                          nu_t.extent(2));
  {
    const auto nu_t_h = Kokkos::create_mirror_view(mean_nu_local);
    Kokkos::deep_copy(stream, nu_t_h, mean_nu_local);
    stream.fence();
    mpipp::reduce(nonstd::span(nu_t_h.data(), nu_t_h.size()),
                  nonstd::span(mean_nu_1.data(), mean_nu_1.size()),
                  mpipp::plus<Real>(),
                  0,
                  mesh.comm().axis_comm[1]);
  }

  // gather in z direction, only execute on one column of processors
  if (!mesh.comm().is_first(1)) return;

  MDView<Real* [2], default_host_memory_pool> const nut_all(
    "nut all", mesh.global_extent(2));
  const auto nz_distribution = mesh.partition().get_distribution(2);
  mpipp::gatherv(nonstd::span(&mesh.zz_h(0), mean_nu_1.extent(0)),
                 &nut_all(0, 0),
                 nonstd::span(nz_distribution),
                 0,
                 mesh.comm().axis_comm[2]);
  mpipp::gatherv(nonstd::span(mean_nu_1.data(), mean_nu_1.extent(0)),
                 &nut_all(0, 1),
                 nonstd::span(nz_distribution),
                 0,
                 mesh.comm().axis_comm[2]);

  if (!mesh.comm().is_first(2)) return;

  const auto         nx_global = mesh.global_extent(0);
  const auto         ny_global = mesh.global_extent(1);
  fmt::memory_buffer str_vec;
  for (int k = 1; k < nut_all.extent_int(0) - 1; ++k) {
    fmt::format_to(std::back_inserter(str_vec),
                   "{:.8f}, {:.3e}, {:.4e}\n",
                   time,
                   nut_all(k, 0),
                   nut_all(k, 1) / nx_global / ny_global);
  }

  auto logger = spdlog::get(logger_name);
  if (!logger) {
    logger = create_logger(logger_name, file_config);
    logger->set_pattern("%v");
  }
  logger->info("{}", fmt::to_string(str_vec));
  logger->flush();
}
} // namespace alps::solver::detail
