//
// Created by xuanx004 on 7/22/24.
//

#include "boussinesq.h"

#include <common/container/view_utils.h>
#include <common/device/device_traits.h>
#include <common/kokkos_abstraction/exec_policy.h>
#include <common/kokkos_abstraction/pool_space.h>
#include <common/real_type.h>
#include <common/runtime/async_utils.h>
#include <decomp/mdcomm.h>
#include <linear_algebra/sum.h>
#include <solvers/mesh/mesh.h>

#include <fmt/format.h>
#include <mpipp/collectives.h>

namespace alps::solver {

BoussinesqOptions BoussinesqOptions::parse_from(ConfigTable const& config)
{
  if (!config.contains(config_key)) {
    return {};
  }
  auto const table = config.extract_table(config_key);

  BoussinesqOptions options;

  options.scalar_id    = table.get_value_or<int>("scalar_id", -1);
  options.scalar_label = table.get_value_or<std::string>("scalar_label", "");

  options.Ri = table.get_value<double>("Ri");

  // use maximum value to indicate that the reference scalar value is not set
  try {
    options.ref_scalar = table.get_value<double>("theta0");
  } catch (std::exception const&) {
    /* do nothing */
  }

  if (options.scalar_id < 0 && options.scalar_label.empty()) {
    throw std::runtime_error("Scalar field id or label must be provided");
  }
  if (options.scalar_id >= 0 && !options.scalar_label.empty()) {
    throw std::runtime_error(
      "Only one of scalar_id and scalar_label should be provided");
  }

  if (auto const abs_Ri = (Real)std::abs(options.Ri);
      abs_Ri < std::numeric_limits<Real>::epsilon() && abs_Ri > 0) {
    throw std::runtime_error(
      fmt::format("Richardson parameter is too small: {}", options.Ri));
  }

  if (options.Ri != 0) {
    options.enabled = true;
  }

  return options;
}

// functions are defined here because they are linked in multiple solvers
void add_buoyancy_force_with_ref_scalar(
  HaloView<Real***> const&             fz,
  HaloView<Real const***> const&       theta,
  Real const                           theta0,
  Real const                           beta,
  Mesh const&                          mesh,
  Kokkos::DefaultExecutionSpace const& stream)
{
  auto const& dzw = mesh.dzw;

  Kokkos::parallel_for(
    "add BoussinesqForce",
    LoopPolicy<3>(stream, local_begins(fz), local_ends(fz)),
    KOKKOS_LAMBDA(int i, int j, int k) {
      auto ratio   = dzw(k - 1) / (dzw(k - 1) + dzw(k));
      auto theta_w = itp2node(theta(i, j, k), theta(i, j, k + 1), ratio);
      fz(i, j, k) -= beta * (theta_w - theta0);
    });
}

void add_buoyancy_force_without_ref_scalar(
  HaloView<Real***> const&             fz,
  HaloView<Real const***> const&       theta,
  Real const                           beta,
  Mesh const&                          mesh,
  Kokkos::DefaultExecutionSpace const& stream)
{
  auto const [nx, ny, nz] = local_extents(theta);

  // compute the average of theta
  // length includes the upper ghost cell
  MDView<Real*, default_memory_pool> local_sum("local_sum " + theta.label(),
                                               nz + 1);

  auto const& local_sum_h = Kokkos::create_mirror_view(
    PoolSpace<Kokkos::SharedHostPinnedSpace>(), local_sum);

  SumXY<Real, Kokkos::DefaultExecutionSpace::memory_space> local_sum_op(
    subview(theta, Kokkos::ALL, Kokkos::ALL, index_range(0, nz + 1)).view(),
    local_sum);
  auto const policy = [stream, teams = nz + 1] {
    using Device = Kokkos::DefaultExecutionSpace;
    if constexpr (is_cuda_execution_space_v<Device>) {
      return GridPolicy<>(stream, teams, 16, 32);
    }
    if constexpr (is_hip_execution_space_v<Device>) {
      return GridPolicy<>(stream, teams, 8, 64);
    }
    return GridPolicy<>(stream, teams, Kokkos::AUTO());
  }();
  local_sum_op.execute(policy);

  Kokkos::deep_copy(stream, local_sum_h, local_sum);

  stream.fence();
  std::vector<Real> local_mean(local_sum_h.size());
  for (std::size_t i = 0; i < local_mean.size(); ++i) {
    local_mean[i] = local_sum_h(i) / mesh.global_extent(0);
    local_mean[i] /= mesh.global_extent(1);
  }
  mpipp::allreduce(nonstd::span(local_mean),
                   nonstd::span(local_sum_h.data(), local_mean.size()),
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
      fz(i, j, k) -= beta * (theta_w - theta0_w);
    });
}

} // namespace alps::solver
