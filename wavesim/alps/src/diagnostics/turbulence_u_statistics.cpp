//
// Created by xuananqing on 4/1/23.
//

#include "field_statistics.h"
#include <common/container/view_types.h>
#include <common/container/view_utils.h>
#include <common/kokkos_abstraction/pool_space.h>
#include <common/real_type.h>
#include <common/runtime/async_utils.h>
#include <decomp/block_partition.h>
#include <linear_algebra/covariance.h>
#include <linear_algebra/mean_variance.h>
#include <solvers/field/flow_field.h>

#include <Kokkos_Core.hpp>
#include <mpipp/collectives.h>
#include <spdlog/fmt/ostr.h>

#include <algorithm>
#include <string>

namespace alps::diagnostics {

class TurbulenceVelocityStatistics : public FieldStatistics
{
 public:
  using base_t = FieldStatistics;

  TurbulenceVelocityStatistics(
    solver::FlowField const&                flow_field_,
    std::unordered_set<FieldStatisticsType> field_statistics_types_);

  void calculate() override;

  std::vector<std::string>
  get_variable_labels(FieldStatisticsType type) const override;

  void write_to_tecplot(std::ostream&       out,
                        FieldStatisticsType type) const override;

  void cleanup() override;

  void calculate_y();

  void calculate_xy();

 private:
  template<typename DT>
  using view_t = Kokkos::View<DT, Kokkos::LayoutLeft, Kokkos::HostSpace>;

  void write_to_tecplot_y(std::ostream& out) const;

  void write_to_tecplot_xy(std::ostream& out) const;

  solver::FlowField const& flow_field;

  view_t<Real**> u_mean_y;
  view_t<Real**> v_mean_y;
  view_t<Real**> w_mean_y;
  view_t<Real**> u_variance_y;
  view_t<Real**> v_variance_y;
  view_t<Real**> w_variance_y;
  view_t<Real**> uw_cov_y;

  view_t<Real*> u_mean_xy;
  view_t<Real*> v_mean_xy;
  view_t<Real*> w_mean_xy;
  view_t<Real*> u_variance_xy;
  view_t<Real*> v_variance_xy;
  view_t<Real*> w_variance_xy;
  view_t<Real*> uw_cov_xy;

  /// Statistics already computed to avoid re-computing them
  std::unordered_set<FieldStatisticsType> calculated_types;
};

TurbulenceVelocityStatistics::TurbulenceVelocityStatistics(
  solver::FlowField const&                flow_field_,
  std::unordered_set<FieldStatisticsType> field_statistics_types_)
  : base_t{std::move(field_statistics_types_)}
  , flow_field{flow_field_}
{}

void TurbulenceVelocityStatistics::calculate()
{
  if (statistic_types.find(FieldStatisticsType::Y) != statistic_types.end()) {
    calculate_y();
  }

  if (statistic_types.find(FieldStatisticsType::XY) != statistic_types.end()) {
    calculate_xy();
  }
}

void TurbulenceVelocityStatistics::cleanup()
{
  u_mean_y      = {};
  v_mean_y      = {};
  w_mean_y      = {};
  u_variance_y  = {};
  v_variance_y  = {};
  w_variance_y  = {};
  uw_cov_y      = {};
  u_mean_xy     = {};
  v_mean_xy     = {};
  w_mean_xy     = {};
  u_variance_xy = {};
  v_variance_xy = {};
  w_variance_xy = {};
  uw_cov_xy     = {};
  calculated_types.clear();
}

std::vector<std::string> TurbulenceVelocityStatistics::get_variable_labels(
  FieldStatisticsType /*type*/) const
{
  return {"U", "V", "W", "uu", "vv", "ww", "uw"};
}

void TurbulenceVelocityStatistics::write_to_tecplot(
  std::ostream&       out,
  FieldStatisticsType type) const
{
  if (type == FieldStatisticsType::Y) {
    write_to_tecplot_y(out);
  } else if (type == FieldStatisticsType::XY) {
    write_to_tecplot_xy(out);
  }
}

void TurbulenceVelocityStatistics::calculate_y()
{
  if (calculated_types.find(FieldStatisticsType::Y) != calculated_types.end()) {
    return;
  }

  auto const  stream1   = get_next_stream();
  auto const  stream2   = get_next_stream();
  auto const& comm      = flow_field.mesh.comm();
  auto const  is_bottom = comm.is_first(2);
  auto const  is_top    = comm.is_last(2);

  // Interpolate w to cell center
  auto const&                                w = flow_field.u.z;
  MDView<Real***, default_memory_pool> const w_center(
    "w_center", local_extent(w, 0), local_extent(w, 1), local_extent(w, 2));

  auto const nz      = local_extent(w, 2);
  auto const z_begin = is_bottom ? 1 : 0;
  auto const z_end   = is_top ? nz - 1 : nz;
  Kokkos::parallel_for(
    "interpolate w",
    LoopPolicy<3>(
      stream1, {0, 0, z_begin}, {local_end(w, 0), local_end(w, 1), z_end}),
    KOKKOS_LAMBDA(int i, int j, int k) {
      w_center(i, j, k) = solver::itp2center(w(i, j, k), w(i, j, k - 1));
    });
  if (is_top) {
    Kokkos::deep_copy(stream2,
                      subview(w_center, ALL, ALL, nz - 1),
                      subview(w, ALL, ALL, nz - 2).view());
  }
  if (is_bottom) {
    Kokkos::deep_copy(
      stream2, subview(w_center, ALL, ALL, 0), subview(w, ALL, ALL, 0).view());
  }
  stream1.fence();
  stream2.fence();

  // The mean and variance are calculated and then stored (copied) in the class
  // Because the statistics are collected on the first processor in the y
  // direction, allocation and copy only on this process
  MeanVariance<Real, default_memory_pool> mean_variance(w_center);
  mean_variance.calculate_local();
  mean_variance.gather(0, comm.axis_comm[1]);
  if (comm.is_first(1)) {
    w_mean_y = Kokkos::create_mirror(Kokkos::HostSpace(), mean_variance.mean);
    w_variance_y = Kokkos::create_mirror(Kokkos::HostSpace(), mean_variance.M2);
    std::swap(w_mean_y, mean_variance.mean);
    std::swap(w_variance_y, mean_variance.M2);
  }

  auto const u    = create_inner_view(flow_field.u.x).view();
  mean_variance.f = u;
  mean_variance.calculate_local();
  mean_variance.gather(0, comm.axis_comm[1]);
  if (comm.is_first(1)) {
    u_mean_y = Kokkos::create_mirror(Kokkos::HostSpace(), mean_variance.mean);
    u_variance_y = Kokkos::create_mirror(Kokkos::HostSpace(), mean_variance.M2);
    std::swap(u_mean_y, mean_variance.mean);
    std::swap(u_variance_y, mean_variance.M2);
  }

  auto const v    = create_inner_view(flow_field.u.y).view();
  mean_variance.f = v;
  mean_variance.calculate_local();
  mean_variance.gather(0, comm.axis_comm[1]);
  if (comm.is_first(1)) {
    v_mean_y = Kokkos::create_mirror(Kokkos::HostSpace(), mean_variance.mean);
    v_variance_y = Kokkos::create_mirror(Kokkos::HostSpace(), mean_variance.M2);
    std::swap(v_mean_y, mean_variance.mean);
    std::swap(v_variance_y, mean_variance.M2);
  }

  Covariance<Real, default_memory_pool> cov_variance(u, w_center);
  cov_variance.calculate_local();
  cov_variance.gather(0, comm.axis_comm[1]);
  if (comm.is_first(1)) {
    uw_cov_y = Kokkos::create_mirror(Kokkos::HostSpace(), cov_variance.M2);
    std::swap(uw_cov_y, cov_variance.M2);
  }

  mpipp::barrier(comm); // Barrier may be optional but is here anyway to avoid
                        // potential synchronization issues

  calculated_types.insert(FieldStatisticsType::Y);
}

void TurbulenceVelocityStatistics::calculate_xy()
{
  // Calculate the mean and variance in the y direction if not already done
  if (calculated_types.find(FieldStatisticsType::Y) == calculated_types.end()) {
    calculate_y();
  }

  auto const& comm = flow_field.mesh.comm();
  if (!comm.is_first(1)) {
    // Only the first processor in the y direction computes the mean and
    // variance in the x direction. Barrier is optional but is set anyway to
    // avoid potential synchronization issues
    mpipp::barrier(comm);
    return;
  }

  auto const nx = local_extent(u_mean_y, 0);
  auto const nz = local_extent(u_mean_y, 1);
  u_mean_xy     = decltype(u_mean_xy)("u_mean_xy", nz);
  v_mean_xy     = decltype(v_mean_xy)("v_mean_xy", nz);
  w_mean_xy     = decltype(w_mean_xy)("w_mean_xy", nz);
  u_variance_xy = decltype(u_variance_xy)("u_var_xy", nz);
  v_variance_xy = decltype(v_variance_xy)("v_var_xy", nz);
  w_variance_xy = decltype(w_variance_xy)("w_var_xy", nz);
  uw_cov_xy     = decltype(uw_cov_xy)("uw_cov_xy", nz);

  // Further do reduction in the x direction
  // Numerically Stable Parallel Computation of (Co-)Variance, Schubert & Gertz
#pragma omp parallel for
  for (int k = 0; k < nz; ++k) {
    for (int i = 0; i < nx; ++i) {
      auto const delta_u = u_mean_y(i, k) - u_mean_xy(k);
      u_mean_xy(k) += delta_u / (i + 1);
      auto const delta_u_new = u_mean_y(i, k) - u_mean_xy(k);
      u_variance_xy(k) += u_variance_y(i, k) + delta_u * delta_u_new;

      auto const delta_v = v_mean_y(i, k) - v_mean_xy(k);
      v_mean_xy(k) += delta_v / (i + 1);
      auto const delta_v_new = v_mean_y(i, k) - v_mean_xy(k);
      v_variance_xy(k) += v_variance_y(i, k) + delta_v * delta_v_new;

      auto const delta_w = w_mean_y(i, k) - w_mean_xy(k);
      w_mean_xy(k) += delta_w / (i + 1);
      auto const delta_w_new = w_mean_y(i, k) - w_mean_xy(k);
      w_variance_xy(k) += w_variance_y(i, k) + delta_w * delta_w_new;

      uw_cov_xy(k) += uw_cov_y(i, k) + delta_u_new * delta_w;
    }
    u_variance_xy(k) /= nx;
    v_variance_xy(k) /= nx;
    w_variance_xy(k) /= nx;
    uw_cov_xy(k) /= nx;
  }

  mpipp::barrier(comm);

  calculated_types.insert(FieldStatisticsType::XY);
}

void TurbulenceVelocityStatistics::write_to_tecplot_xy(std::ostream& out) const
{
  auto const& mesh = flow_field.mesh;
  auto const& comm = flow_field.mesh.comm();
  if (!comm.is_first(1)) {
    // Only the first processor in the y direction has the statistics. Barrier
    // is set to avoid potential synchronization issues
    mpipp::barrier(comm);
    return;
  }

  // collect in the z direction
  auto const nz_distribution = mesh.partition().get_distribution(2);
  auto const z_offsets       = mesh.partition().get_offsets(2);

  auto const                             nz = mesh.global_extent(2);
  MDView<Real*, Kokkos::HostSpace> const gathered("gathered", nz);

  for (auto const& local : {u_mean_xy,
                            v_mean_xy,
                            w_mean_xy,
                            u_variance_xy,
                            v_variance_xy,
                            w_variance_xy,
                            uw_cov_xy}) {
    mpipp::gatherv(nonstd::span(local.data(), local.extent(0)),
                   gathered.data(),
                   nonstd::span(nz_distribution),
                   nonstd::span(z_offsets),
                   0,
                   comm.axis_comm[2]);

    // Write out the gathered variable on the root processor
    if (!comm.is_first(2)) continue;

    for (int k = 0; k < nz; ++k) {
      fmt::print(out, "{} ", gathered(k));
    }
    fmt::print(out, "\n");
  }
  mpipp::barrier(comm);
}

void TurbulenceVelocityStatistics::write_to_tecplot_y(std::ostream& out) const
{
  auto const& mesh = flow_field.mesh;
  auto const& comm = flow_field.mesh.comm();
  if (!comm.is_first(1)) {
    // Only the first processor in the y direction has the statistics. Barrier
    // is set to avoid potential synchronization issues
    mpipp::barrier(comm);
    return;
  }

  // collect in the z direction
  auto const nz_distribution = mesh.partition().get_distribution(2);
  auto const z_offsets       = mesh.partition().get_offsets(2);

  auto const                              nx = mesh.global_extent(0);
  auto const                              nz = mesh.global_extent(2);
  MDView<Real**, Kokkos::HostSpace> const gathered("gathered", nx, nz);

  // The number of elements and offsets on each processor need to be
  // multiplied by nx
  std::vector<int> distribution;
  std::transform(nz_distribution.cbegin(),
                 nz_distribution.cend(),
                 std::back_inserter(distribution),
                 [&](int n) { return n * nx; });
  std::vector<int> offsets;
  std::transform(z_offsets.cbegin(),
                 z_offsets.cend(),
                 std::back_inserter(offsets),
                 [&](int n) { return n * nx; });
  for (auto const& local : {u_mean_y,
                            v_mean_y,
                            w_mean_y,
                            u_variance_y,
                            v_variance_y,
                            w_variance_y,
                            uw_cov_y}) {
    mpipp::gatherv(nonstd::span(local.data(), local.span()),
                   gathered.data(),
                   nonstd::span(distribution),
                   nonstd::span(offsets),
                   0,
                   comm.axis_comm[2]);

    if (!comm.is_first(2)) continue;

    for (int k = 0; k < nz; ++k) {
      for (int i = 0; i < nx; ++i) {
        fmt::print(out, "{} ", gathered(i, k));
      }
      fmt::print(out, "\n");
    }
  }
  mpipp::barrier(comm);
}

std::unique_ptr<FieldStatistics> create_velocity_statistics(
  solver::FlowField const&                flow_field,
  std::unordered_set<FieldStatisticsType> field_statistics_types)
{
  return std::make_unique<TurbulenceVelocityStatistics>(flow_field,
                                                        field_statistics_types);
}
} // namespace alps::diagnostics
