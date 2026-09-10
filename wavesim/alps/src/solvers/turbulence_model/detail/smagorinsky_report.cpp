#include "smagorinsky_report.h"

#include <common/base/logging.h>
#include <common/kokkos_abstraction/pool_space.h>
#include <decomp/block_partition.h>
#include <solvers/mesh/mesh.h>
#include <solvers/turbulence_model/dynamic_smagorinsky.h>
#include <solvers/turbulence_model/smagorinsky.h>

#include <mpipp/collectives.h>

namespace alps::solver::detail {

void report_smagorinsky_C0(ConstantSmagorinsky const&              sgs_model,
                           MDView<Real*, Kokkos::HostSpace> const& scaled_delta,
                           double                                  time,
                           Mesh const&                             mesh,
                           std::string const&                      logger_name,
                           RotatingFileSinkConfig                  file_config)
{
  auto constexpr PI = Kokkos::numbers::pi_v<double>;
  // The filter width considers the effective resolution from dealias
  auto const  dxf  = Real(PI / (double)mesh.pex / int(mesh.global_extent(0) / 2)
                        * sgs_model.options_.horizontal_filter_width_scale);
  auto const  dyf  = Real(PI / (double)mesh.pey / int(mesh.global_extent(1) / 2)
                        * sgs_model.options_.horizontal_filter_width_scale);
  auto const  C0   = sgs_model.options_.C0;
  auto const  nz   = mesh.extent(2);
  auto const& comm = mesh.comm();

  if (!comm.is_first(1)) return;

  std::vector<Real> zw_all(mesh.global_extent(2));
  std::vector<Real> c0_all(mesh.global_extent(2));
  const auto        nz_distribution = mesh.partition().get_distribution(2);
  mpipp::gatherv(nonstd::span(scaled_delta.data(), nz),
                 c0_all.data(),
                 nonstd::span(nz_distribution),
                 0,
                 comm.axis_comm[2]);
  mpipp::gatherv(nonstd::span(&mesh.zw_h(0), nz),
                 zw_all.data(),
                 nonstd::span(nz_distribution),
                 0,
                 comm.axis_comm[2]);
  if (!comm.is_first(2)) return;

  // Write all data to a string buffer then log it to prevent one instant split
  // over different files
  fmt::memory_buffer str_vec;
  for (std::size_t k = 1; k < c0_all.size() - 1; ++k) {
    const auto dzf   = (zw_all[k] - zw_all[k - 1]) * mesh.hbar;
    const auto delta = Kokkos::cbrt(dxf * dyf * dzf);
    const auto zz    = itp2center(zw_all[k], zw_all[k - 1]);
    fmt::format_to(std::back_inserter(str_vec),
                   "{:.8f}, {:.3e}, {:.4e}\n",
                   time,
                   zz,
                   (Real)C0 * c0_all[k] / delta);
  }
  auto logger = spdlog::get(logger_name);
  if (!logger) {
    logger = create_logger(logger_name, file_config);
    logger->set_pattern("%v");
  }
  logger->info("{}", fmt::to_string(str_vec));
  logger->flush();
}

namespace {
void report_dynamic_smagorinsky_C0_impl(
  MDView<Real*> const&   C0_Delta2,
  double const           horizontal_filter_width_scale,
  double const           time,
  Mesh const&            mesh,
  std::string const&     logger_name,
  RotatingFileSinkConfig file_config)
{
  auto const nz   = mesh.extent(2);
  auto const comm = mesh.comm();
  if (!comm.is_first(1)) return; // only execute on one column of processors

  auto constexpr PI = Kokkos::numbers::pi_v<double>;
  const auto dxf = Real(PI / (double)mesh.pex / int(mesh.global_extent(0) / 2)
                        * horizontal_filter_width_scale);
  const auto dyf = Real(PI / (double)mesh.pey / int(mesh.global_extent(1) / 2)
                        * horizontal_filter_width_scale);
  const auto C0_Delta2_h =
    Kokkos::create_mirror_view(Kokkos::HostSpace(), C0_Delta2);
  Kokkos::deep_copy(C0_Delta2_h, C0_Delta2);

  std::vector<Real> zw_all(mesh.global_extent(2));
  std::vector<Real> c0_all(mesh.global_extent(2));
  const auto        nz_distribution = mesh.partition().get_distribution(2);
  mpipp::gatherv(nonstd::span(&C0_Delta2_h(0), nz),
                 c0_all.data(),
                 nonstd::span(nz_distribution),
                 0,
                 comm.axis_comm[2]);
  mpipp::gatherv(nonstd::span(&mesh.zw_h(0), nz),
                 zw_all.data(),
                 nonstd::span(nz_distribution),
                 0,
                 comm.axis_comm[2]);

  if (!comm.is_first(2)) return;

  // Write all data to a string buffer then log it to prevent one instant split
  // over different files
  fmt::memory_buffer str_vec;
  for (std::size_t k = 1; k < c0_all.size() - 1; ++k) {
    const auto dzf   = (zw_all[k] - zw_all[k - 1]) * mesh.hbar;
    const auto delta = Kokkos::cbrt(dxf * dyf * dzf);
    const auto zz    = itp2center(zw_all[k], zw_all[k - 1]);
    fmt::format_to(std::back_inserter(str_vec),
                   "{:.8f}, {:.3e}, {:.4e}\n",
                   time,
                   zz,
                   Kokkos::sqrt(c0_all[k]) / delta);
  }
  auto logger = spdlog::get(logger_name);
  if (!logger) {
    logger = create_logger(logger_name, file_config);
    logger->set_pattern("%v");
  }
  logger->info("{}", fmt::to_string(str_vec));
  logger->flush();
}
} // namespace

void report_dynamic_smagorinsky_C0(DynamicSmagorinsky const& sgs_model,
                                   double                    time,
                                   Mesh const&               mesh,
                                   std::string const&        logger_name,
                                   RotatingFileSinkConfig    file_config)
{
  report_dynamic_smagorinsky_C0_impl(
    sgs_model.C0_Delta2,
    sgs_model.options_.horizontal_filter_width_scale,
    time,
    mesh,
    logger_name,
    file_config);
  mesh.comm().comm.barrier();
}

void report_dynamic_smagorinsky_C0(DynamicSmagorinskyScalar const& sgs_model,
                                   double                          time,
                                   Mesh const&                     mesh,
                                   std::string const&              logger_name,
                                   RotatingFileSinkConfig          file_config)
{
  report_dynamic_smagorinsky_C0_impl(
    sgs_model.C0_Delta2,
    sgs_model.sgs_model_->options_.horizontal_filter_width_scale,
    time,
    mesh,
    logger_name,
    file_config);
  mesh.comm().comm.barrier();
}
} // namespace alps::solver::detail
