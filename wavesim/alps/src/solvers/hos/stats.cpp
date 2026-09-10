//
// Created by xuanx004 on 5/25/24.
//

#include "stats.h"

#include <common/base/logging.h>
#include <common/container/view_utils.h>
#include <common/device/device_traits.h>
#include <common/kokkos_abstraction/pool_space.h>
#include <common/math.h>
#include <common/runtime/async_utils.h>
#include <linear_algebra/mean_variance.h>
#include <solvers/hos/field.h>
#include <spectral/spectral.h>

#include <Kokkos_Core.hpp>
#include <Kokkos_StdAlgorithms.hpp>
#include <mpipp/collectives.h>
#include <spdlog/sinks/rotating_file_sink.h>

namespace alps::solver::hos {

bool check_validity(MDView<Real const**> const& eta, Grid const& grid)
{
  if (!eta.span_is_contiguous()) {
    throw std::runtime_error("eta is not contiguous");
  }
  MDView<Real const*> eta_flat(eta.data(), eta.size());

  auto local_validity = Kokkos::Experimental::none_of(
    "check inf or nan",
    Kokkos::DefaultExecutionSpace(),
    eta_flat,
    KOKKOS_LAMBDA(Real val) {
      return Kokkos::isinf(val) || Kokkos::isnan(val);
    });

  // Use int for portability
  int global_validity = 0;
  int validity_int    = local_validity ? 1 : 0;
  mpipp::allreduce(
    validity_int, global_validity, mpipp::min<int>(), grid.comm());
  return global_validity == 1;
}

struct WaveStatsFileLoggerConfig
{
  static constexpr auto max_size    = (std::size_t)1024 * 1024 * 10; // 10MB
  static constexpr auto max_files   = 2;
  static constexpr auto filename    = "wave_stats.log";
  static constexpr auto logger_name = "wave_stats_logger";
};

namespace {
[[nodiscard]] Logger get_hos_wave_stats_file_logger()
{
  auto logger = spdlog::get(WaveStatsFileLoggerConfig::logger_name);
  if (!logger) {
    auto new_logger =
      spdlog::rotating_logger_mt(WaveStatsFileLoggerConfig::logger_name,
                                 WaveStatsFileLoggerConfig::filename,
                                 WaveStatsFileLoggerConfig::max_size,
                                 WaveStatsFileLoggerConfig::max_files);
    new_logger->set_pattern("[%T.%e] %v");
    new_logger->flush_on(spdlog::level::info);
    return new_logger;
  }
  return logger;
}
} // namespace

void report_wave_stats(HOSField const& field, double time)
{
  auto eta_rms = compute_eta_rms(field.eta(), field.grid);

  auto comm = field.grid.comm();
  if (comm.rank() == 0) {
    get_hos_wave_stats_file_logger()->info("t={}\neta_rms={}", time, eta_rms);
  }
  mpipp::barrier(comm);
}

Real compute_eta_rms(MDView<Real const**> const& eta, Grid const& grid)
{
  auto tile_size = []() -> Kokkos::Array<long, 2> {
    using Device = Kokkos::DefaultExecutionSpace;
    if constexpr (is_cuda_execution_space_v<Device>) return {16, 16};
    if constexpr (is_hip_execution_space_v<Device>) return {};
    return {};
  }();
  auto policy = LoopPolicy<2>(get_next_stream(),
                              {0, 0},
                              {eta.extent_int(0), eta.extent_int(1)},
                              tile_size);
  Real eta2{};
  Kokkos::parallel_reduce(
    "sum eta^2",
    policy,
    KOKKOS_LAMBDA(int i, int j, Real& sum) { sum += eta(i, j) * eta(i, j); },
    eta2);
  policy.space().fence();

  Real eta2_global{};
  mpipp::allreduce(eta2, eta2_global, mpipp::plus<Real>(), grid.comm());

  return (Real)std::sqrt((double)eta2_global / (double)grid.global_extent(0)
                         / (double)grid.global_extent(1));
}

MDView<Real**, Kokkos::HostSpace>
compute_eta_density_2d(Kokkos::View<Real const**> const& eta, Grid const& grid)
{
  if (grid.extent(2) == 1) {
    throw std::runtime_error("compute_eta_spectrum_2d: 2D grid is required");
  }

  MDView<Real const** [1]> eta_1(eta.data(), eta.extent(0), eta.extent(1));
  MDView<Real** [1], default_memory_pool> eta_kxy(
    "eta kxy",
    grid.get_r2c_xy_output_layout<Real, Kokkos::DefaultExecutionSpace>());

  auto stream = get_next_stream();
  spectral::fft_r2c_xy(eta_kxy, eta_1, grid, stream);

  auto kx0      = grid.pex;
  auto ky0      = grid.pey;
  auto max_kx   = grid.global_extent(0) / 2;
  auto max_ky   = grid.global_extent(1) / 2;
  auto offset_0 = grid.offset(0, Pencil::Y);
  auto offset_1 = grid.offset(1, Pencil::Y);

  MDView<Real**, default_memory_pool> Ek_local(
    "Ek", grid.global_extent(0) / 2, grid.global_extent(1));
  auto Ek_local_host =
    create_mirror_view(PoolSpace<Kokkos::SharedHostPinnedSpace>(), Ek_local);
  Kokkos::parallel_for(
    "compute local Ek",
    LoopPolicy<2>(
      stream, {0, 0}, {eta_kxy.extent_int(0), eta_kxy.extent_int(1)}),
    KOKKOS_LAMBDA(int i, int j) {
      auto i_ky = (offset_0 + i) / 2;
      auto i_kx = (offset_1 + j) / 2;
      if (i_ky < max_ky && i_kx < max_kx) {
        auto E = square(eta_kxy(i, j, 0) / Real(max_kx * 2) / Real(max_ky * 2));
        E *= (i_kx != 0) ? 2 : 1;
        E *= (i_ky != 0) ? 2 : 1;
        Kokkos::atomic_add(&Ek_local(i_kx, i_ky), E / (Real)kx0 / (Real)ky0);
      }
    });

  Kokkos::deep_copy(stream, Ek_local_host, Ek_local);
  stream.fence();
  eta_kxy = {}; // release memory

  auto Ek_global = create_mirror(PoolSpace<Kokkos::HostSpace>(), Ek_local_host);
  mpipp::allreduce(nonstd::span(Ek_local_host.data(), Ek_local_host.span()),
                   nonstd::span(Ek_global.data(), Ek_global.span()),
                   mpipp::plus<Real>(),
                   grid.comm());

  return Ek_global;
}

} // namespace alps::solver::hos
