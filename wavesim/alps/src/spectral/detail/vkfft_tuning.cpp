#include <spectral/detail/vkfft_app_utils.h>
#include <spectral/detail/vkfft_tuning.h>

#include <common/async/event_cuda.h>
#include <common/utils/to_lower_case.h>
#include <fft/vkfft.h>
#include <mpipp/comm.h>
#include <mpipp/environment.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <utility>

template<>
struct fmt::formatter<alps::spectral::detail::VkFFTSchedulerParams>
{
  constexpr auto parse(format_parse_context& ctx) { return ctx.begin(); }

  template<typename FormatContext>
  auto format(alps::spectral::detail::VkFFTSchedulerParams const& params,
              FormatContext&                                      ctx) const
  {
    return fmt::format_to(ctx.out(),
                          "threshold={}, disable_refine_batch={}, "
                          "aim_threads={}, grouped_batch0={}",
                          params.register_threshold,
                          params.disable_refine_batch,
                          params.aim_threads,
                          params.grouped_batch0);
  }
};

namespace alps::spectral::detail {

VkFFTTuningConfig make_vkfft_tuning_config(VkFFTTuningControl control)
{
  auto config = VkFFTTuningConfig{};
  if (control == VkFFTTuningControl::Enabled) {
    config.enabled = true;
    return config;
  }
  if (control == VkFFTTuningControl::Disabled) {
    config.enabled = false;
    return config;
  }
  char const* const skip_env =
    std::getenv("ALPS_SKIP_VKFFT_TUNE"); // NOLINT(concurrency-mt-unsafe)
  if (skip_env != nullptr) {
    auto const skip_value = alps::to_lower_case(skip_env);
    if (skip_value == "1" || skip_value == "true") {
      config.enabled = false;
      return config;
    }
  }
  config.enabled = true;
  return config;
}

VkFFTConfiguration with_scheduler_params(VkFFTConfiguration const&   config,
                                         VkFFTSchedulerParams const& params)
{
  auto tuned_cfg                      = config;
  tuned_cfg.custom_register_threshold = params.register_threshold;
  tuned_cfg.disable_refine_batch      = params.disable_refine_batch;
  if (params.aim_threads > 0) {
    tuned_cfg.aimThreads = static_cast<pfUINT>(params.aim_threads);
  }
  if (params.grouped_batch0 > 0) {
    tuned_cfg.groupedBatch[0] = static_cast<pfUINT>(params.grouped_batch0);
  }
  return tuned_cfg;
}

namespace {
void log_vkfft_tuning_notice_once(Logger const& logger)
{
  static std::once_flag once;
  std::call_once(once, [&logger]() {
    if (mpipp::initialized() && mpipp::COMM_WORLD().rank() != 0) {
      return;
    }
    logger->info(
      "VkFFT scheduler autotuning is enabled: plan creation will benchmark "
      "candidate scheduler settings. Disable with ALPS_SKIP_VKFFT_TUNE=1.");
  });
}

double median(std::vector<double> values)
{
  if (values.empty()) {
    return std::numeric_limits<double>::infinity();
  }

  std::sort(values.begin(), values.end());
  auto const n = values.size();
  if ((n % 2) == 0) {
    return 0.5 * (values[n / 2 - 1] + values[n / 2]);
  }
  return values[n / 2];
}

std::vector<VkFFTSchedulerParams>
build_scheduler_candidates(VkFFTTuningConfig const& config)
{
  if (!config.enabled) return {};

  std::vector<VkFFTSchedulerParams> candidates;
  candidates.push_back({}); // Start with the default config as a candidate
  for (int threshold : config.thresholds) {
    candidates.push_back({threshold, 0, 0, 0});
  }

  // Expand with refine_batch sweep if enabled
  if (config.sweep_refine_batch) {
    auto const base_size = candidates.size();
    for (std::size_t i = 0; i < base_size; ++i) {
      auto params                 = candidates[i];
      params.disable_refine_batch = 1;
      candidates.push_back(params);
    }
  }

  // Expand with aim_threads sweep if enabled
  if (config.sweep_aim_threads) {
    auto const base_size = candidates.size();
    for (std::size_t i = 0; i < base_size; ++i) {
      for (int aim_threads : config.aim_threads_values) {
        auto params        = candidates[i];
        params.aim_threads = aim_threads;
        candidates.push_back(params);
      }
    }
  }

  if (config.sweep_grouped_batch0) {
    auto const base_size = candidates.size();
    for (std::size_t i = 0; i < base_size; ++i) {
      for (int grouped_batch0 : config.grouped_batch0_values) {
        auto params           = candidates[i];
        params.grouped_batch0 = grouped_batch0;
        candidates.push_back(params);
      }
    }
  }

  std::sort(
    candidates.begin(), candidates.end(), [](auto const& lhs, auto const& rhs) {
      if (lhs.register_threshold != rhs.register_threshold) {
        return lhs.register_threshold < rhs.register_threshold;
      }
      if (lhs.disable_refine_batch != rhs.disable_refine_batch) {
        return lhs.disable_refine_batch < rhs.disable_refine_batch;
      }
      if (lhs.aim_threads != rhs.aim_threads) {
        return lhs.aim_threads < rhs.aim_threads;
      }
      return lhs.grouped_batch0 < rhs.grouped_batch0;
    });
  candidates.erase(
    std::unique(candidates.begin(),
                candidates.end(),
                [](auto const& lhs, auto const& rhs) {
                  return lhs.register_threshold == rhs.register_threshold
                      && lhs.disable_refine_batch == rhs.disable_refine_batch
                      && lhs.aim_threads == rhs.aim_threads
                      && lhs.grouped_batch0 == rhs.grouped_batch0;
                }),
    candidates.end());

  return candidates;
}

void append_code0_fingerprint(std::string&     out,
                              VkFFTPlan const* plan,
                              bool&            has_code0)
{
  if (plan == nullptr) return;

  for (int dim = 0; dim < 4; ++dim) {
    auto const n_uploads = static_cast<int>(plan->numAxisUploads[dim]);
    for (int u = 0; u < n_uploads; ++u) {
      out += fmt::format("[dim={}, upload={}]\n", dim, u);

      auto const* code0 = plan->axes[dim][u].specializationConstants.code0;
      if (code0 != nullptr) {
        has_code0 = true;
        out += code0;
      } else {
        out += "<missing-code0>";
      }
      out += "\n<kernel-delimiter>\n";
    }
  }
}

void sync_stream(VkFFTExecStream execute_stream)
{
#if defined(KOKKOS_ENABLE_CUDA)
  cudaStreamSynchronize(execute_stream);
#elif defined(KOKKOS_ENABLE_HIP)
  hipStreamSynchronize(execute_stream);
#endif
}

void record_event(async::queue_event<Kokkos::DefaultExecutionSpace>& event,
                  VkFFTExecStream execute_stream)
{
#if defined(KOKKOS_ENABLE_CUDA)
  cudaError_t status = cudaEventRecord(event.get(), execute_stream);
  if (status != cudaSuccess) {
    throw std::runtime_error("Failed to record CUDA event for VkFFT tuning");
  }
#elif defined(KOKKOS_ENABLE_HIP)
  hipError_t status = hipEventRecord(event.get(), execute_stream);
  if (status != hipSuccess) {
    throw std::runtime_error("Failed to record HIP event for VkFFT tuning");
  }
#endif
}

double wait_and_time_event(
  async::queue_event<Kokkos::DefaultExecutionSpace>& start_event,
  async::queue_event<Kokkos::DefaultExecutionSpace>& end_event)
{
  wait_for(end_event);
  float milliseconds = 0.0f;
#if defined(KOKKOS_ENABLE_CUDA)
  cudaEventElapsedTime(&milliseconds, start_event.get(), end_event.get());
#elif defined(KOKKOS_ENABLE_HIP)
  hipEventElapsedTime(&milliseconds, start_event.get(), end_event.get());
#endif
  return static_cast<double>(milliseconds) * 1000.0; // convert to micro
}

int estimate_launch_count_for_window(double estimated_per_launch_us,
                                     double target_window_us)
{
  if (!std::isfinite(estimated_per_launch_us)
      || (estimated_per_launch_us <= 0.0) || !std::isfinite(target_window_us)
      || (target_window_us <= 0.0)) {
    return 1;
  }

  auto const launch_count =
    static_cast<int>(std::ceil(target_window_us / estimated_per_launch_us));
  return std::max(launch_count, 1);
}

} // namespace

VkFFTPlanFingerprint fingerprint_vkfft_plan(VkFFTApplication const* app)
{
  VkFFTPlanFingerprint result{};
  if (app == nullptr) {
    return result;
  }

  bool has_code0           = false;
  result.code0_fingerprint = "[forward]\n";
  append_code0_fingerprint(
    result.code0_fingerprint, app->localFFTPlan, has_code0);
  result.code0_fingerprint += "[inverse]\n";
  append_code0_fingerprint(
    result.code0_fingerprint, app->localFFTPlan_inverse, has_code0);
  result.dedup_eligible = has_code0;
  return result;
}

VkApp tune_scheduler_params(
  std::function<VkFFTConfiguration()> const& make_base_config,
  std::function<void(VkFFTApplication*, VkFFTExecStream)> const& execute_plan,
  VkFFTExecStream                                                execute_stream,
  VkFFTTuningConfig const&                                       config,
  Logger                                                         logger)
{
  log_vkfft_tuning_notice_once(logger);
  auto const candidates = build_scheduler_candidates(config);
  if (candidates.empty()) {
    logger->debug("VkFFT tuning: no candidates generated; using base config");
    return make_vkfft_app(make_base_config());
  }

#if defined(KOKKOS_ENABLE_CUDA)
  constexpr unsigned int event_flag = cudaEventDefault;
#elif defined(KOKKOS_ENABLE_HIP)
  constexpr unsigned int event_flag = hipEventDefault;
#endif
  auto start_evt =
    async::queue_event<Kokkos::DefaultExecutionSpace>{event_flag};
  auto end_evt = async::queue_event<Kokkos::DefaultExecutionSpace>{event_flag};

  double               best_median = std::numeric_limits<double>::infinity();
  VkFFTSchedulerParams best_params{};
  VkApp                best_app{nullptr, &deleteVkFFT};

  std::vector<VkFFTPlanFingerprint> seen_plans;
  int                               duplicate_count   = 0;
  int                               benchmarked_count = 0;

  // One warmup launch is sufficient; the pilot batch provides additional
  // thermal/clock stabilization before the timed benchmark samples.
  constexpr int warmup_count  = 1;
  auto const    bench_samples = std::max(config.benchmark_samples, 1);

  auto const min_measurable_batch_time_us =
    std::max(config.min_measurable_batch_time_us, 10.0);
  auto const target_batch_time_us =
    std::max(config.target_batch_time_us, min_measurable_batch_time_us);
  auto const pilot_batch_time_us =
    std::max(config.pilot_batch_time_us, min_measurable_batch_time_us);

  for (auto const params : candidates) {
    auto const active_config =
      with_scheduler_params(make_base_config(), params);

    auto app = make_vkfft_app([=]() -> VkFFTConfiguration {
      auto cfg                = active_config;
      cfg.keepShaderCode      = 1; // Keep shader code around for fingerprinting
      cfg.skip_kernel_compile = 1; // Skip compile for fast fingerprinting
      return cfg;
    }());

    auto const fp = fingerprint_vkfft_plan(app.get());
    if (fp.dedup_eligible) {
      auto const is_duplicate =
        std::find(seen_plans.begin(), seen_plans.end(), fp) != seen_plans.end();
      if (is_duplicate) {
        logger->trace("VkFFT tuning: skipping duplicate plan for {}", params);
        ++duplicate_count;
        continue;
      }
      seen_plans.push_back(fp);
    }

    // Re-initialize the app
    app = make_vkfft_app(active_config);
    ++benchmarked_count;

    sync_stream(execute_stream);
    for (int i = 0; i < warmup_count; ++i) {
      execute_plan(app.get(), execute_stream);
    }
    sync_stream(execute_stream);

    record_event(start_evt, execute_stream);
    execute_plan(app.get(), execute_stream);
    record_event(end_evt, execute_stream);
    auto single_launch_time_us = wait_and_time_event(start_evt, end_evt);

    int    pilot_launch_count      = 0;
    double estimated_per_launch_us = single_launch_time_us;
    if (single_launch_time_us <= min_measurable_batch_time_us) {
      pilot_launch_count = estimate_launch_count_for_window(
        single_launch_time_us, pilot_batch_time_us);

      sync_stream(execute_stream);
      record_event(start_evt, execute_stream);
      for (int i = 0; i < pilot_launch_count; ++i) {
        execute_plan(app.get(), execute_stream);
      }
      record_event(end_evt, execute_stream);
      auto const pilot_time_us = wait_and_time_event(start_evt, end_evt);
      // pilot_launch_count is always >= 1 (guaranteed by
      // estimate_launch_count_for_window)
      estimated_per_launch_us =
        pilot_time_us / static_cast<double>(pilot_launch_count);
    }

    auto const bench_count = estimate_launch_count_for_window(
      estimated_per_launch_us, target_batch_time_us);

    std::vector<double> timings;
    timings.reserve(static_cast<std::size_t>(bench_samples));
    for (int k = 0; k < bench_samples; ++k) {
      sync_stream(execute_stream);
      record_event(start_evt, execute_stream);
      for (int i = 0; i < bench_count; ++i) {
        execute_plan(app.get(), execute_stream);
      }
      record_event(end_evt, execute_stream);
      auto const duration = wait_and_time_event(start_evt, end_evt)
                          / static_cast<double>(bench_count);
      timings.push_back(duration);
    }

    auto const candidate_median = median(std::move(timings));
    logger->trace(
      "VkFFT tuning: {}, median_time={:.2f} us", params, candidate_median);

    if (candidate_median < best_median) {
      best_median = candidate_median;
      best_params = params;
      best_app    = std::move(app);
    }
  }

  logger->debug(
    "VkFFT tuning: {} candidates benchmarked, {} duplicates skipped",
    benchmarked_count,
    duplicate_count);

  logger->debug("VkFFT tuning: best {}", best_params);
  return best_app;
}

} // namespace alps::spectral::detail
