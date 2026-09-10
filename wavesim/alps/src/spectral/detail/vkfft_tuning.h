#pragma once

#include <common/base/logging.h>
#include <fft/vkfft_structs.h>
#include <spectral/detail/vkfft_app_utils.h>
#include <spectral/spectral_fwd.h>

#include <Kokkos_Macros.hpp>

#include <functional>
#include <string>
#include <vector>

namespace alps::spectral::detail {

struct VkFFTTuningConfig
{
  bool             enabled                      = false;
  std::vector<int> thresholds                   = {8, 12, 16, 24};
  bool             sweep_refine_batch           = true;
  std::vector<int> aim_threads_values           = {32, 64, 128, 256};
  bool             sweep_aim_threads            = true;
  std::vector<int> grouped_batch0_values        = {0, 1};
  bool             sweep_grouped_batch0         = true;
  int              benchmark_samples            = 5;
  double           min_measurable_batch_time_us = 20.0;
  double           target_batch_time_us         = 3000.0;
  double           pilot_batch_time_us          = 250.0;
};

/// Build a VkFFT tuning config from a programmatic control value and the
/// \c ALPS_SKIP_VKFFT_TUNE environment variable.
///
/// Precedence (highest to lowest):
///  1. Programmatic \c Enabled or \c Disabled — always wins.
///  2. \c ALPS_SKIP_VKFFT_TUNE == "1" or a case-insensitive "true" →
///     disabled.
///  3. Default → enabled.
VkFFTTuningConfig make_vkfft_tuning_config(VkFFTTuningControl control);

struct VkFFTSchedulerParams
{
  int register_threshold   = 0;
  int disable_refine_batch = 0;
  int aim_threads          = 0;
  int grouped_batch0       = 0;
};

#if defined(KOKKOS_ENABLE_CUDA)
using VkFFTExecStream = cudaStream_t;
#elif defined(KOKKOS_ENABLE_HIP)
using VkFFTExecStream = hipStream_t;
#else
#error "Unsupported execution space for VkFFT tuning"
#endif

VkApp tune_scheduler_params(
  std::function<VkFFTConfiguration()> const& make_base_config,
  std::function<void(VkFFTApplication*, VkFFTExecStream)> const& execute_plan,
  VkFFTExecStream                                                execute_stream,
  VkFFTTuningConfig const&                                       config,
  Logger                                                         logger);

VkFFTConfiguration with_scheduler_params(VkFFTConfiguration const&   config,
                                         VkFFTSchedulerParams const& params);

struct VkFFTPlanFingerprint
{
  std::string code0_fingerprint;
  bool        dedup_eligible = false;
};

inline bool operator==(VkFFTPlanFingerprint const& lhs,
                       VkFFTPlanFingerprint const& rhs)
{
  return lhs.dedup_eligible == rhs.dedup_eligible
      && lhs.code0_fingerprint == rhs.code0_fingerprint;
}

VkFFTPlanFingerprint fingerprint_vkfft_plan(VkFFTApplication const* app);

} // namespace alps::spectral::detail
