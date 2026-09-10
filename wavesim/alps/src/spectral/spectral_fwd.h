#pragma once

#include <fft/fftplan_fwd.h>

namespace alps {
namespace spectral {

template<class T, class Device>
class SpectralPlanBase;

template<class T, class Device, class Backend>
class SpectralPlan;

struct SpectralGrid;

enum class SpectralPostOp
{
  AssignAfterTranspose = 0,
  AddAfterTranspose    = 1
};

/// Programmatic control over VkFFT scheduler tuning.
/// When \c Auto the effective state is determined by the
/// \c ALPS_SKIP_VKFFT_TUNE environment variable (absent or non-matching
/// values → enabled; "1" or case-insensitive "true" → disabled).
/// Explicit \c Enabled and \c Disabled always override the environment.
/// Pass via \c SpectralOptions to \c SpectralGrid constructors.
enum class VkFFTTuningControl
{
  Auto     = 0,
  Enabled  = 1,
  Disabled = 2,
};

} // namespace spectral

using Grid = spectral::SpectralGrid;

} // namespace alps
