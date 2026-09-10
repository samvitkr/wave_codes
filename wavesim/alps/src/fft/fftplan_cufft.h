#pragma once

#include <common/base/logging_fwd.h>

#include <cufft.h>

namespace alps::fft {

/** Class FFTPlan wraps FFT handle and methods. */
template<class T, class TransformType, class Backend>
class FFTPlan;

class CUFFT;

class R2C;
class C2R;

/** \brief Class for managing CUFFT plan and work area.
 *  This class manages the lifetime of a CUFFT plan and the work area of the
 * plan. When taking over an external CUFFT plan, it is set to user-managed work
 * area. The work area is created on-the-fly using the cufftSetWorkArea function
 * at the time of the plan execution. The plan is executed in a private stream.
 */
template<class T, class TransformType>
class FFTPlan<T, TransformType, CUFFT>
{
 public:
  using handle_t = cufftHandle;

  /// Create an empty and uninitialized plan.
  FFTPlan() = default;

  // Create a plan from a raw CUFFT handle and assign a work space.
  explicit FFTPlan(handle_t cufft_handle);

  FFTPlan(FFTPlan&& src) noexcept;

  FFTPlan& operator=(FFTPlan&& src) noexcept;

  FFTPlan(const FFTPlan&) = delete;

  FFTPlan& operator=(const FFTPlan&) = delete;

  /// Queue the execution of the cufft plan with respect to a given stream.
  void run(void const* input, void* output, cudaStream_t space) const;

  /// Return the raw plan.
  handle_t raw_handle() const noexcept { return plan_; }

  /// Return the work size of the plan.
  std::size_t work_size() const noexcept;

  ~FFTPlan();

  /// Check if there is an associated plan.
  explicit operator bool() const noexcept { return initialized_; }

 private:
  Logger logger_{get_logger("fft")};

  handle_t plan_{};

  bool initialized_{false};

  cudaStream_t stream_{nullptr};

  std::size_t work_size_{};
};

} // namespace alps::fft
