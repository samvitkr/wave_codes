#pragma once

namespace alps {

/// @brief Get the number of devices seen by this process
int get_device_count();

/// @brief Set the active device
/// @note This function is intended to be used to set the device for threads
/// other than the main thread; normally, it is not necessary to call this
/// function directly
void set_device(int device_id);

/// @brief Start profiling, works with, e.g. range capture
void start_profiling();

/// @brief Stop profiling, works with, e.g. range capture
void stop_profiling();

namespace detail {
/// check last error
void check_last_device_error(char const* file, int line);
} // namespace detail

// NOLINTNEXTLINE(cppcoreguidelines-macro-usage)
#define ALPS_CHECK_LAST_DEVICE_ERROR() \
  ::alps::detail::check_last_device_error(__FILE__, __LINE__)

} // namespace alps
