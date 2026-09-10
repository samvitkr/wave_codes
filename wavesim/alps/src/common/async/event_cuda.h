#pragma once

#include <Kokkos_Macros.hpp>

#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
#if defined(KOKKOS_ENABLE_CUDA)
#include <cuda_runtime_api.h>
#elif defined(KOKKOS_ENABLE_HIP)
#include <hip/hip_runtime.h>
#endif

#include <utility>

namespace Kokkos {
#if defined(KOKKOS_ENABLE_CUDA)
class Cuda;
using DefaultExecutionSpace = Cuda;
#elif defined(KOKKOS_ENABLE_HIP)
class HIP;
using DefaultExecutionSpace = HIP;
#endif
class OpenMP;
} // namespace Kokkos

namespace alps::async {

template<class ExecutionSpace>
class event_pool; // forward declaration

template<typename ExecutionSpace>
class queue_event; // forward declaration

template<>
class queue_event<Kokkos::DefaultExecutionSpace>
{
 public:
  using execution_space = Kokkos::DefaultExecutionSpace;
#if defined(KOKKOS_ENABLE_CUDA)
  using handle_t = CUevent;
#elif defined(KOKKOS_ENABLE_HIP)
  using handle_t = hipEvent_t;
#endif
  using pool_t = event_pool<execution_space>;

  /// Wrap an event in a queue_event
  explicit queue_event(handle_t event)
    : event_{event}
  {}

  /// Wrap an event from a pool in a queue_event
  explicit queue_event(handle_t event, pool_t* pool)
    : event_{event}
    , pool_{pool}
  {}

  /// Create a new event
  queue_event();

  /// Create a new event with the given flag
  explicit queue_event(unsigned int flag);

  /// Get the underlying event
  handle_t get() const noexcept { return event_; }

  /// Release the ownership of the event
  handle_t release() noexcept
  {
    pool_ = nullptr;
    return std::exchange(event_, nullptr);
  }

  ~queue_event();

  // move only
  queue_event(queue_event&&) noexcept;

  queue_event& operator=(queue_event&&) noexcept;

  queue_event(queue_event const&) = delete;

  queue_event& operator=(queue_event const&) = delete;

 private:
  handle_t event_{};
  pool_t*  pool_{};

#if defined(KOKKOS_ENABLE_CUDA)
  static constexpr unsigned int DEFAULT_FLAG =
    cudaEventDisableTiming | cudaEventBlockingSync;
#elif defined(KOKKOS_ENABLE_HIP)
  static constexpr unsigned int DEFAULT_FLAG =
    hipEventDisableTiming | hipEventBlockingSync;
#endif
};
} // namespace alps::async

namespace alps {
void enqueue(async::queue_event<Kokkos::DefaultExecutionSpace>& event,
             Kokkos::DefaultExecutionSpace const&               space);

void wait_for(async::queue_event<Kokkos::DefaultExecutionSpace>& event);

void wait_for(async::queue_event<Kokkos::DefaultExecutionSpace>& event,
              Kokkos::DefaultExecutionSpace const&               space);

void wait_for(async::queue_event<Kokkos::DefaultExecutionSpace>& event,
              Kokkos::OpenMP const& /*space*/);

bool is_complete(async::queue_event<Kokkos::DefaultExecutionSpace>& event);
} // namespace alps
#endif
