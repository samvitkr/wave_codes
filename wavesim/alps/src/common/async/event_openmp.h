#pragma once

// forward declaration
namespace Kokkos {
class OpenMP;
} // namespace Kokkos

namespace alps::async {

template<class ExecutionSpace>
class event_pool; // forward declaration

template<typename ExecutionSpace>
class queue_event; // forward declaration

template<>
class queue_event<Kokkos::OpenMP>
{
 public:
  using execution_space = Kokkos::OpenMP;
  using handle_t        = void*;
  using pool_t          = event_pool<execution_space>;

  /// Wrap an event in a queue_event
  explicit queue_event(handle_t /*event*/) { /* no-op */ }

  /// Wrap an event from a pool in a queue_event (also reset the status of the
  /// event)
  explicit queue_event(handle_t /*event*/, pool_t* pool)
    : pool_{pool}
  {}

  /// Create a new event
  queue_event();

  /// Create a new event with the given flag
  explicit queue_event(unsigned int flag);

  /// Release the ownership of the event
  handle_t release() noexcept
  {
    pool_ = nullptr;
    return nullptr;
  }

  ~queue_event();

  // move only
  queue_event(queue_event&&) noexcept;

  queue_event& operator=(queue_event&&) noexcept;

  queue_event(queue_event const&) = delete;

  queue_event& operator=(queue_event const&) = delete;

 private:
  pool_t* pool_{};
};
} // namespace alps::async

namespace alps {

void enqueue(async::queue_event<Kokkos::OpenMP>& event,
             Kokkos::OpenMP const& /*space*/);

void wait_for(async::queue_event<Kokkos::OpenMP>& event);

void wait_for(async::queue_event<Kokkos::OpenMP>& event,
              Kokkos::OpenMP const& /*space*/);

inline bool is_complete(async::queue_event<Kokkos::OpenMP>& /*event*/)
{
  return true;
}
} // namespace alps
