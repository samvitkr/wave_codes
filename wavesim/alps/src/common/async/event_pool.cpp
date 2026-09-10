#include "event_pool.h"

#include <common/base/logging.h>

namespace alps::async {

// internal, not thread-safe
template<class ExecutionSpace>
void event_pool<ExecutionSpace>::impl_add(std::size_t n)
{
  if (n == 0) return;
  events_.reserve(total_size_ + n);
  for (size_t i = 0; i < n; ++i) {
    auto new_event = event_t();
    events_.push_back(new_event.release());
    ++total_size_;
  }

  logger_->debug("Created/Added {} events.", n);
}

template<class ExecutionSpace>
void event_pool<ExecutionSpace>::resize(std::size_t new_size)
{
  std::lock_guard<std::mutex> lock(mutex_);
  if (new_size > total_size_) {
    impl_add(new_size - total_size_);
  }
}

template<class ExecutionSpace>
void event_pool<ExecutionSpace>::clear() noexcept
{
  logger_->trace("Clear events.");
  std::lock_guard<std::mutex> lock(mutex_);
  while (!events_.empty()) {
    event_t new_event(events_.back()); // destroyed automatically
    events_.pop_back();
    --total_size_;
  }
}

template<class ExecutionSpace>
event_pool<ExecutionSpace>::~event_pool()
{
  if (total_size_ != events_.size()) {
    logger_->warn("Event pool is being destroyed with {} events still in use.",
                  total_size_ - events_.size());
  }
  clear();
}

template<class ExecutionSpace>
[[nodiscard]] typename event_pool<ExecutionSpace>::event_t
event_pool<ExecutionSpace>::acquire() noexcept
{
  std::lock_guard<std::mutex> lock(mutex_);
  if (events_.empty()) {
    impl_add(default_incremental_size);
  }

  auto e = events_.back();
  events_.pop_back();
  logger_->trace("Release event object ({}) from the pool, {} available now.",
                 fmt::ptr(e),
                 events_.size());
  return event_t(e, this);
}

template<class ExecutionSpace>
void event_pool<ExecutionSpace>::release(handle_t e) noexcept
{
  std::lock_guard<std::mutex> lock(mutex_);
  events_.push_back(e);
  logger_->trace("Event object ({}) returned to the pool, {} available now.",
                 fmt::ptr(e),
                 events_.size());
}

// Specialize for different ExecutionSpace
#define SPECIALIZE_EVENT_POOL_CONSTRUCTOR(EXEC_SPACE)   \
  template<>                                            \
  event_pool<Kokkos::EXEC_SPACE>::event_pool() noexcept \
    : logger_{get_logger(#EXEC_SPACE "_event_pool")}    \
  {                                                     \
    logger_->debug("Event pool created.");              \
  }

SPECIALIZE_EVENT_POOL_CONSTRUCTOR(OpenMP)
template class event_pool<Kokkos::OpenMP>;

#if defined(KOKKOS_ENABLE_CUDA)
SPECIALIZE_EVENT_POOL_CONSTRUCTOR(Cuda)
template class event_pool<Kokkos::Cuda>;
#elif defined(KOKKOS_ENABLE_HIP)
SPECIALIZE_EVENT_POOL_CONSTRUCTOR(HIP)
template class event_pool<Kokkos::HIP>;
#endif

#undef SPECIALIZE_EVENT_POOL_CONSTRUCTOR

} // namespace alps::async
