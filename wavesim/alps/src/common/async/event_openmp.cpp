#include "event_openmp.h"

#include "event_pool.h"

namespace alps::async {
queue_event<Kokkos::OpenMP>::queue_event() = default;

queue_event<Kokkos::OpenMP>::queue_event(unsigned int /*flag*/)
{ /* no-op */
}

queue_event<Kokkos::OpenMP>::queue_event(queue_event&& other) noexcept
  : pool_{other.pool_}
{
  other.pool_ = nullptr;
}

queue_event<Kokkos::OpenMP>&
queue_event<Kokkos::OpenMP>::operator=(queue_event&& rhs) noexcept
{
  if (this->pool_ != nullptr) {
    this->pool_->release(nullptr);
  }
  pool_     = rhs.pool_;
  rhs.pool_ = nullptr;
  return *this;
}

queue_event<Kokkos::OpenMP>::~queue_event()
{
  if (pool_ != nullptr) {
    pool_->release(nullptr);
  }
}
} // namespace alps::async

namespace alps {
void enqueue(async::queue_event<Kokkos::OpenMP>& /*event*/,
             Kokkos::OpenMP const& /*space*/)
{}

void wait_for(async::queue_event<Kokkos::OpenMP>& /*event*/) {}

void wait_for(async::queue_event<Kokkos::OpenMP>& /*event*/,
              Kokkos::OpenMP const& /*space*/)
{}
} // namespace alps
