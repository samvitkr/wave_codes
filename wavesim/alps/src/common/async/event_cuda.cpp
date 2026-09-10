#include "event_cuda.h"

#include "event_pool.h"

#include <Kokkos_Core.hpp>

namespace alps::async {
queue_event<Kokkos::DefaultExecutionSpace>::queue_event()
  : queue_event(DEFAULT_FLAG)
{}

queue_event<Kokkos::DefaultExecutionSpace>::queue_event(unsigned int flag)
{
#if defined(KOKKOS_ENABLE_CUDA)
  cudaEventCreateWithFlags(&event_, flag);
#elif defined(KOKKOS_ENABLE_HIP)
  KOKKOS_IMPL_HIP_SAFE_CALL(hipEventCreateWithFlags(&event_, flag));
#endif
}

queue_event<Kokkos::DefaultExecutionSpace>::queue_event(
  queue_event&& other) noexcept
  : event_{other.event_}
  , pool_{other.pool_}
{
  other.event_ = nullptr;
  other.pool_  = nullptr;
}

queue_event<Kokkos::DefaultExecutionSpace>&
queue_event<Kokkos::DefaultExecutionSpace>::operator=(
  queue_event&& rhs) noexcept
{
  if (this->event_ != nullptr) {
    if (this->pool_ != nullptr) {
      this->pool_->release(this->event_);
    } else {
#if defined(KOKKOS_ENABLE_CUDA)
      cudaEventDestroy(event_);
#elif defined(KOKKOS_ENABLE_HIP)
      KOKKOS_IMPL_HIP_SAFE_CALL(hipEventDestroy(event_));
#endif
    }
  }
  event_     = rhs.event_;
  pool_      = rhs.pool_;
  rhs.event_ = nullptr;
  rhs.pool_  = nullptr;
  return *this;
}

queue_event<Kokkos::DefaultExecutionSpace>::~queue_event()
{
  if (this->event_ == nullptr) {
    return;
  }
  if (pool_ != nullptr) {
    pool_->release(event_);
  } else {
#if defined(KOKKOS_ENABLE_CUDA)
    cudaEventDestroy(event_);
#elif defined(KOKKOS_ENABLE_HIP)
    KOKKOS_IMPL_HIP_SAFE_CALL(hipEventDestroy(event_));
#endif
  }
}
} // namespace alps::async

namespace alps {
void enqueue(async::queue_event<Kokkos::DefaultExecutionSpace>& event,
             Kokkos::DefaultExecutionSpace const&               space)
{
#if defined(KOKKOS_ENABLE_CUDA)
  cudaEventRecord(event.get(), space.cuda_stream());
#elif defined(KOKKOS_ENABLE_HIP)
  KOKKOS_IMPL_HIP_SAFE_CALL(hipEventRecord(event.get(), space.hip_stream()));
#endif
}

void wait_for(async::queue_event<Kokkos::DefaultExecutionSpace>& event)
{
#if defined(KOKKOS_ENABLE_CUDA)
  cudaEventSynchronize(event.get());
#elif defined(KOKKOS_ENABLE_HIP)
  KOKKOS_IMPL_HIP_SAFE_CALL(hipEventSynchronize(event.get()));
#endif
}

void wait_for(async::queue_event<Kokkos::DefaultExecutionSpace>& event,
              Kokkos::DefaultExecutionSpace const&               space)
{
#if defined(KOKKOS_ENABLE_CUDA)
  cudaStreamWaitEvent(space.cuda_stream(), event.get());
#elif defined(KOKKOS_ENABLE_HIP)
  KOKKOS_IMPL_HIP_SAFE_CALL(
    hipStreamWaitEvent(space.hip_stream(), event.get(), 0));
#endif
}

void wait_for(async::queue_event<Kokkos::DefaultExecutionSpace>& event,
              Kokkos::OpenMP const& /*space*/)
{
#if defined(KOKKOS_ENABLE_CUDA)
  cudaEventSynchronize(event.get());
#elif defined(KOKKOS_ENABLE_HIP)
  KOKKOS_IMPL_HIP_SAFE_CALL(hipEventSynchronize(event.get()));
#endif
}

bool is_complete(async::queue_event<Kokkos::DefaultExecutionSpace>& event)
{
#if defined(KOKKOS_ENABLE_CUDA)
  cudaError_t status = cudaEventQuery(event.get());
  if (status == cudaSuccess) return true;
  if (status == cudaErrorNotReady) return false;
  throw std::runtime_error("cudaEventQuery failed");
#elif defined(KOKKOS_ENABLE_HIP)
  hipError_t status = hipEventQuery(event.get());
  if (status == hipSuccess) return true;
  if (status == hipErrorNotReady) return false;
  throw std::runtime_error("hipEventQuery failed");
#else
  return true;
#endif
}
} // namespace alps
