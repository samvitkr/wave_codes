#include "streams.h"

#include <common/base/logging.h>

#include <Kokkos_Core.hpp>

#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
#include <algorithm>

#if defined(KOKKOS_ENABLE_CUDA)
#include <cuda.h>
#elif defined(KOKKOS_ENABLE_HIP)
#include <hip/hip_runtime.h>
#endif

namespace alps {

namespace detail {
void clear_cuda_instances(std::vector<Kokkos::DefaultExecutionSpace>& streams)
{
#if defined(KOKKOS_ENABLE_CUDA)
  while (!streams.empty()) {
    streams.erase(std::remove_if(streams.begin(),
                                 streams.end(),
                                 [](const Kokkos::Cuda& space) {
                                   return cudaSuccess
                                       == cudaStreamQuery(space.cuda_stream());
                                 }),
                  streams.end());
  }
#elif defined(KOKKOS_ENABLE_HIP)
  while (!streams.empty()) {
    streams.erase(std::remove_if(streams.begin(),
                                 streams.end(),
                                 [](const Kokkos::HIP& space) {
                                   auto err =
                                     hipStreamQuery(space.hip_stream());
                                   return err == hipSuccess;
                                 }),
                  streams.end());
  }
#endif
}
} // namespace detail

StreamPool<Kokkos::DefaultExecutionSpace>::StreamPool(int pool_size) noexcept
  : size{pool_size}
  , logger_(get_logger("cuda_stream_pool"))
{
  logger_->debug("Create stream pool with {} streams.", size);
}

void StreamPool<Kokkos::DefaultExecutionSpace>::initialize() const
{
  logger_->debug("Initialize stream pool with {} streams.", size);
  streams_.reserve(size);
  for (int s = 0; s < size; ++s) {
#if defined(KOKKOS_ENABLE_CUDA)
    cudaStream_t new_stream{};
    KOKKOS_IMPL_CUDA_SAFE_CALL(cudaStreamCreate(&new_stream));
    // assign managed Cuda instances
    streams_.emplace_back(new_stream, true);
#elif defined(KOKKOS_ENABLE_HIP)
    hipStream_t new_stream{};
    KOKKOS_IMPL_HIP_SAFE_CALL(hipStreamCreate(&new_stream));
    // assign managed HIP instances
    streams_.emplace_back(new_stream, true);
#endif
  }
}

void StreamPool<Kokkos::DefaultExecutionSpace>::clear() const
{
  logger_->trace("Clean up stream pool with {} streams.", streams_.size());
  detail::clear_cuda_instances(streams_);
}

StreamPool<Kokkos::DefaultExecutionSpace>::~StreamPool()
{
  clear();
}

} // namespace alps
#endif
