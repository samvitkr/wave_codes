#include "pooled_memory_resource.h"

#include "dynamic_size_pool.h"

#include <Kokkos_Core.hpp>

namespace alps::memory {

PooledMemoryResources::PooledMemoryResources()
{
  constexpr std::size_t Alignment =
    512u; // an alignment consistent with recent Cuda behaviour

  using DevMemSpace = Kokkos::DefaultExecutionSpace::memory_space;
  pool_             = std::make_unique<DynamicSizePool<DevMemSpace>>(Alignment);
  if constexpr (has_device_v) {
    host_pool_ =
      std::make_unique<DynamicSizePool<Kokkos::HostSpace>>(Alignment);
    shared_pool_ =
      std::make_unique<DynamicSizePool<Kokkos::SharedSpace>>(Alignment);
    pinned_pool_ =
      std::make_unique<DynamicSizePool<Kokkos::SharedHostPinnedSpace>>(
        Alignment);
  }
}

void PooledMemoryResources::release_all()
{
  if (pool_) pool_->release();
  if (host_pool_) host_pool_->release();
  if (shared_pool_) shared_pool_->release();
  if (pinned_pool_) pinned_pool_->release();
}

PooledMemoryResources::~PooledMemoryResources() = default;

} // namespace alps::memory
