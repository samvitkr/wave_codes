#pragma once

#include <Kokkos_Core_fwd.hpp>

#include <memory>
#include <type_traits>

namespace alps::memory {

// Forward declaration
template<class BaseSpace>
class DynamicSizePool;

class PooledMemoryResources
{
 private:
  static constexpr bool has_device_v =
    !std::is_same_v<Kokkos::DefaultExecutionSpace,
                    Kokkos::DefaultHostExecutionSpace>;

 public:
  PooledMemoryResources();
  ~PooledMemoryResources();

  template<typename Space>
  [[nodiscard]] auto& get() const noexcept
  {
    if constexpr (has_device_v) {
      if constexpr (std::is_same_v<Space, Kokkos::HostSpace>) {
        return *host_pool_;
      } else if constexpr (std::is_same_v<Space, Kokkos::SharedSpace>) {
        return *shared_pool_;
      } else if constexpr (std::is_same_v<Space,
                                          Kokkos::SharedHostPinnedSpace>) {
        return *pinned_pool_;
      } else {
        // Default device space (CudaSpace, HIPSpace, etc.)
        return *pool_;
      }
    } else {
      // Host-only
      return *pool_;
    }
  }

  /// Free all releasable memory in all pools
  void release_all();

 private:
  template<typename S>
  using resource_t = std::unique_ptr<DynamicSizePool<S>>;

#if defined(KOKKOS_ENABLE_CUDA)
  resource_t<Kokkos::CudaSpace> pool_;
#elif defined(KOKKOS_ENABLE_HIP)
  resource_t<Kokkos::HIPSpace> pool_;
#else
  resource_t<Kokkos::HostSpace> pool_;
#endif
  resource_t<Kokkos::HostSpace>             host_pool_;
  resource_t<Kokkos::SharedSpace>           shared_pool_;
  resource_t<Kokkos::SharedHostPinnedSpace> pinned_pool_;
};
} // namespace alps::memory
