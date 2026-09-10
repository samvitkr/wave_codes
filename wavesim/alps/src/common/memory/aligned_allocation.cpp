#include "aligned_allocation.h"

#include <cstdint>

#include <Kokkos_Macros.hpp>

#include <common/base/logging.h>

#if defined(KOKKOS_ENABLE_CUDA)
#include <cuda_runtime_api.h>
#elif defined(KOKKOS_ENABLE_HIP)
#include <hip/hip_runtime.h>
#endif

namespace Kokkos {
class HostSpace;
#if defined(KOKKOS_ENABLE_CUDA)
class CudaSpace;
class CudaUVMSpace;
class CudaHostPinnedSpace;
#elif defined(KOKKOS_ENABLE_HIP)
class HIPSpace;
class HIPManagedSpace;
class HIPHostPinnedSpace;
#endif
} // namespace Kokkos

namespace alps::memory {

template<class BaseSpace>
class base_allocator
{
 public:
  [[nodiscard]] static void* allocate(std::size_t size) noexcept;
  static void                deallocate(void* ptr) noexcept;
};

template<class BaseAllocator>
AlignedAllocation<BaseAllocator>::AlignedAllocation(
  std::size_t alignment) noexcept
  : alignment_{alignment}
{
  static_assert(sizeof(std::size_t) == 8);
  // Round up alignment to next power of 2
  --alignment_;
  alignment_ |= alignment_ >> 1;
  alignment_ |= alignment_ >> 2;
  alignment_ |= alignment_ >> 4;
  alignment_ |= alignment_ >> 8;
  alignment_ |= alignment_ >> 16;
  alignment_ |= alignment_ >> 32;
  ++alignment_;
}

template<class BaseAllocator>
[[nodiscard]] void* AlignedAllocation<BaseAllocator>::allocate(std::size_t size)
{
  static_assert(sizeof(uintptr_t) == sizeof(void*),
                "uintptr_t is not the same size of a pointer.");

  void* ptr   = BaseAllocator::allocate(size + alignment_);
  auto  i_ptr = (uintptr_t)ptr; // NOLINT
  auto* aligned_ptr =
    (void*)((i_ptr - 1u + alignment_) & -alignment_); // NOLINT

  ptr_records.emplace(aligned_ptr, ptr);
  return aligned_ptr;
}

template<class BaseAllocator>
void AlignedAllocation<BaseAllocator>::deallocate(void* ptr)
{
  auto* base_ptr = ptr_records.at(ptr);
  ptr_records.erase(ptr);
  return BaseAllocator::deallocate(base_ptr);
}

template<>
[[nodiscard]] void*
base_allocator<Kokkos::HostSpace>::allocate(std::size_t size) noexcept
{
  return std::malloc(size); // NOLINT
}
template<>
void base_allocator<Kokkos::HostSpace>::deallocate(void* ptr) noexcept
{
  std::free(ptr); // NOLINT
}

#if defined(KOKKOS_ENABLE_CUDA)
template<>
[[nodiscard]] void*
base_allocator<Kokkos::CudaSpace>::allocate(std::size_t size) noexcept
{
  void*             ptr{nullptr};
  cudaError_t const error = ::cudaMalloc(&ptr, size);
  if (error != cudaSuccess) {
    alps::get_logger("memory")->error(
      "cudaMalloc ({}B) failed because: {}", size, cudaGetErrorString(error));
    return nullptr;
  }
  return ptr;
}

template<>
void base_allocator<Kokkos::CudaSpace>::deallocate(void* ptr) noexcept
{
  cudaError_t const error = ::cudaFree(ptr);
  if (error != cudaSuccess) {
    alps::get_logger("memory")->error("cudaFree (ptr: {}) failed because: {}",
                                      fmt::ptr(ptr),
                                      cudaGetErrorString(error));
  }
}

template<>
[[nodiscard]] void*
base_allocator<Kokkos::CudaUVMSpace>::allocate(std::size_t size) noexcept
{
  void*             ptr{nullptr};
  cudaError_t const error = ::cudaMallocManaged(&ptr, size);
  if (error != cudaSuccess) {
    alps::get_logger("memory")->error(
      "cudaMallocManaged ({}B) failed because: {}",
      size,
      cudaGetErrorString(error));
    return nullptr;
  }
  return ptr;
}

template<>
void base_allocator<Kokkos::CudaUVMSpace>::deallocate(void* ptr) noexcept
{
  cudaError_t const error = ::cudaFree(ptr);
  if (error != cudaSuccess) {
    alps::get_logger("memory")->error("cudaFree (ptr: {}) failed because: {}",
                                      fmt::ptr(ptr),
                                      cudaGetErrorString(error));
  }
}

template<>
[[nodiscard]] void*
base_allocator<Kokkos::CudaHostPinnedSpace>::allocate(std::size_t size) noexcept
{
  void*             ptr{nullptr};
  cudaError_t const error = ::cudaMallocHost(&ptr, size);
  if (error != cudaSuccess) {
    alps::get_logger("memory")->error("cudaMallocHost ({}B) failed because: {}",
                                      size,
                                      cudaGetErrorString(error));
    return nullptr;
  }
  return ptr;
}

template<>
void base_allocator<Kokkos::CudaHostPinnedSpace>::deallocate(void* ptr) noexcept
{
  cudaError_t const error = ::cudaFreeHost(ptr);
  if (error != cudaSuccess) {
    alps::get_logger("memory")->error("cudaFree (ptr: {}) failed because: {}",
                                      fmt::ptr(ptr),
                                      cudaGetErrorString(error));
  }
}
#endif

#if defined(KOKKOS_ENABLE_HIP)
template<>
[[nodiscard]] void*
base_allocator<Kokkos::HIPSpace>::allocate(std::size_t size) noexcept
{
  void*            ptr{nullptr};
  hipError_t const error = ::hipMalloc(&ptr, size);
  if (error != hipSuccess) {
    alps::get_logger("memory")->error(
      "hipMalloc ({}B) failed because: {}", size, hipGetErrorString(error));
    return nullptr;
  }
  return ptr;
}

template<>
void base_allocator<Kokkos::HIPSpace>::deallocate(void* ptr) noexcept
{
  hipError_t const error = ::hipFree(ptr);
  if (error != hipSuccess) {
    alps::get_logger("memory")->error("hipFree (ptr: {}) failed because: {}",
                                      fmt::ptr(ptr),
                                      hipGetErrorString(error));
  }
}

template<>
[[nodiscard]] void*
base_allocator<Kokkos::HIPManagedSpace>::allocate(std::size_t size) noexcept
{
  void*            ptr{nullptr};
  hipError_t const error = ::hipMallocManaged(&ptr, size);
  if (error != hipSuccess) {
    alps::get_logger("memory")->error(
      "hipMallocManaged ({}B) failed because: {}",
      size,
      hipGetErrorString(error));
    return nullptr;
  }
  return ptr;
}

template<>
void base_allocator<Kokkos::HIPManagedSpace>::deallocate(void* ptr) noexcept
{
  hipError_t const error = ::hipFree(ptr);
  if (error != hipSuccess) {
    alps::get_logger("memory")->error("hipFree (ptr: {}) failed because: {}",
                                      fmt::ptr(ptr),
                                      hipGetErrorString(error));
  }
}

template<>
[[nodiscard]] void*
base_allocator<Kokkos::HIPHostPinnedSpace>::allocate(std::size_t size) noexcept
{
  void*      ptr{nullptr};
  auto const error = ::hipHostMalloc(&ptr, size, hipHostMallocNonCoherent);
  if (error != hipSuccess) {
    alps::get_logger("memory")->error(
      "hipHostMalloc ({}B) failed because: {}", size, hipGetErrorString(error));
    return nullptr;
  }
  return ptr;
}

template<>
void base_allocator<Kokkos::HIPHostPinnedSpace>::deallocate(void* ptr) noexcept
{
  hipError_t const error = ::hipHostFree(ptr);
  if (error != hipSuccess) {
    alps::get_logger("memory")->error(
      "hipHostFree (ptr: {}) failed because: {}",
      fmt::ptr(ptr),
      hipGetErrorString(error));
  }
}
#endif

template class AlignedAllocation<base_allocator<Kokkos::HostSpace>>;
#if defined(KOKKOS_ENABLE_CUDA)
template class AlignedAllocation<base_allocator<Kokkos::CudaSpace>>;
template class AlignedAllocation<base_allocator<Kokkos::CudaUVMSpace>>;
template class AlignedAllocation<base_allocator<Kokkos::CudaHostPinnedSpace>>;
#elif defined(KOKKOS_ENABLE_HIP)
template class AlignedAllocation<base_allocator<Kokkos::HIPSpace>>;
template class AlignedAllocation<base_allocator<Kokkos::HIPManagedSpace>>;
template class AlignedAllocation<base_allocator<Kokkos::HIPHostPinnedSpace>>;
#endif

} // namespace alps::memory
