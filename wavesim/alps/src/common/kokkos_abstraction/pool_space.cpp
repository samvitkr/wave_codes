#include "pool_space.h"

#include <common/memory/dynamic_size_pool.h>
#include <common/runtime/manager.h>

#include <Kokkos_Core.hpp>

#include <string>

namespace alps {
template<class BaseSpace>
typename PoolSpace<BaseSpace>::pool_t& PoolSpace<BaseSpace>::get_allocator()
{
  // return pool_t::global_instance();
  return RuntimeManager::instance().memory_pool<BaseSpace>();
}

template<class BaseSpace>
void* PoolSpace<BaseSpace>::allocate(const execution_space& exec_space,
                                     size_t arg_alloc_size) const
{
  return allocate(exec_space, "[unlabeled]", arg_alloc_size);
}

template<class BaseSpace>
void* PoolSpace<BaseSpace>::allocate(const execution_space& exec_space,
                                     const char*            arg_label,
                                     size_t                 arg_alloc_size,
                                     size_t arg_logical_size) const
{
  return impl_allocate(exec_space, arg_label, arg_alloc_size, arg_logical_size);
}

template<class BaseSpace>
void* PoolSpace<BaseSpace>::allocate(size_t arg_alloc_size) const
{
  return allocate("[unlabeled]", arg_alloc_size);
}

template<class BaseSpace>
void* PoolSpace<BaseSpace>::allocate(const char* arg_label,
                                     size_t      arg_alloc_size,
                                     size_t      arg_logical_size) const
{
  return impl_allocate(arg_label, arg_alloc_size, arg_logical_size);
}

template<class BaseSpace>
void PoolSpace<BaseSpace>::deallocate(void*  arg_alloc_ptr,
                                      size_t arg_alloc_size) const
{
  deallocate("[unlabeled]", arg_alloc_ptr, arg_alloc_size);
}

template<class BaseSpace>
void PoolSpace<BaseSpace>::deallocate(const char* arg_label,
                                      void*       arg_alloc_ptr,
                                      size_t      arg_alloc_size,
                                      size_t      arg_logical_size) const
{
  impl_deallocate(arg_label, arg_alloc_ptr, arg_alloc_size, arg_logical_size);
}

template<class BaseSpace>
void* PoolSpace<BaseSpace>::impl_allocate(
  const execution_space&     exec_space,
  const char*                arg_label,
  size_t                     arg_alloc_size,
  size_t                     arg_logical_size,
  Kokkos::Tools::SpaceHandle arg_handle) const
{
// pooled memory always uses the default device, check if the same device
#if defined(KOKKOS_ENABLE_CUDA)
  if constexpr (std::is_same_v<execution_space, Kokkos::Cuda>) {
    if (exec_space.cuda_device() != execution_space().cuda_device()) {
      throw std::runtime_error(
        "Error: Attempting to allocate pooled memory on a different CUDA "
        "device than the default device.");
    }
  }
#elif defined(KOKKOS_ENABLE_HIP)
  if constexpr (std::is_same_v<execution_space, Kokkos::HIP>) {
    if (exec_space.hip_device() != execution_space().hip_device()) {
      throw std::runtime_error(
        "Error: Attempting to allocate pooled memory on a different HIP "
        "device than the default device.");
    }
  }
#endif
  (void)exec_space; // suppress unused variable warning
  return impl_allocate(arg_label, arg_alloc_size, arg_logical_size, arg_handle);
}

template<class BaseSpace>
void* PoolSpace<BaseSpace>::impl_allocate(
  const char*                arg_label,
  size_t                     arg_alloc_size,
  size_t                     arg_logical_size,
  Kokkos::Tools::SpaceHandle arg_handle) const
{
  const size_t reported_size =
    (arg_logical_size > 0) ? arg_logical_size : arg_alloc_size;

  void* ptr = get_allocator().allocate(arg_alloc_size);
  if (Kokkos::Profiling::profileLibraryLoaded()) {
    Kokkos::Profiling::allocateData(arg_handle, arg_label, ptr, reported_size);
  }
  return ptr;
}

template<class BaseSpace>
void PoolSpace<BaseSpace>::impl_deallocate(
  const char*                arg_label,
  void*                      arg_alloc_ptr,
  size_t                     arg_alloc_size,
  size_t                     arg_logical_size,
  Kokkos::Tools::SpaceHandle arg_handle) const
{
  if (arg_alloc_ptr != nullptr) {
    size_t const reported_size =
      (arg_logical_size > 0) ? arg_logical_size : arg_alloc_size;
    if (Kokkos::Profiling::profileLibraryLoaded()) {
      Kokkos::Profiling::deallocateData(
        arg_handle, arg_label, arg_alloc_ptr, reported_size);
    }
    get_allocator().deallocate(arg_alloc_ptr);
  }
}

template class PoolSpace<Kokkos::HostSpace>;
#if defined(KOKKOS_ENABLE_CUDA)
template class PoolSpace<Kokkos::CudaSpace>;
template class PoolSpace<Kokkos::CudaUVMSpace>;
template class PoolSpace<Kokkos::CudaHostPinnedSpace>;
#elif defined(KOKKOS_ENABLE_HIP)
template class PoolSpace<Kokkos::HIPSpace>;
template class PoolSpace<Kokkos::HIPHostPinnedSpace>;
template class PoolSpace<Kokkos::HIPManagedSpace>;
#endif
} // namespace alps
