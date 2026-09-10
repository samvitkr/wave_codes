#include "pool_space.h"

#include <Kokkos_Core.hpp>

#include <string>

// Instantiate SharedAllocationRecord and associated classes of PoolSpace
#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE
#include <impl/Kokkos_SharedAlloc_timpl.hpp>
#undef KOKKOS_IMPL_PUBLIC_INCLUDE
#endif

namespace Kokkos::Impl {
#if defined(KOKKOS_ENABLE_CUDA)
template<>
HostInaccessibleSharedAllocationRecordCommon<
  alps::PoolSpace<Kokkos::CudaSpace>>::
  HostInaccessibleSharedAllocationRecordCommon(
    alps::PoolSpace<Kokkos::CudaSpace> const&         space,
    std::string const&                                label,
    std::size_t                                       alloc_size,
    SharedAllocationRecord<void, void>::function_type dealloc)
  : SharedAllocationRecord<void, void>(
#ifdef KOKKOS_ENABLE_DEBUG
      &s_root_record,
#endif
      checked_allocation_with_header(space, label, alloc_size),
      sizeof(SharedAllocationHeader) + alloc_size,
      dealloc,
      label)
  , m_space(space)
{
  SharedAllocationHeader header;

  fill_host_accessible_header_info(this, header, label);

  // Different from the original implementation, we don't copy the header to
  // device memory here unless the bounds check is enabled
#if defined(KOKKOS_ENABLE_DEBUG_BOUNDS_CHECK)
  using MemorySpace = alps::PoolSpace<Kokkos::CudaSpace>;
  typename MemorySpace::execution_space exec;
  Kokkos::Impl::DeepCopy<MemorySpace, HostSpace>(
    exec,
    SharedAllocationRecord<void, void>::m_alloc_ptr,
    &header,
    sizeof(SharedAllocationHeader));
  exec.fence(std::string("SharedAllocationRecord<Kokkos::")
             + MemorySpace::name()
             + "Space, void>::SharedAllocationRecord(): "
               "fence after copying header from HostSpace");
#endif
}
#endif

#if defined(KOKKOS_ENABLE_HIP)
template<>
HostInaccessibleSharedAllocationRecordCommon<
  alps::PoolSpace<Kokkos::HIPSpace>>::
  HostInaccessibleSharedAllocationRecordCommon(
    alps::PoolSpace<Kokkos::HIPSpace> const&          space,
    std::string const&                                label,
    std::size_t                                       alloc_size,
    SharedAllocationRecord<void, void>::function_type dealloc)
  : SharedAllocationRecord<void, void>(
#ifdef KOKKOS_ENABLE_DEBUG
      &s_root_record,
#endif
      checked_allocation_with_header(space, label, alloc_size),
      sizeof(SharedAllocationHeader) + alloc_size,
      dealloc,
      label)
  , m_space(space)
{
  SharedAllocationHeader header;

  fill_host_accessible_header_info(this, header, label);

  // Different from the original implementation, we don't copy the header to
  // device memory here unless the bounds check is enabled
#if defined(KOKKOS_ENABLE_DEBUG_BOUNDS_CHECK)
  using MemorySpace = alps::PoolSpace<Kokkos::HIPSpace>;
  typename MemorySpace::execution_space exec;
  Kokkos::Impl::DeepCopy<MemorySpace, HostSpace>(
    exec,
    SharedAllocationRecord<void, void>::m_alloc_ptr,
    &header,
    sizeof(SharedAllocationHeader));
  exec.fence(std::string("SharedAllocationRecord<Kokkos::")
             + MemorySpace::name()
             + "Space, void>::SharedAllocationRecord(): "
               "fence after copying header from HostSpace");
#endif
}
#endif
} // namespace Kokkos::Impl

KOKKOS_IMPL_SHARED_ALLOCATION_RECORD_EXPLICIT_INSTANTIATION(
  alps::PoolSpace<Kokkos::HostSpace>);
#if defined(KOKKOS_ENABLE_CUDA)
KOKKOS_IMPL_HOST_INACCESSIBLE_SHARED_ALLOCATION_RECORD_EXPLICIT_INSTANTIATION(
  alps::PoolSpace<Kokkos::CudaSpace>);
KOKKOS_IMPL_SHARED_ALLOCATION_RECORD_EXPLICIT_INSTANTIATION(
  alps::PoolSpace<Kokkos::CudaUVMSpace>);
KOKKOS_IMPL_SHARED_ALLOCATION_RECORD_EXPLICIT_INSTANTIATION(
  alps::PoolSpace<Kokkos::CudaHostPinnedSpace>);
#elif defined(KOKKOS_ENABLE_HIP)
KOKKOS_IMPL_HOST_INACCESSIBLE_SHARED_ALLOCATION_RECORD_EXPLICIT_INSTANTIATION(
  alps::PoolSpace<Kokkos::HIPSpace>);
KOKKOS_IMPL_SHARED_ALLOCATION_RECORD_EXPLICIT_INSTANTIATION(
  alps::PoolSpace<Kokkos::HIPManagedSpace>);
KOKKOS_IMPL_SHARED_ALLOCATION_RECORD_EXPLICIT_INSTANTIATION(
  alps::PoolSpace<Kokkos::HIPHostPinnedSpace>);
#endif
