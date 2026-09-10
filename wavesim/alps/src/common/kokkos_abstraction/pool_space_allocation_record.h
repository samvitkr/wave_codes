#pragma once

#include <Kokkos_Core_fwd.hpp>
#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE
#include <KokkosCore_Config_DeclareBackend.hpp>
#undef KOKKOS_IMPL_PUBLIC_INCLUDE
#endif

#include <string>

namespace Kokkos {
namespace Impl {

// A specialization of HostInaccessibleSharedAllocationRecordCommon for
// alps::PoolSpace<MEMORY_SPACE>. This is to avoid copying the header to device,
// which would incur a fence.
#if defined(KOKKOS_ENABLE_CUDA)
template<>
template<class ExecutionSpace>
HostInaccessibleSharedAllocationRecordCommon<
  alps::PoolSpace<Kokkos::CudaSpace>>::
  HostInaccessibleSharedAllocationRecordCommon(
    ExecutionSpace const&                     exec,
    alps::PoolSpace<Kokkos::CudaSpace> const& space,
    std::string const&                        label,
    std::size_t                               alloc_size,
    record_base_t::function_type              dealloc)
  : SharedAllocationRecord<void, void>(
#ifdef KOKKOS_ENABLE_DEBUG
      &s_root_record,
#endif
      checked_allocation_with_header(exec, space, label, alloc_size),
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
  Kokkos::Impl::DeepCopy<alps::PoolSpace<Kokkos::CudaSpace>, HostSpace>(
    exec,
    SharedAllocationRecord<void, void>::m_alloc_ptr,
    &header,
    sizeof(SharedAllocationHeader));
#endif
}
#endif

#if defined(KOKKOS_ENABLE_HIP)
template<>
template<class ExecutionSpace>
HostInaccessibleSharedAllocationRecordCommon<
  alps::PoolSpace<Kokkos::HIPSpace>>::
  HostInaccessibleSharedAllocationRecordCommon(
    ExecutionSpace const&                    exec,
    alps::PoolSpace<Kokkos::HIPSpace> const& space,
    std::string const&                       label,
    std::size_t                              alloc_size,
    record_base_t::function_type             dealloc)
  : SharedAllocationRecord<void, void>(
#ifdef KOKKOS_ENABLE_DEBUG
      &s_root_record,
#endif
      checked_allocation_with_header(exec, space, label, alloc_size),
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
  Kokkos::Impl::DeepCopy<alps::PoolSpace<Kokkos::HIPSpace>, HostSpace>(
    exec,
    SharedAllocationRecord<void, void>::m_alloc_ptr,
    &header,
    sizeof(SharedAllocationHeader));
#endif
}
#endif
} // namespace Impl
} // namespace Kokkos
