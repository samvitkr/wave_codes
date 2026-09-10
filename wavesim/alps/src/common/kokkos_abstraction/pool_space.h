#pragma once

#include <Kokkos_Core_fwd.hpp>
#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE
#include <KokkosCore_Config_DeclareBackend.hpp>
#include <Kokkos_Concepts.hpp>
#include <Kokkos_CopyViews.hpp>
#undef KOKKOS_IMPL_PUBLIC_INCLUDE
#endif

#include <string>

namespace alps {
namespace memory {
// forward declaration
template<class BaseSpace>
class DynamicSizePool;
} // namespace memory

/// PoolSpace class wraps a dynamic memory pool for Kokkos use
template<class BaseSpace>
class PoolSpace
{
  static_assert(::Kokkos::is_memory_space<BaseSpace>::value,
                "template arguments must be memory spaces");

  /// Concatenate BaseSpace's name with "_pool"
  template<typename>
  struct Name;
  template<std::size_t... I1>
  struct Name<std::index_sequence<I1...>>
  {
    static constexpr const char
      value[]{BaseSpace::name()[I1]..., '_', 'p', 'o', 'o', 'l', 0};
  };

  using pool_t = ::alps::memory::DynamicSizePool<BaseSpace>;

 public:
  using base_space   = BaseSpace;
  using memory_space = PoolSpace<BaseSpace>;
  using size_type    = typename BaseSpace::size_type;

  using execution_space = typename BaseSpace::execution_space;

  using device_type = Kokkos::Device<execution_space, memory_space>;

  PoolSpace() = default;

  /**\brief  Allocate untracked memory in the space */
  void* allocate(const execution_space& exec_space,
                 size_t                 arg_alloc_size) const;
  void* allocate(const execution_space& exec_space,
                 const char*            arg_label,
                 size_t                 arg_alloc_size,
                 size_t                 arg_logical_size = 0) const;
  void* allocate(size_t arg_alloc_size) const;
  void* allocate(const char* arg_label,
                 size_t      arg_alloc_size,
                 size_t      arg_logical_size = 0) const;

  /**\brief  Deallocate untracked memory in the space */
  void deallocate(void* arg_alloc_ptr, size_t arg_alloc_size) const;
  void deallocate(const char* arg_label,
                  void*       arg_alloc_ptr,
                  size_t      arg_alloc_size,
                  size_t      arg_logical_size = 0) const;

  /**\brief Return Name of the MemorySpace */
  static constexpr const char* name()
  {
    return Name<std::make_index_sequence<std::char_traits<char>::length(
      BaseSpace::name())>>::value;
  }

  static pool_t& get_allocator();

 private:
  template<class>
  friend class PoolSpace;
  friend class Kokkos::Impl::SharedAllocationRecord<memory_space, void>;

  void* impl_allocate(const execution_space&     exec_space,
                      const char*                arg_label,
                      size_t                     arg_alloc_size,
                      size_t                     arg_logical_size = 0,
                      Kokkos::Tools::SpaceHandle arg_handle =
                        Kokkos::Tools::make_space_handle(name())) const;
  void* impl_allocate(const char*                arg_label,
                      size_t                     arg_alloc_size,
                      size_t                     arg_logical_size = 0,
                      Kokkos::Tools::SpaceHandle arg_handle =
                        Kokkos::Tools::make_space_handle(name())) const;
  void  impl_deallocate(const char*                arg_label,
                        void*                      arg_alloc_ptr,
                        size_t                     arg_alloc_size,
                        size_t                     arg_logical_size = 0,
                        Kokkos::Tools::SpaceHandle arg_handle =
                          Kokkos::Tools::make_space_handle(name())) const;
};

template<class>
struct is_pool_space : public std::false_type
{};

template<class S>
struct is_pool_space<PoolSpace<S>> : public std::true_type
{};

template<class S>
inline constexpr bool is_pool_space_v = is_pool_space<S>::value;

template<class Space, class Enable = void>
struct pool_base_space
{};

template<class Space>
struct pool_base_space<
  Space,
  typename std::enable_if_t<Kokkos::is_memory_space_v<Space>
                            && !alps::is_pool_space_v<Space>>>
{
  using type = Space;
};

template<class Space>
struct pool_base_space<
  PoolSpace<Space>,
  typename std::enable_if_t<Kokkos::is_memory_space_v<Space>>>
{
  using type = Space;
};

template<class Space>
using pool_base_space_t = typename pool_base_space<Space>::type;

template<class BaseMemorySpace>
using memory_pool = alps::PoolSpace<BaseMemorySpace>;

using default_memory_pool =
  memory_pool<Kokkos::DefaultExecutionSpace::memory_space>;

using default_host_memory_pool = memory_pool<Kokkos::HostSpace>;

} // namespace alps

//----------------------------------------------------------------------------

namespace Kokkos {

template<class BaseSpace>
struct is_memory_space<alps::PoolSpace<BaseSpace>> : is_memory_space<BaseSpace>
{};

namespace Impl {

#define ALPS_DECLARE_SPACE_ACCESS(MEMSPACE)                      \
  template<class BaseSpace>                                      \
  struct MemorySpaceAccess<alps::PoolSpace<BaseSpace>, MEMSPACE> \
  {                                                              \
    enum                                                         \
    {                                                            \
      assignable = MemorySpaceAccess < BaseSpace,                \
      MEMSPACE > ::assignable,                                   \
    };                                                           \
    enum                                                         \
    {                                                            \
      accessible = MemorySpaceAccess < BaseSpace,                \
      MEMSPACE > ::accessible,                                   \
    };                                                           \
    enum                                                         \
    {                                                            \
      deepcopy = MemorySpaceAccess < BaseSpace,                  \
      MEMSPACE > ::deepcopy                                      \
    };                                                           \
  };                                                             \
  template<class BaseSpace>                                      \
  struct MemorySpaceAccess<MEMSPACE, alps::PoolSpace<BaseSpace>> \
  {                                                              \
    enum                                                         \
    {                                                            \
      assignable = MemorySpaceAccess < MEMSPACE,                 \
      BaseSpace > ::assignable,                                  \
    };                                                           \
    enum                                                         \
    {                                                            \
      accessible = MemorySpaceAccess < MEMSPACE,                 \
      BaseSpace > ::accessible,                                  \
    };                                                           \
    enum                                                         \
    {                                                            \
      deepcopy = MemorySpaceAccess < MEMSPACE,                   \
      BaseSpace > ::deepcopy                                     \
    };                                                           \
  };

ALPS_DECLARE_SPACE_ACCESS(Kokkos::HostSpace)
#if defined(KOKKOS_ENABLE_CUDA)
ALPS_DECLARE_SPACE_ACCESS(Kokkos::CudaSpace)
ALPS_DECLARE_SPACE_ACCESS(Kokkos::CudaHostPinnedSpace)
ALPS_DECLARE_SPACE_ACCESS(Kokkos::CudaUVMSpace)
#elif defined(KOKKOS_ENABLE_HIP)
ALPS_DECLARE_SPACE_ACCESS(Kokkos::HIPSpace)
ALPS_DECLARE_SPACE_ACCESS(Kokkos::HIPHostPinnedSpace)
ALPS_DECLARE_SPACE_ACCESS(Kokkos::HIPManagedSpace)
#endif

#undef ALPS_DECLARE_SPACE_ACCESS

template<class DstSpace, class SrcSpace, class ExecutionSpace>
struct DeepCopy<DstSpace,
                SrcSpace,
                ExecutionSpace,
                std::enable_if_t<::alps::is_pool_space_v<DstSpace>
                                 || ::alps::is_pool_space_v<SrcSpace>>>
  : public DeepCopy<typename ::alps::pool_base_space_t<DstSpace>,
                    typename ::alps::pool_base_space_t<SrcSpace>,
                    ExecutionSpace>
{
  using DeepCopy<typename ::alps::pool_base_space_t<DstSpace>,
                 typename ::alps::pool_base_space_t<SrcSpace>,
                 ExecutionSpace>::DeepCopy;
};
} // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------

#include "pool_space_allocation_record.h"

KOKKOS_IMPL_SHARED_ALLOCATION_SPECIALIZATION(
  alps::PoolSpace<Kokkos::HostSpace>);
#if defined(KOKKOS_ENABLE_CUDA)
KOKKOS_IMPL_HOST_INACCESSIBLE_SHARED_ALLOCATION_SPECIALIZATION(
  alps::PoolSpace<Kokkos::CudaSpace>);
KOKKOS_IMPL_SHARED_ALLOCATION_SPECIALIZATION(
  alps::PoolSpace<Kokkos::CudaUVMSpace>);
KOKKOS_IMPL_SHARED_ALLOCATION_SPECIALIZATION(
  alps::PoolSpace<Kokkos::CudaHostPinnedSpace>);
#elif defined(KOKKOS_ENABLE_HIP)
KOKKOS_IMPL_HOST_INACCESSIBLE_SHARED_ALLOCATION_SPECIALIZATION(
  alps::PoolSpace<Kokkos::HIPSpace>);
KOKKOS_IMPL_SHARED_ALLOCATION_SPECIALIZATION(
  alps::PoolSpace<Kokkos::HIPHostPinnedSpace>);
KOKKOS_IMPL_SHARED_ALLOCATION_SPECIALIZATION(
  alps::PoolSpace<Kokkos::HIPManagedSpace>);
#endif
