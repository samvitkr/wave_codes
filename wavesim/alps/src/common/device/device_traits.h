#pragma once

#include <Kokkos_Core_fwd.hpp>

#include <type_traits>

namespace alps {
/// checks if Kokkos is configured with a heterogeneous execution space
struct has_device
  : std::negation<std::is_same<Kokkos::DefaultExecutionSpace,
                               Kokkos::DefaultHostExecutionSpace>>
{};

/// if Kokkos is configured with a heterogeneous execution space
inline constexpr auto has_device_v = has_device::value;

/// checks if an execution space is the default execution space
template<class ExecSpace>
struct is_default_execution_space
  : std::is_same<Kokkos::DefaultExecutionSpace, std::decay_t<ExecSpace>>
{};

/// if an execution space is the default execution space
template<class ExecSpace>
inline constexpr auto is_default_execution_space_v =
  is_default_execution_space<ExecSpace>::value;

/// checks if an execution space is the host execution space
template<class ExecSpace>
struct is_host_execution_space
  : std::is_same<Kokkos::DefaultHostExecutionSpace, std::decay_t<ExecSpace>>
{};

/// if an execution space is the host execution space
template<class ExecSpace>
inline constexpr auto is_host_execution_space_v =
  is_host_execution_space<ExecSpace>::value;

template<class ExecSpace>
struct is_openmp_execution_space
  : std::is_same<Kokkos::OpenMP, std::decay_t<ExecSpace>>
{};
template<class ExecSpace>
inline constexpr auto is_openmp_execution_space_v =
  is_openmp_execution_space<ExecSpace>::value;

#if defined(KOKKOS_ENABLE_CUDA)
template<class ExecSpace>
struct is_cuda_execution_space
  : std::is_same<Kokkos::Cuda, std::decay_t<ExecSpace>>
{};
#else
template<class ExecSpace>
struct is_cuda_execution_space : std::false_type
{};
#endif
template<class ExecSpace>
inline constexpr auto is_cuda_execution_space_v =
  is_cuda_execution_space<ExecSpace>::value;

#if defined(KOKKOS_ENABLE_HIP)
template<class ExecSpace>
struct is_hip_execution_space
  : std::is_same<Kokkos::HIP, std::decay_t<ExecSpace>>
{};
#else
template<class ExecSpace>
struct is_hip_execution_space : std::false_type
{};
#endif
template<class ExecSpace>
inline constexpr auto is_hip_execution_space_v =
  is_hip_execution_space<ExecSpace>::value;

#if (ALPS_MPI_IS_GPU_AWARE != 0)
struct mpi_is_gpu_aware : std::true_type
{};
#else
struct mpi_is_gpu_aware : std::false_type
{};
#endif
inline constexpr auto mpi_is_gpu_aware_v = mpi_is_gpu_aware::value;

template<class MemorySpace>
struct mpi_can_access : std::false_type
{};

template<>
struct mpi_can_access<Kokkos::HostSpace> : std::true_type
{};

#if defined(KOKKOS_ENABLE_CUDA)
template<>
struct mpi_can_access<Kokkos::CudaSpace> : mpi_is_gpu_aware
{};
#elif defined(KOKKOS_ENABLE_HIP)
template<>
struct mpi_can_access<Kokkos::HIPSpace> : mpi_is_gpu_aware
{};
#endif

// Kokkos::SharedSpace and Kokkos::SharedHostPinnedSpace point to
// Kokkos::HostSpace when compiled for host only therefore the following
// definitions need to be guareded by KOKKOS_ENABLE_CUDA or KOKKOS_ENABLE_HIP
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
template<>
struct mpi_can_access<Kokkos::SharedHostPinnedSpace> : std::true_type
{};

template<>
struct mpi_can_access<Kokkos::SharedSpace> : std::true_type
{};
#endif

template<class BaseSpace>
class PoolSpace;

template<class BaseSpace>
struct mpi_can_access<PoolSpace<BaseSpace>> : mpi_can_access<BaseSpace>
{};

template<class MemorySpace>
inline constexpr auto mpi_can_access_v = mpi_can_access<MemorySpace>::value;

} // namespace alps
