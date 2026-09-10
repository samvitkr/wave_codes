#include "devices.h"

#include <common/base/logging.h>

#include <Kokkos_Core.hpp>

#include <cstdlib>

#if defined(KOKKOS_ENABLE_CUDA)
#include <cuda_profiler_api.h>
#elif defined(KOKKOS_ENABLE_HIP)
#include <hip/hip_runtime.h>
#endif

namespace alps {

// From Kokkos_Core.cpp
int get_device_count()
{
#if defined(KOKKOS_ENABLE_CUDA)
  int count;
  KOKKOS_IMPL_CUDA_SAFE_CALL(cudaGetDeviceCount(&count));
  return count;
#elif defined(KOKKOS_ENABLE_HIP)
  int count;
  KOKKOS_IMPL_HIP_SAFE_CALL(hipGetDeviceCount(&count));
  return count;
#elif defined(KOKKOS_ENABLE_SYCL)
  return sycl::device::get_devices(sycl::info::device_type::gpu).size();
#elif defined(KOKKOS_ENABLE_OPENACC)
  return acc_get_num_devices(
    Kokkos::Experimental::Impl::OpenACC_Traits::dev_type);
#elif defined(KOKKOS_ENABLE_OPENMPTARGET)
  return omp_get_num_devices();
#else
  Kokkos::abort("implementation bug");
  return -1;
#endif
}

void set_device(int device_id)
{
  if (device_id < 0) return;
#if defined(KOKKOS_ENABLE_CUDA)
  if (auto ret = cudaSetDevice(device_id); ret != cudaSuccess) {
    throw std::runtime_error(std::string("Failed to set device: ")
                             + cudaGetErrorString(ret));
  }
#elif defined(KOKKOS_ENABLE_HIP)
  if (auto ret = hipSetDevice(device_id); ret != hipSuccess) {
    throw std::runtime_error(std::string("Failed to set device: ")
                             + hipGetErrorString(ret));
  }
#endif
}

void start_profiling()
{
#if defined(KOKKOS_ENABLE_CUDA)
  KOKKOS_IMPL_CUDA_SAFE_CALL(cudaProfilerStart());
#endif
}

void stop_profiling()
{
#if defined(KOKKOS_ENABLE_CUDA)
  KOKKOS_IMPL_CUDA_SAFE_CALL(cudaProfilerStop());
#endif
}

namespace detail {
void check_last_device_error(char const* const file, int const line)
{
#if defined(KOKKOS_ENABLE_CUDA)
  cudaError_t const err{cudaGetLastError()};
  if (err != cudaSuccess) {
    spdlog::error("CUDA error captured: {} at {}:{}. The reported source "
                  "location may not be the actual source of the error. The "
                  "error likely occurred before this location.",
                  cudaGetErrorString(err),
                  file,
                  line);
  }
#elif defined(KOKKOS_ENABLE_HIP)
  hipError_t const err{hipGetLastError()};
  if (err != hipSuccess) {
    spdlog::error("HIP error captured: {} at {}:{}. The reported source "
                  "location may not be the actual source of the error. The "
                  "error likely occurred before this location.",
                  hipGetErrorString(err),
                  file,
                  line);
  }
#endif
  (void)file;
  (void)line;
}
} // namespace detail

} // namespace alps
