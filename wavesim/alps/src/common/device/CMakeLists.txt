set(HEADERS devices.h device_traits.h)
set(SOURCES devices.cpp)
add_library(alps_device STATIC ${HEADERS} ${SOURCES})
target_link_libraries(
  alps_device
  PUBLIC Kokkos::kokkoscore
  INTERFACE alps_src_root
  PRIVATE alps_base project_warnings
)
if(ALPS_MPI_IS_GPU_AWARE)
  target_compile_definitions(alps_device PUBLIC ALPS_MPI_IS_GPU_AWARE=1)
else()
  target_compile_definitions(alps_device PUBLIC ALPS_MPI_IS_GPU_AWARE=0)
endif()
if(ALPS_ENABLE_CUDA)
  target_link_libraries(alps_device PRIVATE CUDA::cudart)
endif()
# roctx tracing (for ranges and markers) needs linking to roctx library we list roctx library as
# a dependency here to ensure it is linked to all downstream targets
if(ALPS_ENABLE_HIP)
  # As of ROCm 6.2, rocprofiler-sdk-roctx supercedes roctracer.
  find_package(rocprofiler-sdk-roctx CONFIG PATHS $ENV{ROCM_PATH})
  if(rocprofiler-sdk-roctx_FOUND)
    target_link_libraries(alps_device PUBLIC rocprofiler-sdk-roctx::rocprofiler-sdk-roctx)
  else()
    find_library(ROCM_ROCTX_LIB roctx64 HINTS $ENV{ROCM_PATH}/lib)
    if(ROCM_ROCTX_LIB)
      target_link_libraries(alps_device PUBLIC ${ROCM_ROCTX_LIB})
    endif()
  endif()
endif()
target_formatter_setup(alps_device)
