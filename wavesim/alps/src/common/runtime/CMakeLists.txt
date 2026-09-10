set(HEADERS manager.h)
set(SOURCES manager.cpp)
add_library(alps_runtime STATIC ${HEADERS} ${SOURCES})
target_link_libraries(
  alps_runtime
  PUBLIC Kokkos::kokkoscore alps_memory alps_async
  INTERFACE alps_src_root
  PRIVATE alps_base
          alps_device
          alps_utils
          mpipp
          project_warnings
)
if(ALPS_ENABLE_CUDA)
  target_link_libraries(alps_runtime PRIVATE CUDA::cudart)
endif()
target_formatter_setup(alps_runtime)
target_static_analysis_setup(alps_runtime)
