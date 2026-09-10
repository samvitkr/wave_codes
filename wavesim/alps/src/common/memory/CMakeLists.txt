set(HEADERS dynamic_size_pool.h aligned_allocation.h pooled_memory_resource.h)
set(SOURCES aligned_allocation.cpp dynamic_size_pool.cpp pooled_memory_resource.cpp)
add_library(alps_memory STATIC ${HEADERS} ${SOURCES})
target_link_libraries(
  alps_memory
  PUBLIC Kokkos::kokkoscore alps_base
  INTERFACE alps_src_root
  PRIVATE project_warnings
)
target_formatter_setup(alps_memory)
