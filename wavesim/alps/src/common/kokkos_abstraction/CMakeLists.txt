set(HEADERS pool_space_allocation_record.h pool_space.h exec_policy.h reducer.h)
set(SOURCES pool_space_allocation_record.cpp pool_space.cpp)
add_library(alps_kokkos_abstraction STATIC ${HEADERS} ${SOURCES})
target_link_libraries(
  alps_kokkos_abstraction
  PUBLIC Kokkos::kokkoscore
  INTERFACE alps_src_root
  PRIVATE alps_memory alps_runtime project_warnings
)
target_formatter_setup(alps_kokkos_abstraction)
