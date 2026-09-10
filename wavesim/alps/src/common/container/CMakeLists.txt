set(HEADERS view_types.h view_utils.h vector_field.h matrix_field.h)
set(SOURCES dummy.cpp)
add_library(alps_container STATIC ${HEADERS} ${SOURCES})
target_link_libraries(
  alps_container
  PUBLIC Kokkos::kokkoscore Kokkos::kokkoscontainers
  INTERFACE alps_src_root
  PRIVATE alps_kokkos_abstraction project_warnings
)
target_formatter_setup(alps_container)
