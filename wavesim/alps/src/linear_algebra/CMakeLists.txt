# Define the header and source files (and directories)
set(HEADERS
    nrminf3d.h
    axpy.h
    sum.h
    mean_variance.h
    covariance.h
)

set(SOURCES nrminf3d.cpp mean_variance.cpp covariance.cpp)

add_library(alps_la ${HEADERS} ${SOURCES})
target_link_libraries(
  alps_la
  PUBLIC alps_base
         alps_container
         alps_kokkos_abstraction
         Kokkos::kokkoscore
         mpipp
  INTERFACE alps_src_root
  PRIVATE alps_real_type project_warnings
)
target_formatter_setup(alps_la)
target_static_analysis_setup(alps_la)
