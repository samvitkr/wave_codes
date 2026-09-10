# Define the header and source files (and directories)
set(HEADERS
    block_partition.h
    mdcomm.h
    pencil_plan.h
    pencil_transpose.h
    ghost_cell_exchange.h
)

set(SOURCES block_partition.cpp mdcomm.cpp pencil_plan.cpp ghost_cell_exchange.cpp)

add_library(alps_decomp STATIC ${HEADERS} ${SOURCES})
target_link_libraries(
  alps_decomp
  PUBLIC alps_common
         alps_transpose
         Kokkos::kokkoscore
         mpipp
         enum.hpp
  INTERFACE alps_src_root
  PRIVATE alps_real_type fmt::fmt project_warnings
)
target_formatter_setup(alps_decomp)
target_static_analysis_setup(alps_decomp)
