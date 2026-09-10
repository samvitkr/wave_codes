# Define the header and source files (and directories)
set(HEADERS field_statistics.h mesh_statistics.h diagnostics.h)

set(SOURCES
    turbulence_u_statistics.cpp
    scalar_statistics.cpp
    field_statistics.cpp
    mesh_statistics.cpp
    diagnostics.cpp
)

add_library(alps_diagnostics ${HEADERS} ${SOURCES})
target_link_libraries(
  alps_diagnostics
  PUBLIC alps_la alps_solvers_mesh alps_solvers_field
  PRIVATE project_warnings
)
target_formatter_setup(alps_diagnostics)
target_static_analysis_setup(alps_diagnostics)
