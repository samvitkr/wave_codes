add_library(alps_solvers_mesh_fwd INTERFACE)
target_link_libraries(alps_solvers_mesh_fwd INTERFACE alps_src_root)

# Define the header and source files (and directories)
set(HEADERS mesh_fwd.h mesh.h curvilinear_mesh.h)

set(SOURCES mesh.cpp curvilinear_mesh.cpp)

add_library(alps_solvers_mesh STATIC ${HEADERS} ${SOURCES})
target_link_libraries(
  alps_solvers_mesh
  PUBLIC alps_solvers_mesh_fwd
         alps_spectral_fwd
         alps_container
         alps_math
         alps_real_type
         alps_io
         Kokkos::kokkoscore
         nonstd::span-lite
  INTERFACE alps_src_root
  PRIVATE alps_kokkos_abstraction
          alps_spectral
          fmt::fmt
          mpipp
          project_warnings
)
target_formatter_setup(alps_solvers_mesh)
target_static_analysis_setup(alps_solvers_mesh)
