# target: alps_solvers_ns_common
add_library(alps_solvers_ns_common STATIC)
target_sources(
  alps_solvers_ns_common
  PRIVATE advection_invoke.h
          channel_options.cpp
          channel_options.h
          convection_invoke.h
          pressure_eqn_coeffs.h
          validate_max_div.h
)
target_link_libraries(
  alps_solvers_ns_common
  # public and self dependencies: channel_options.h channel_options.cpp
  PUBLIC alps_program_options alps_solvers_source_terms alps_solvers_turb_model
  # dependency for consumers of header implementations: validate_max_div.h etc.
  INTERFACE alps_common
            alps_decomp
            alps_spectral
            alps_solvers_mesh
            Kokkos::kokkoscore
            mpipp
            alps_src_root
  PRIVATE project_warnings
)
target_formatter_setup(alps_solvers_ns_common)
target_static_analysis_setup(alps_solvers_ns_common)

add_subdirectory(channel)
add_subdirectory(curvilinear_common)
add_subdirectory(flow_over_wave)
add_subdirectory(free_surface)
