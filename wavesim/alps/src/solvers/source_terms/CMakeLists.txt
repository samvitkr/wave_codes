# target: alps_solvers_source_terms
add_library(alps_solvers_source_terms STATIC)
target_sources(
  alps_solvers_source_terms
  PRIVATE body_force_types.h
          boussinesq_options.h
          boussinesq.cpp
          boussinesq.h
          coriolis_options.cpp
          coriolis_options.h
          coriolis.h
          pressure_grad_options.cpp
          pressure_grad_options.h
          pressure_grad.h
          rayleigh_damp_options.cpp
          rayleigh_damp_options.h
          rayleigh_damp.h
          scalar_source_types.h
)

target_link_libraries(
  alps_solvers_source_terms
  PUBLIC alps_src_root fmt::fmt Kokkos::kokkoscore
  PRIVATE alps_common
          alps_decomp
          alps_la
          alps_solvers_mesh
          alps_solvers_field
          fmt::fmt
          mpipp
          project_warnings
)
target_formatter_setup(alps_solvers_source_terms)
target_static_analysis_setup(alps_solvers_source_terms)
