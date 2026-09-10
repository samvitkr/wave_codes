# target: alps_solvers_fs
add_library(alps_solvers_fs STATIC)
target_sources(
  alps_solvers_fs
  PRIVATE bc.cpp
          bc.h
          boussinesq.cpp
          coriolis.cpp
          options.cpp
          options.h
          pressure_grad.cpp
          rayleigh_damp.cpp
          scalar_bc.cpp
          scalar.cpp
          solver_correct.cpp
          solver_io.cpp
          solver_rhs.cpp
          solver_uhat.cpp
          solver.cpp
          solver.h
)

target_link_libraries(
  alps_solvers_fs
  PUBLIC alps_src_root
         alps_common
         alps_solvers_field_fs
         alps_solvers_ns_common
         alps_solvers_source_terms
         alps_solvers_turb_model
  PRIVATE alps_common
          alps_decomp
          alps_spectral
          alps_solvers_mesh
          alps_solvers_operators
          alps_solvers_curve_ns_common
          alps_solvers_turb_model_fs
          project_warnings
)
target_formatter_setup(alps_solvers_fs)
target_static_analysis_setup(alps_solvers_fs)
