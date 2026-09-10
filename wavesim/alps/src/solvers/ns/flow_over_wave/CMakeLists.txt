# target: alps_solvers_fow
add_library(alps_solvers_fow STATIC)
target_sources(
  alps_solvers_fow
  PRIVATE bc.cpp
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
          solver_uhat_ab2cn.cpp
          solver.cpp
          solver.h
)

target_link_libraries(
  alps_solvers_fow
  PUBLIC alps_src_root
         alps_common
         alps_solvers_field_fow
         alps_solvers_ns_common
         alps_solvers_source_terms
         alps_solvers_turb_model
  PRIVATE alps_common
          alps_decomp
          alps_spectral
          alps_solvers_mesh
          alps_solvers_operators
          alps_solvers_curve_ns_common
          alps_solvers_turb_model_fow
          project_warnings
)
target_formatter_setup(alps_solvers_fow)
target_static_analysis_setup(alps_solvers_fow)
