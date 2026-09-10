# target: alps_solvers_channel
add_library(alps_solvers_channel STATIC)
target_sources(
  alps_solvers_channel
  PRIVATE bc.cpp
          bc.h
          boussinesq.cpp
          convection_diffusion.cpp
          convection_diffusion.h
          coriolis.cpp
          diffusion_cn.cpp
          diffusion_cn.h
          pressure_grad.cpp
          pressure.cpp
          pressure.h
          rayleigh_damp.cpp
          scalar_advection_diffusion.cpp
          scalar_advection_diffusion.h
          scalar_bc.cpp
          scalar.cpp
          solver_ab2_io.cpp
          solver_ab2_uhat.cpp
          solver_ab2.cpp
          solver_ab2.h
          solver_ab2cn_uhat.cpp
          solver_ab2cn.cpp
          solver_ab2cn.h
          solver_projection.cpp
          solver.h
          strain_rate.cpp
          strain_rate.h
)
target_link_libraries(
  alps_solvers_channel
  # solver.h dependencies
  PUBLIC alps_src_root
         alps_base
         alps_solvers_field
         alps_solvers_ns_common
         alps_solvers_source_terms
         alps_solvers_turb_model
  PRIVATE alps_common
          alps_decomp
          alps_spectral
          alps_la
          alps_solvers_mesh
          alps_solvers_operators
          alps_solvers_poisson
          alps_solvers_field
          alps_solvers_source_terms
          alps_solvers_turb_model
          alps_solvers_ns_common
          alps_io
          Kokkos::kokkoscore
          fmt::fmt
          mpipp
          project_warnings
)
target_formatter_setup(alps_solvers_channel)
target_static_analysis_setup(alps_solvers_channel)
