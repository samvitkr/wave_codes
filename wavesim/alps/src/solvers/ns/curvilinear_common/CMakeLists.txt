# target: alps_solvers_curve_ns_common
add_library(alps_solvers_curve_ns_common STATIC)
target_sources(
  alps_solvers_curve_ns_common
  PRIVATE add_boussinesq_buoyancy.cpp
          add_boussinesq_buoyancy.h
          add_coriolis_force.cpp
          add_coriolis_force.h
          add_pressure_grad.cpp
          add_pressure_grad.h
          advection_functor.h
          apply_bc.h
          convection_functor.h
          diffusion_cn_zeta/reduced_common.h
          diffusion_cn_zeta/jacobi.cpp
          diffusion_cn_zeta/jacobi.h
          diffusion_cn_zeta/pcr.cpp
          diffusion_cn_zeta/pcr.h
          diffusion_cn_zeta/tdma.cpp
          diffusion_cn_zeta/tdma.h
          diffusion_cn_zeta_coeff.h
          diffusion_cn_zeta_eqn.cpp
          diffusion_cn_zeta_eqn.h
          pressure.h
          rayleigh_damp-inl.h
          scalar_bc.cpp
          scalar_bc.h
          strain_rate.h
          viscous.h
)

if(ALPS_POISSON_MIXED_PRECISION_IR)
  target_compile_definitions(
    alps_solvers_curve_ns_common PUBLIC ALPS_POISSON_MIXED_PRECISION_IR=1
  )
else()
  target_compile_definitions(
    alps_solvers_curve_ns_common PUBLIC ALPS_POISSON_MIXED_PRECISION_IR=0
  )
endif()
target_link_libraries(
  alps_solvers_curve_ns_common
  # dependencies for public headers, e.g., add_pressure_grad.h
  PUBLIC alps_container alps_real_type alps_solvers_mesh alps_solvers_bc_types
  # dependencies for consumers of header implementations, e.g. advection_functor.h,
  # rayleigh_damp-inl.h, viscous.h, strain_rate.h, pressure.h, apply_bc.h, etc.
  INTERFACE alps_common
            alps_decomp
            alps_spectral
            alps_solvers_ns_common
            alps_solvers_operators
            alps_solvers_operators_curvilinear
            alps_solvers_mesh
            alps_solvers_field
            alps_solvers_poisson
            alps_solvers_source_terms
            Kokkos::kokkoscore
            mpipp
            alps_src_root
  PRIVATE mpipp project_warnings
)
target_formatter_setup(alps_solvers_curve_ns_common)
target_static_analysis_setup(alps_solvers_curve_ns_common)
