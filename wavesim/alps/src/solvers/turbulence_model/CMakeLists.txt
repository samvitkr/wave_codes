# target: alps_solvers_turb_model
add_library(alps_solvers_turb_model STATIC)
target_sources(
  alps_solvers_turb_model
  PRIVATE anisotropic_minimum_dissipation_options.h
          anisotropic_minimum_dissipation_scalar.cpp
          anisotropic_minimum_dissipation.h
          anisotropic_minimum_dissipation.cpp
          constant_schmidt_options.h
          constant_schmidt.cpp
          constant_schmidt.h
          dynamic_smagorinsky_options.h
          dynamic_smagorinsky_scalar.cpp
          dynamic_smagorinsky.cpp
          dynamic_smagorinsky.h
          log_law_wall_model.cpp
          log_law_wall_model.h
          models_options.cpp
          models_options.h
          models.h
          sgs_flux_curve-inl.h
          sgs_flux.cpp
          sgs_flux.h
          sgs_scalar_flux_curve-inl.h
          sgs_scalar_flux.cpp
          sgs_scalar_flux.h
          smagorinsky_options.h
          smagorinsky.cpp
          smagorinsky.h
          wall_models.cpp
          wall_models.h
          detail/dynamic_smagorinsky_functors.cpp
          detail/dynamic_smagorinsky_functors.h
          detail/dynamic_smagorinsky_scalar_functors.cpp
          detail/dynamic_smagorinsky_scalar_functors.h
          detail/nut_report.cpp
          detail/nut_report.h
          detail/smagorinsky_report.cpp
          detail/smagorinsky_report.h
)

if(ALPS_SGS_USE_LOWER_PRECISION)
  target_compile_definitions(alps_solvers_turb_model PUBLIC ALPS_SGS_USE_LOWER_PRECISION=1)
else()
  target_compile_definitions(alps_solvers_turb_model PUBLIC ALPS_SGS_USE_LOWER_PRECISION=0)
endif()
target_link_libraries(
  alps_solvers_turb_model
  PUBLIC alps_src_root alps_common alps_solvers_field alps_solvers_source_terms
  PRIVATE alps_common
          alps_decomp
          alps_spectral
          alps_la
          alps_solvers_operators
          alps_solvers_mesh
          mpipp
          fmt::fmt
          Kokkos::kokkoscore
          project_warnings
)
target_formatter_setup(alps_solvers_turb_model)
target_static_analysis_setup(alps_solvers_turb_model)

# target: alps_solvers_turb_model_fow
add_library(alps_solvers_turb_model_fow STATIC)
target_sources(alps_solvers_turb_model_fow PRIVATE flow_over_wave.cpp)
target_link_libraries(
  alps_solvers_turb_model_fow
  PUBLIC alps_solvers_turb_model alps_src_root
  PRIVATE alps_solvers_field_fow
          alps_common
          alps_decomp
          alps_spectral
          alps_solvers_operators
          alps_solvers_curve_ns_common
          project_warnings
)
target_formatter_setup(alps_solvers_turb_model_fow)
target_static_analysis_setup(alps_solvers_turb_model_fow)

# alps_solvers_turb_model_fs
add_library(alps_solvers_turb_model_fs STATIC)
target_sources(alps_solvers_turb_model_fs PRIVATE free_surface.cpp)
target_link_libraries(
  alps_solvers_turb_model_fs
  PUBLIC alps_solvers_turb_model alps_src_root
  PRIVATE alps_solvers_field_fs
          alps_common
          alps_decomp
          alps_spectral
          alps_solvers_operators
          alps_solvers_curve_ns_common
          project_warnings
)
