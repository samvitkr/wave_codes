# =================================================================================
add_library(alps_solvers_bc_types OBJECT)
target_sources(
  alps_solvers_bc_types PRIVATE bc_types.cpp bc_types.h scalar_bc_types.cpp scalar_bc_types.h
)
target_link_libraries(
  alps_solvers_bc_types
  PUBLIC alps_program_options alps_real_type Kokkos::kokkoscore
  INTERFACE alps_src_root
  PRIVATE fmt::fmt project_warnings
)
target_formatter_setup(alps_solvers_bc_types)
target_static_analysis_setup(alps_solvers_bc_types)

# =================================================================================
# target: alps_solvers_field
add_library(alps_solvers_field STATIC)
target_sources(
  alps_solvers_field
  PRIVATE $<TARGET_OBJECTS:alps_solvers_bc_types>
          flow_field.cpp
          flow_field.h
          scalar_fields.cpp
          scalar_fields.h
          traits.h
)
target_link_libraries(
  alps_solvers_field
  PUBLIC alps_src_root
         alps_common
         alps_io
         alps_solvers_mesh
         alps_solvers_bc_types
  INTERFACE alps_src_root
  PRIVATE project_warnings
)
target_formatter_setup(alps_solvers_field)
target_static_analysis_setup(alps_solvers_field)

# =================================================================================
# target: alps_solvers_field_fow
add_library(alps_solvers_field_fow STATIC)
target_sources(
  alps_solvers_field_fow
  PRIVATE bc_waves.cpp
          bc_waves.h
          prescribed_waves.cpp
          prescribed_waves.h
          flow_over_wave_field.cpp
          flow_over_wave_field.h
)
target_link_libraries(
  alps_solvers_field_fow
  PUBLIC alps_solvers_field
         alps_solvers_bc_types
         alps_solvers_mesh
         alps_common
         alps_io
  INTERFACE alps_src_root
  PRIVATE alps_common alps_decomp fmt::fmt project_warnings
)
target_formatter_setup(alps_solvers_field_fow)
target_static_analysis_setup(alps_solvers_field_fow)

# =================================================================================
# target: alps_solvers_field_fs
add_library(alps_solvers_field_fs STATIC)
target_sources(
  alps_solvers_field_fs PRIVATE free_surface_flow_field.cpp free_surface_flow_field.h
)
target_link_libraries(
  alps_solvers_field_fs
  PUBLIC alps_solvers_field
         alps_solvers_bc_types
         alps_solvers_mesh
         alps_common
         alps_io
  INTERFACE alps_src_root
  PRIVATE alps_common alps_decomp alps_solvers_fs project_warnings
)
target_formatter_setup(alps_solvers_field_fs)
target_static_analysis_setup(alps_solvers_field_fs)
