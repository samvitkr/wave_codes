# =================================================================================
add_library(alps_solvers_operators STATIC)
target_sources(alps_solvers_operators PRIVATE div.cpp div.h grad.cpp grad.h)
target_link_libraries(
  alps_solvers_operators
  PUBLIC alps_src_root alps_solvers_mesh alps_container Kokkos::kokkoscore
  PRIVATE alps_common alps_decomp alps_spectral project_warnings
)
target_formatter_setup(alps_solvers_operators)
target_static_analysis_setup(alps_solvers_operators)

# =================================================================================
add_library(alps_solvers_operators_curvilinear STATIC)
target_sources(
  alps_solvers_operators_curvilinear
  PRIVATE div_curvilinear.cpp
          div_curvilinear.h
          grad_curvilinear.cpp
          grad_curvilinear.h
          lap_curvilinear.cpp
          lap_curvilinear.h
)
target_link_libraries(
  alps_solvers_operators_curvilinear
  PUBLIC alps_src_root alps_solvers_mesh alps_container
  PRIVATE alps_solvers_operators
          alps_common
          alps_decomp
          alps_spectral
          project_warnings
)
target_formatter_setup(alps_solvers_operators_curvilinear)
target_static_analysis_setup(alps_solvers_operators_curvilinear)
