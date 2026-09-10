add_library(alps_solvers_hos STATIC)
target_link_libraries(
  alps_solvers_hos
  PUBLIC alps_src_root
  PRIVATE project_warnings
)

target_sources(
  alps_solvers_hos
  PRIVATE field.cpp
          field.h
          rhs.cpp
          rhs.h
          series.cpp
          series.h
          smooth.cpp
          smooth.h
          solver_rk_stage.cpp
          solver_rk_stage.h
          solver.cpp
          solver.h
          stats.cpp
          stats.h
)

target_link_libraries(
  alps_solvers_hos
  PUBLIC alps_common alps_io Kokkos::kokkoscore std::filesystem
  PRIVATE alps_spectral
          alps_la
          Kokkos::kokkoscore
          Kokkos::kokkosalgorithms
          fmt::fmt
          spdlog::spdlog
          mpipp
)
target_formatter_setup(alps_solvers_hos)
target_static_analysis_setup(alps_solvers_hos)
