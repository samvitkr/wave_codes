add_library(alps_solvers_poisson STATIC)

# public interface (base class header)
target_sources(alps_solvers_poisson PRIVATE tridiagonal_solver.h)
target_link_libraries(
  alps_solvers_poisson PUBLIC alps_src_root alps_base Kokkos::kokkoscore mpipp
)

# implementation
target_sources(
  alps_solvers_poisson
  PRIVATE detail/tridiagonal_wang_backward_alltoall.cpp
          detail/tridiagonal_wang_tuning.cpp
          tridiagonal_mpi_lu.cpp
          tridiagonal_mpi_lu.h
          tridiagonal_single_cusparse.cpp
          tridiagonal_single_cusparse.h
          tridiagonal_single.cpp
          tridiagonal_single.h
          tridiagonal_solver.cpp
          tridiagonal_wang.cpp
          tridiagonal_wang.h
)
target_link_libraries(
  alps_solvers_poisson
  PRIVATE alps_common
          alps_decomp
          Kokkos::kokkoscore
          mpipp
          project_warnings
)
if(${PROJECT_NAME_UPPERCASE}_ENABLE_CUDA)
  target_link_libraries(alps_solvers_poisson PRIVATE CUDA::cusparse)
elseif(${PROJECT_NAME_UPPERCASE}_ENABLE_HIP)
  find_package(rocsparse REQUIRED CONFIG)
  target_link_libraries(alps_solvers_poisson PRIVATE roc::rocsparse)
endif()

target_formatter_setup(alps_solvers_poisson)
target_static_analysis_setup(alps_solvers_poisson)
