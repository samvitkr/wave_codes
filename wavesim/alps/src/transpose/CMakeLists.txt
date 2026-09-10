add_library(alps_transpose STATIC)

# Define the header and source files (and directories)
target_sources(
  alps_transpose
  PRIVATE transpose.h
          transposer_base.cpp
          transposer_base.h
          transposer_mpi_all2all.cpp
          transposer_mpi_all2all.h
          transposer_mpi_p2p.cpp
          transposer_mpi_p2p.h
          transposer_mpi_p2pshm.cpp
          transposer_mpi_p2pshm.h
          transposer_single_openmp.cpp
          transposer_single.h
          transposer_pool.cpp
          transposer_pool.h
)

if(${PROJECT_NAME_UPPERCASE}_ENABLE_CUDA OR ${PROJECT_NAME_UPPERCASE}_ENABLE_HIP)
  target_sources(alps_transpose PRIVATE transposer_single_cuda.cpp)
endif()

target_link_libraries(
  alps_transpose
  PUBLIC alps_common Kokkos::kokkoscore mpipp enum.hpp
  INTERFACE alps_src_root
  PRIVATE nonstd::span-lite fmt::fmt project_warnings
)
target_formatter_setup(alps_transpose)
target_static_analysis_setup(alps_transpose)
