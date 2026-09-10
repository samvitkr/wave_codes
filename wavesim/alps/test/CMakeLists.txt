set(TEST_ROOT_DIR ${CMAKE_CURRENT_SOURCE_DIR})

add_library(catch_test OBJECT catch_main.cpp mpi_reporter.cpp mpi_reporter.h)
target_link_libraries(catch_test PUBLIC alps_common Catch2 fmt::fmt MPI::MPI_CXX)
target_include_directories(catch_test PUBLIC ${TEST_ROOT_DIR})

add_subdirectory(transpose)
add_subdirectory(decomp)
add_subdirectory(mpipp)
add_subdirectory(spectral)
add_subdirectory(solvers)
