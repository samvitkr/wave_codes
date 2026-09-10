cmake_minimum_required(VERSION 3.21)

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules" ${CMAKE_MODULE_PATH})
if (DEFINED ENV{BC_HOST})
if ($ENV{BC_HOST} STREQUAL "narwhal")
  # On Narwhal, when craype-accel-nvidia70 is loaded, the C compiler automatically links gtl library, which requires linking C++ library and as a result does not work properly
  # message(WARNING "Skipping C compiler check.")
  # set(CMAKE_C_COMPILER_WORKS 1)
  set(CMAKE_C_FLAGS_INIT "-lstdc++")
endif()
endif()

project(alps NONE)
string(TOLOWER ${PROJECT_NAME} PROJECT_NAME_LOWERCASE)
string(TOUPPER ${PROJECT_NAME} PROJECT_NAME_UPPERCASE)

option(${PROJECT_NAME_UPPERCASE}_BUILD_HDF5 "Download and build HDF5" OFF)

option(${PROJECT_NAME_UPPERCASE}_ENABLE_CUDA "Enable the CUDA backend of ${PROJECT_NAME}." OFF)
option(${PROJECT_NAME_UPPERCASE}_MPI_IS_GPU_AWARE "Whether MPI is GPU aware." ON)
option(${PROJECT_NAME_UPPERCASE}_BUILD_EXAMPLES "Build example main programs for ${PROJECT_NAME}" ON)
option(${PROJECT_NAME_UPPERCASE}_BUILD_TESTS "Build tests for ${PROJECT_NAME}" OFF)
option(${PROJECT_NAME_UPPERCASE}_BUILD_BENCHMARKS "Build benchmarks for ${PROJECT_NAME}" OFF)
option(${PROJECT_NAME_UPPERCASE}_ALLOW_OPENMPI_OVERSUBSCRIBE "Allow OpenMPI to oversubscribe CPUs when executing tests." OFF)
option(${PROJECT_NAME_UPPERCASE}_DEFAULT_SINGLE_PRECISION "Set default Real type to float." OFF)
option(${PROJECT_NAME_UPPERCASE}_MIXED_PRECISION "Allow mixed precision." OFF)
option(${PROJECT_NAME_UPPERCASE}_SGS_USE_LOWER_PRECISION "Allow certain SGS computation to use lower precision." ${${PROJECT_NAME_UPPERCASE}_MIXED_PRECISION})
option(${PROJECT_NAME_UPPERCASE}_POISSON_MIXED_PRECISION_IR "Allow iterative Poisson solver to use mixed precision." ON)

option(${PROJECT_NAME_UPPERCASE}_ENABLE_NORMAL_DEBUG_INFO
       "Enable normal debug info (-g). When OFF, uses minimal debug info for backtraces only (-g1)." OFF
)
option(${PROJECT_NAME_UPPERCASE}_ENABLE_DEBUG_FOR_DEPENDENCIES
       "Enable debug info for dependencies." ON
)
option(${PROJECT_NAME_UPPERCASE}_ENABLE_CLANG_TIDY "Enable static analysis with Clang-Tidy." OFF)
option(${PROJECT_NAME_UPPERCASE}_ENABLE_CPPCHECK "Enable static analysis with Cppcheck." OFF)
option(${PROJECT_NAME_UPPERCASE}_ENABLE_CLANG_FORMAT "Adding custom targets to apply clang-format." OFF)
option(${PROJECT_NAME_UPPERCASE}_WARNINGS_AS_ERRORS "Treat compiler warnings as errors." OFF)

if(NOT CMAKE_BUILD_TYPE)
  message(WARNING "No CMAKE_BUILD_TYPE Specified. Default to Release.")
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: Debug, Release, RelWithDebInfo and MinSizeRel." FORCE)
endif()

# Forbid in-source builds to prevent source tree corruption.
if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
  message(FATAL_ERROR "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there.\n")
endif()

# Detect if we are running in a CI environment
if(DEFINED ENV{CI} AND "$ENV{CI}" STREQUAL "true")
  set(${PROJECT_NAME_UPPERCASE}_CI_BUILD TRUE)
  message(STATUS "CI environment detected.")
else()
  set(${PROJECT_NAME_UPPERCASE}_CI_BUILD FALSE)
endif()

include(cmake/CPM.cmake)

include(cmake/CompilerSelection.cmake) # Config CMAKE_CXX_COMPILER
enable_language(CXX C)

include(cmake/CompileOptions.cmake)
include(cmake/CompilerWarnings.cmake)
include(cmake/StaticAnalyzers.cmake)
include(cmake/Formatter.cmake)

include(cmake/Dependencies.cmake)

set_global_debug_flags()

add_subdirectory(src)

if(${PROJECT_NAME_UPPERCASE}_BUILD_EXAMPLES)
  add_subdirectory(examples)
endif()

if(${PROJECT_NAME_UPPERCASE}_BUILD_TESTS)
  enable_testing()
  add_subdirectory(test)
endif()

if(${PROJECT_NAME_UPPERCASE}_BUILD_BENCHMARKS)
  add_subdirectory(benchmark)
endif()
