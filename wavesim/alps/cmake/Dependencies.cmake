find_package(MPI COMPONENTS CXX)
find_package(OpenMP REQUIRED)

find_package(
  Filesystem
  COMPONENTS Final
  REQUIRED
)

find_package(FFTW COMPONENTS double float)

if(${PROJECT_NAME_UPPERCASE}_ENABLE_CUDA)
  find_package(CUDAToolkit REQUIRED)
  if(NOT TARGET CUDA::cufft)
    message(FATAL_ERROR "CUFFT not found.")
  endif()
  if(NOT TARGET CUDA::nvrtc)
    message(FATAL_ERROR "NVRTC not found.")
  endif()
endif()

if(${PROJECT_NAME_UPPERCASE}_BUILD_TESTS)
  CPMAddPackage(
    NAME catch2
    URL https://github.com/catchorg/Catch2/archive/refs/tags/v3.11.0/Catch2-v3.11.0.tar.gz
  )
  if(TARGET Catch2)
    target_compile_host_only(Catch2)
  endif()
  if(TARGET Catch2WithMain)
    target_compile_host_only(Catch2WithMain)
  endif()
endif()

include(${CMAKE_CURRENT_LIST_DIR}/deps/hdf5.cmake)

# Enable debug flags for subsequent dependencies (called after catch2 and HDF5 to exclude them)
if(${PROJECT_NAME_UPPERCASE}_ENABLE_DEBUG_FOR_DEPENDENCIES)
  set_global_debug_flags()
endif()

include(${CMAKE_CURRENT_LIST_DIR}/deps/kokkos.cmake)

include(${CMAKE_CURRENT_LIST_DIR}/deps/vkfft.cmake)

include(${CMAKE_CURRENT_LIST_DIR}/deps/fmt_spdlog.cmake)

CPMAddPackage(
  NAME HighFive
  URL https://github.com/highfive-devs/highfive/archive/refs/tags/v3.2.0/highfive-v3.2.0.tar.gz
  URL_HASH SHA256=01ea2eed7dbce1cf5dfff59476cfa113a7822b641aecbd99c674592fe7a4e630
  DOWNLOAD_ONLY TRUE
)
add_library(HighFive INTERFACE)
target_include_directories(HighFive INTERFACE ${HighFive_SOURCE_DIR}/include)
target_link_libraries(HighFive INTERFACE HDF5::HDF5 MPI::MPI_CXX)
target_compile_definitions(HighFive INTERFACE HIGHFIVE_CXX_STD=${CMAKE_CXX_STANDARD})

CPMAddPackage(
  NAME cxxopts
  GITHUB_REPOSITORY jarro2783/cxxopts
  VERSION 3.3.1
  OPTIONS "CXXOPTS_BUILD_EXAMPLES Off" "CXXOPTS_BUILD_TESTS Off" "CXXOPTS_ENABLE_INSTALL Off"
          "CXXOPTS_ENABLE_WARNINGS Off"
)

CPMAddPackage(
  NAME tomlplusplus
  GITHUB_REPOSITORY marzer/tomlplusplus
  GIT_TAG v3.4.0
)

add_subdirectory(${PROJECT_SOURCE_DIR}/thirdparty/span-lite)
add_subdirectory(${PROJECT_SOURCE_DIR}/thirdparty/enum)

set(GIT_FAIL_IF_NONZERO_EXIT false)
add_subdirectory(${PROJECT_SOURCE_DIR}/thirdparty/cmake-git-version-tracking)

# A sanity check to ensure that -ccbin is not overridden by dependencies
if(${PROJECT_NAME_UPPERCASE}_CXX_COMPILER_IS_NVCC_WRAPPER AND (NOT "$CACHE{CMAKE_CXX_FLAGS}"
                                                               MATCHES "-ccbin")
)
  message(WARNING "-ccbin in CMAKE_CXX_FLAGS was overridden by dependencies.")
endif()
