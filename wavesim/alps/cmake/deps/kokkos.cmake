option(Kokkos_ENABLE_SERIAL "Enable Serial backend for Kokkos." ON)
option(Kokkos_ENABLE_OPENMP "Enable OpenMP backend for Kokkos." ON)
option(Kokkos_ENABLE_CUDA "Enable CUDA backend for Kokkos."
       ${${PROJECT_NAME_UPPERCASE}_ENABLE_CUDA}
)
option(Kokkos_ENABLE_HIP "Enable HIP backend for Kokkos."
       ${${PROJECT_NAME_UPPERCASE}_ENABLE_HIP}
)
# 2023-11-22 Currently, not all UCX versions available support async memory allocation.
option(Kokkos_ENABLE_IMPL_CUDA_MALLOC_ASYNC "Enable asynchronous CUDA memory allocation." OFF)
option(Kokkos_ENABLE_IMPL_HIP_MALLOC_ASYNC "Enable asynchronous HIP memory allocation." OFF)
# Get all arch variables and forward to kokkos
get_cmake_property(_variable_names VARIABLES)
list(
  FILTER
  _variable_names
  INCLUDE
  REGEX
  "^${PROJECT_NAME_UPPERCASE}_ARCH_[A-Za-z0-9_]+"
)
foreach(_arch_var_name ${_variable_names})
  string(REGEX REPLACE "^${PROJECT_NAME_UPPERCASE}" "Kokkos" _kokkos_arch_var_name
                       ${_arch_var_name}
  )
  option(${_kokkos_arch_var_name} "" ${${_arch_var_name}})
endforeach()

CPMAddPackage(
  NAME kokkos
  VERSION 4.7.02
  URL https://github.com/kokkos/kokkos/releases/download/4.7.02/kokkos-4.7.02.tar.gz
  URL_HASH SHA256=a81826ac0a167933d13506bc2a986fb5517038df9abb780fe9bb2c1d4e80803b
  DOWNLOAD_ONLY TRUE
)

if(CPM_PACKAGE_kokkos_VERSION VERSION_LESS "5.0")
  set(_tmp_cmake_cmp0126 ${CMAKE_POLICY_DEFAULT_CMP0126})
  set(CMAKE_POLICY_DEFAULT_CMP0126 OLD) # Kokkos seems to rely on OLD behaviour of CMP0126 to
                                        # correctly detect nvcc_wrapper as of v4.0.0
endif()

add_subdirectory(${kokkos_SOURCE_DIR} ${kokkos_BINARY_DIR})

if(CPM_PACKAGE_kokkos_VERSION VERSION_LESS "5.0")
  set(CMAKE_POLICY_DEFAULT_CMP0126 ${_tmp_cmake_cmp0126})
endif()

get_target_property(_KOKKOS_INCLUDES kokkoscore INTERFACE_INCLUDE_DIRECTORIES)
set_target_properties(
  kokkoscore PROPERTIES INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${_KOKKOS_INCLUDES}"
)
get_target_property(_KOKKOS_INCLUDES kokkoscontainers INTERFACE_INCLUDE_DIRECTORIES)
set_target_properties(
  kokkoscontainers PROPERTIES INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${_KOKKOS_INCLUDES}"
)
