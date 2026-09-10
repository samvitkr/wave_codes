# Detect C++ compiler
include(CheckLanguage)
check_language(CXX)

if(CMAKE_CXX_COMPILER)
  # Preliminary check if the compiler is an executable
  execute_process(
    COMMAND "${CMAKE_CXX_COMPILER}" --version
    RESULT_VARIABLE _ret
    OUTPUT_VARIABLE _out
    ERROR_VARIABLE _err
  )
  if(NOT _ret EQUAL 0)
    message(
      FATAL_ERROR "C++ compiler '${CMAKE_CXX_COMPILER}' errors when invoked with --version.\n"
                  "Output: ${_out}\nError: ${_err}"
    )
  endif()
else()
  message(FATAL_ERROR "Failed to find a working C++ compiler.")
endif()

function(alps_identify_nvcc_wrapper varName)
  # Check if the --nvcc-wrapper-show option is supported
  execute_process(
    COMMAND "${ARGN}" --nvcc-wrapper-show
    RESULT_VARIABLE _ret
    OUTPUT_QUIET ERROR_QUIET
  )
  # cmake-format: off
  if(NOT _ret EQUAL 0)
    set(${varName} false PARENT_SCOPE)
  else()
    set(${varName} true PARENT_SCOPE)
  endif()
  # cmake-format: on
endfunction()

macro(alps_enable_cuda_language)
  # Clang 17/18 produces errors about the unsupported float128 type when used as the CUDA
  # compiler. Workaround by explicitly adding a -std=c++11 flag before enabling CUDA.
  set(_append_flag FALSE)
  set(_std_flag "-std=c++11")
  get_property(_enabled GLOBAL PROPERTY ENABLED_LANGUAGES)
  if(NOT _enabled MATCHES "CUDA")
    # Determine if a "-std=" or "-std " flag is already set in CMAKE_CUDA_FLAGS
    set(_flag_regex "-std[= ]")
    # Where to append the standard flag depends on whether a cached CMAKE_CUDA_FLAGS exists
    if(DEFINED CACHE{CMAKE_CUDA_FLAGS})
      # set local CMAKE_CUDA_FLAGS to influence enable_language. Overrides any existing
      # CMAKE_CUDA_FLAGS with the cached value
      if(NOT CMAKE_CUDA_FLAGS MATCHES "${_flag_regex}")
        set(CMAKE_CUDA_FLAGS "$CACHE{CMAKE_CUDA_FLAGS} ${_std_flag}")
        set(_append_flag TRUE)
      endif()
    elseif(NOT ("${CMAKE_CUDA_FLAGS_INIT}" MATCHES "${_flag_regex}"
                OR ("$ENV{CUDAFLAGS}" MATCHES "${_flag_regex}"))
    )
      set(CMAKE_CUDA_FLAGS_INIT "${CMAKE_CUDA_FLAGS_INIT} ${_std_flag}")
      set(_append_flag TRUE)
    endif()
    unset(_flag_regex)

    if(NOT CMAKE_CUDA_HOST_COMPILER)
      set(CMAKE_CUDA_HOST_COMPILER "${CMAKE_CXX_COMPILER}")
    endif()
  endif()
  enable_language(CUDA)
  if(_append_flag)
    # Remove the appended flag to avoid side effects
    set(CMAKE_CUDA_FLAGS "") # ensure unsetting the CMAKE_CUDA_FLAGS in the local scope
    string(REPLACE "${_std_flag}" "" CMAKE_CUDA_FLAGS "$CACHE{CMAKE_CUDA_FLAGS}")
    set(CMAKE_CUDA_FLAGS
        "${CMAKE_CUDA_FLAGS}"
        CACHE STRING "CUDA compile flags" FORCE
    )
    string(REPLACE "${_std_flag}" "" CMAKE_CUDA_FLAGS_INIT "${CMAKE_CUDA_FLAGS_INIT}")
  endif()
  unset(_append_flag)
  unset(_std_flag)
endmacro()

macro(alps_generate_nvcc_wrapper_script)
  if(NOT CMAKE_CUDA_HOST_COMPILER)
    set(_cudahostcxx "${CMAKE_CXX_COMPILER}")
  else()
    set(_cudahostcxx "${CMAKE_CUDA_HOST_COMPILER}")
  endif()
  CPMAddPackage(
    NAME nvcc_wrapper
    URL https://github.com/kokkos/kokkos/raw/5.0.0/bin/nvcc_wrapper
    VERSION 5.0.0
    DOWNLOAD_ONLY TRUE
    DOWNLOAD_NO_EXTRACT TRUE
  )

  set(_script_dir "${PROJECT_BINARY_DIR}/CMakeFiles")
  set(_script "${_script_dir}/nvcc_wrapper")
  # cmake-format: off
  file(
    COPY "${nvcc_wrapper_SOURCE_DIR}/nvcc_wrapper"
    DESTINATION "${_script_dir}"
    FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
  )
  # cmake-format: on

  set(patch_file ${_script_dir}/nvcc_wrapper.patch)
  configure_file(${CMAKE_CURRENT_LIST_DIR}/deps/nvcc_wrapper.patch.in ${patch_file} @ONLY)

  find_package(Patch REQUIRED)
  execute_process(
    COMMAND ${Patch_EXECUTABLE} -p1 -i ${patch_file}
    WORKING_DIRECTORY ${_script_dir}
    RESULT_VARIABLE ret
    ERROR_VARIABLE err
    TIMEOUT 5
  )
  if(NOT ret EQUAL 0)
    message(FATAL_ERROR "Failed to patch nvcc_wrapper: ${err}")
  endif()
  unset(patch_file)
  unset(ret)
  unset(err)

  set(${PROJECT_NAME_UPPERCASE}_GENERATED_NVCC_WRAPPER
      "${_script}"
      CACHE INTERNAL "Path to the generated nvcc_wrapper script"
  )
  set(${PROJECT_NAME_UPPERCASE}_NVCC_WRAPPER_HOST_COMPILER
      "${_cudahostcxx}"
      CACHE INTERNAL "Host compiler used by the generated nvcc_wrapper script"
  )
  message(STATUS "Generated nvcc_wrapper script: ${_script}")
  message(STATUS "nvcc_wrapper invokes:")
  message(STATUS "  NVCC compiler: ${CMAKE_CUDA_COMPILER}")
  message(STATUS "  Host compiler: ${_cudahostcxx}")
  unset(_cudahostcxx)
  unset(_script)
  unset(_script_dir)
endmacro()

# CUDA
if(${PROJECT_NAME_UPPERCASE}_ENABLE_CUDA
   AND (NOT DEFINED CACHE{${PROJECT_NAME_UPPERCASE}_CXX_COMPILER_IS_NVCC_WRAPPER})
)
  alps_identify_nvcc_wrapper(_cxx_is_nvcc_wrapper ${CMAKE_CXX_COMPILER})
  set(${PROJECT_NAME_UPPERCASE}_CXX_COMPILER_IS_NVCC_WRAPPER
      ${_cxx_is_nvcc_wrapper}
      CACHE INTERNAL "NVCC wrapper script is used as CXX compiler"
  )
  unset(_cxx_is_nvcc_wrapper)
  if(${PROJECT_NAME_UPPERCASE}_CXX_COMPILER_IS_NVCC_WRAPPER
     AND (NOT DEFINED ${PROJECT_NAME_UPPERCASE}_GENERATED_NVCC_WRAPPER)
  )
    # Got here because nvcc_wrapper is specified as CXX compiler by the user.
    message(FATAL_ERROR "nvcc_wrapper is not intended to be set directly.")
  endif()

  # the detected CXX compiler is probably a host compiler, let CMake determine
  alps_enable_cuda_language()

  # A workaround to add missing CMAKE_CUDA_COMPILER_TOOLKIT_VERSION for CMake < 3.23
  if(NOT CMAKE_CUDA_COMPILER_TOOLKIT_VERSION)
    if(CMAKE_CUDA_COMPILER_ID STREQUAL "Clang")
      execute_process(
        COMMAND ${_CUDA_NVCC_EXECUTABLE} "--version"
        OUTPUT_VARIABLE CMAKE_CUDA_COMPILER_ID_OUTPUT
      )
    endif()

    if(CMAKE_CUDA_COMPILER_ID_OUTPUT MATCHES [=[V([0-9]+\.[0-9]+\.[0-9]+)]=])
      set(CMAKE_CUDA_COMPILER_TOOLKIT_VERSION "${CMAKE_MATCH_1}")
    endif()
  endif()

  if(CMAKE_CUDA_COMPILER_ID STREQUAL "NVIDIA")
    alps_generate_nvcc_wrapper_script()
    set(${PROJECT_NAME_UPPERCASE}_CXX_COMPILER_IS_NVCC_WRAPPER
        TRUE
        CACHE INTERNAL "NVCC wrapper script is used as CXX compiler"
    )

    # Prepend -ccbin to CMAKE_CXX_FLAGS to ensure nvcc_wrapper uses the correct host compiler
    # Although the host compiler is already specified in the nvcc_wrapper script, adding -ccbin
    # here helps ccache to pick up the host compiler.
    if("${CMAKE_CXX_FLAGS_INIT} ${CMAKE_CXX_FLAGS} $ENV{CXXFLAGS}" MATCHES "-ccbin[= ]")
      message(
        FATAL_ERROR
          "-ccbin flag should not be set in CMAKE_CXX_FLAGS, CMAKE_CXX_FLAGS_INIT or CXXFLAGS environment variable."
      )
    endif()
    set(_prepend_flag "-ccbin ${${PROJECT_NAME_UPPERCASE}_NVCC_WRAPPER_HOST_COMPILER}")
    if(DEFINED CACHE{CMAKE_CXX_FLAGS})
      # If user has set CMAKE_CXX_FLAGS in cache, prepend -ccbin to it
      if(NOT "$CACHE{CMAKE_CXX_FLAGS}" STREQUAL "${CMAKE_CXX_FLAGS}")
        message(
          WARNING
            "CMAKE_CXX_FLAGS has been modified from the cached value and is being overridden."
        )
      endif()
      set(CMAKE_CXX_FLAGS "${_prepend_flag} $CACHE{CMAKE_CXX_FLAGS}")
      set(CMAKE_CXX_FLAGS
          "${CMAKE_CXX_FLAGS}"
          CACHE STRING "CXX compile flags" FORCE
      )
    else()
      # otherwise prepend -ccbin to CMAKE_CXX_FLAGS_INIT
      set(CMAKE_CXX_FLAGS_INIT "${_prepend_flag} ${CMAKE_CXX_FLAGS_INIT}")
    endif()
    unset(_prepend_flag)

    set(CMAKE_CXX_COMPILER "${${PROJECT_NAME_UPPERCASE}_GENERATED_NVCC_WRAPPER}")
    set(CMAKE_CXX_COMPILER
        "${${PROJECT_NAME_UPPERCASE}_GENERATED_NVCC_WRAPPER}"
        CACHE STRING "CXX compiler" FORCE
    )
    message(STATUS "Replace CXX compiler using nvcc_wrapper.")
  endif()
endif()
if(NOT DEFINED CACHE{${PROJECT_NAME_UPPERCASE}_CXX_COMPILER_IS_NVCC_WRAPPER})
  set(${PROJECT_NAME_UPPERCASE}_CXX_COMPILER_IS_NVCC_WRAPPER
      FALSE
      CACHE INTERNAL "NVCC wrapper script is used as CXX compiler"
  )
endif()

include(${CMAKE_CURRENT_LIST_DIR}/NVCCWrapperUtils.cmake)
