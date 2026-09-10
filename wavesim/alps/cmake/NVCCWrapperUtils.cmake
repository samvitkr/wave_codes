function(target_compile_host_only target)
  if(${PROJECT_NAME_UPPERCASE}_ENABLE_CUDA
     AND ${PROJECT_NAME_UPPERCASE}_CXX_COMPILER_IS_NVCC_WRAPPER
  )
    target_compile_options(${target} PRIVATE --host-only)
  endif()
endfunction()

# Usage: source_compile_host_only(source1 [source2 ...])
function(source_compile_host_only)
  if(${PROJECT_NAME_UPPERCASE}_ENABLE_CUDA
     AND ${PROJECT_NAME_UPPERCASE}_CXX_COMPILER_IS_NVCC_WRAPPER
  )
    foreach(source_file IN LISTS ARGN)
      get_source_file_property(old_compile_options ${source_file} COMPILE_OPTIONS)
      if(NOT old_compile_options)
        set(old_compile_options)
      endif()

      list(FIND old_compile_options --host-only host_only_index)
      if(host_only_index EQUAL -1)
        set_source_files_properties(
          ${source_file} PROPERTIES COMPILE_OPTIONS "${old_compile_options};--host-only"
        )
      endif()
    endforeach()
  endif()
endfunction()

macro(target_compile_options_nvcc_only)
  if(${PROJECT_NAME_UPPERCASE}_ENABLE_CUDA
     AND ${PROJECT_NAME_UPPERCASE}_CXX_COMPILER_IS_NVCC_WRAPPER
  )
    target_compile_options(${ARGV})
  endif()
endmacro()
