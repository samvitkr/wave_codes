set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)
if(${PROJECT_NAME_UPPERCASE}_ENABLE_CUDA)
  set(CMAKE_CUDA_STANDARD 17)
  set(CMAKE_CUDA_EXTENSIONS OFF)
endif()
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Add debug flags for all languages to allow backtraces when not running in CI. Debug symbols
# are omitted in CI builds to reduce build time and artifact size.
macro(set_global_debug_flags)
  if(NOT DEFINED alps_global_debug_flags_set)
    set(alps_global_debug_flags_set FALSE)
  endif()

  if((NOT ${PROJECT_NAME_UPPERCASE}_CI_BUILD) AND (NOT alps_global_debug_flags_set))
    set(local_c_debug_flag "-g")
    set(local_cxx_debug_flag "-g")
    if(NOT ${PROJECT_NAME_UPPERCASE}_ENABLE_NORMAL_DEBUG_INFO)
      # if CXX or C compiler is one of GNU, Clang, IntelLLVM, add -g1
      if(CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang|IntelLLVM)$")
        set(local_cxx_debug_flag "-g1")
      endif()
      if(CMAKE_C_COMPILER_ID MATCHES "^(GNU|Clang|IntelLLVM)$")
        set(local_c_debug_flag "-g1")
      endif()
    endif()
    add_compile_options("$<$<COMPILE_LANGUAGE:CXX>:${local_cxx_debug_flag}>")
    add_compile_options("$<$<COMPILE_LANGUAGE:C>:${local_c_debug_flag}>")
    if(${local_cxx_debug_flag} STREQUAL "-g")
      # no special handling needed for "-g" flag in CUDA
      add_compile_options("$<$<COMPILE_LANGUAGE:CUDA>:${local_cxx_debug_flag}>")
    else()
      # assume nvcc host compiler has the same capabilities as CXX compiler
      add_compile_options(
        "$<$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>:SHELL:-Xcompiler ${local_cxx_debug_flag}>"
      )
      add_compile_options("$<$<COMPILE_LANG_AND_ID:CUDA,Clang>:${local_cxx_debug_flag}>")
    endif()
    unset(local_cxx_debug_flag)
    unset(local_c_debug_flag)

    # add -gz for GNU and IntelLLVM compilers to compress debug info
    add_compile_options("$<$<COMPILE_LANG_AND_ID:C,GNU,IntelLLVM>:-gz>")
    add_compile_options("$<$<COMPILE_LANG_AND_ID:CXX,GNU,IntelLLVM>:-gz>")
    # assume nvcc host compiler is the same as CXX compiler
    add_compile_options(
      "$<$<AND:$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>,$<CXX_COMPILER_ID:GNU,Clang,IntelLLVM>>:SHELL:-Xcompiler -gz>"
    )
    # for Clang, add -Xarch_host -gz to limit -gz to host code only
    add_compile_options("$<$<COMPILE_LANG_AND_ID:C,Clang>:SHELL:-Xarch_host -gz>")
    add_compile_options("$<$<COMPILE_LANG_AND_ID:CXX,Clang>:SHELL:-Xarch_host -gz>")
    add_compile_options("$<$<COMPILE_LANG_AND_ID:CUDA,Clang>:SHELL:-Xarch_host -gz>")

    # for more reliable backtraces, add asynchronous unwind tables
    set(_arg "-fasynchronous-unwind-tables")
    set(_ids "GNU,Clang,IntelLLVM") # compilers that support this flag
    add_compile_options("$<$<COMPILE_LANG_AND_ID:CXX,${_ids}>:${_arg}>")
    add_compile_options("$<$<COMPILE_LANG_AND_ID:C,${_ids}>:${_arg}>")
    # assume nvcc host compiler is the same as CXX compiler
    add_compile_options(
      "$<$<AND:$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>,$<CXX_COMPILER_ID:${_ids}>>:SHELL:-Xcompiler ${_arg}>"
    )
    add_compile_options("$<$<COMPILE_LANG_AND_ID:CUDA,Clang>:${_arg}>")
    unset(_ids)
    unset(_arg)

    # add -rdynamic to executables, which exposes more function symbols to dynamic loaders and
    # debugging tools for improved symbol resolution in backtraces
    set(_arg "-rdynamic")
    set(_ids "GNU,Clang,IntelLLVM") # compilers that support this flag
    set(_is_exe "$<STREQUAL:$<TARGET_PROPERTY:TYPE>,EXECUTABLE>")
    add_link_options("$<$<AND:${_is_exe},$<LINK_LANG_AND_ID:CXX,${_ids}>>:${_arg}>")
    add_link_options("$<$<AND:${_is_exe},$<LINK_LANG_AND_ID:C,${_ids}>>:${_arg}>")
    add_link_options(
      "$<$<AND:${_is_exe},$<LINK_LANG_AND_ID:CUDA,NVIDIA>,$<CXX_COMPILER_ID:${_ids}>>:SHELL:-Xcompiler ${_arg}>"
    )
    add_link_options("$<$<AND:${_is_exe},$<LINK_LANG_AND_ID:CUDA,Clang>>:${_arg}>")
    unset(_ids)
    unset(_arg)
    unset(_is_exe)

    set(alps_global_debug_flags_set TRUE)
  endif()
endmacro()

# Thin archives makes static archives that only link to backing object files, which saves some
# IO and space.
option(${PROJECT_NAME_UPPERCASE}_ENABLE_THIN_ARCHIVES
       "Enable thin archives for static libraries" ON
)
if(${PROJECT_NAME_UPPERCASE}_ENABLE_THIN_ARCHIVES)
  set(thin_archives_supported TRUE)
  if(${PROJECT_NAME_UPPERCASE}_CXX_COMPILER_IS_NVCC_WRAPPER)
    # NVCC device linker does not support thin archives
    set(thin_archives_supported FALSE)
  endif()

  # Determine whether the archiver in use supports thin archives.
  separate_arguments(AR_VERSION_COMMAND UNIX_COMMAND "${CMAKE_AR} -V")
  execute_process(
    COMMAND ${AR_VERSION_COMMAND}
    OUTPUT_VARIABLE AR_VERSION
    RESULT_VARIABLE AR_STATUS
    ERROR_QUIET
  )
  if(AR_STATUS EQUAL 0)
    if(NOT "${AR_VERSION}" MATCHES "^GNU ar |^LLVM ")
      set(thin_archives_supported FALSE)
    endif()
  else()
    set(thin_archives_supported FALSE)
  endif()
  unset(AR_VERSION_COMMAND)
  unset(AR_STATUS)
  unset(AR_VERSION)

  if(thin_archives_supported)
    message(STATUS "Enabling thin archives (static libraries will not be relocatable)")
    set(CMAKE_CXX_ARCHIVE_CREATE "<CMAKE_AR> crT <TARGET> <LINK_FLAGS> <OBJECTS>")
    set(CMAKE_C_ARCHIVE_CREATE "<CMAKE_AR> crT <TARGET> <LINK_FLAGS> <OBJECTS>")
    set(CMAKE_CXX_ARCHIVE_APPEND "<CMAKE_AR> qT <TARGET> <LINK_FLAGS> <OBJECTS>")
    set(CMAKE_C_ARCHIVE_APPEND "<CMAKE_AR> qT <TARGET> <LINK_FLAGS> <OBJECTS>")
  else()
    message(WARNING "Archiver or linker does not support thin archives")
  endif()
endif()

option(${PROJECT_NAME_UPPERCASE}_ENABLE_FATBIN_COMPRESS_ALL
       "Force fatbin to compress all CUDA functions" ON
)
if(${PROJECT_NAME_UPPERCASE}_ENABLE_CUDA)
  if(${PROJECT_NAME_UPPERCASE}_CXX_COMPILER_IS_NVCC_WRAPPER)
    add_compile_options(
      "$<$<OR:$<COMPILE_LANGUAGE:CXX>,$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>>:SHELL:-Xfatbin --compress-all>"
    )
  else()
    set(_args
        "--start-no-unused-arguments -Xcuda-fatbinary --compress-all --end-no-unused-arguments"
    )
    add_compile_options(
      "$<$<OR:$<COMPILE_LANG_AND_ID:CXX,Clang>,$<COMPILE_LANG_AND_ID:CUDA,Clang>>:SHELL:${_args}>"
    )
    unset(_args)
  endif()
endif()
