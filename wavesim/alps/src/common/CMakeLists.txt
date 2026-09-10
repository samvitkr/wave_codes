add_subdirectory(base)
add_subdirectory(hash)
add_subdirectory(program_options)
add_subdirectory(memory)
add_subdirectory(container)
add_subdirectory(async)
add_subdirectory(device)
add_subdirectory(utils)
add_subdirectory(runtime)
add_subdirectory(kokkos_abstraction)

add_library(alps_real_type INTERFACE)
target_sources(alps_real_type PUBLIC real_type.h)
target_link_libraries(alps_real_type INTERFACE alps_src_root)
if(ALPS_DEFAULT_SINGLE_PRECISION)
  target_compile_definitions(alps_real_type INTERFACE ALPS_REAL_TYPE=float)
else()
  target_compile_definitions(alps_real_type INTERFACE ALPS_REAL_TYPE=double)
endif()
target_formatter_setup(alps_real_type)

add_library(alps_math INTERFACE)
target_sources(alps_math PUBLIC math.h)
target_link_libraries(alps_math INTERFACE alps_src_root)
target_formatter_setup(alps_math)

add_library(alps_common INTERFACE ${HEADERS})
target_link_libraries(
  alps_common
  INTERFACE alps_base
            alps_container
            alps_async
            alps_memory
            alps_hash
            alps_device
            alps_runtime
            alps_kokkos_abstraction
            alps_program_options
            alps_utils
            alps_math
            alps_real_type
)
target_formatter_setup(alps_common)
target_static_analysis_setup(alps_common)
