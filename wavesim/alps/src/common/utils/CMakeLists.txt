# Define the header and source files (and directories)
set(HEADERS to_lower_case.h timeout_mpi_barrier.h bench.h kernel_cache.h)

set(SOURCES to_lower_case.cpp timeout_mpi_barrier.cpp bench.cpp kernel_cache.cpp)

add_library(alps_utils STATIC ${HEADERS} ${SOURCES})

target_include_directories(alps_utils PRIVATE ${PROJECT_SOURCE_DIR}/thirdparty/xxHash)

source_compile_host_only(to_lower_case.cpp timeout_mpi_barrier.cpp bench.cpp kernel_cache.cpp)

target_link_libraries(
  alps_utils
  PUBLIC mpipp
  INTERFACE alps_src_root
  PRIVATE alps_base project_warnings
)
# target_compile_host_only(alps_utils)
target_formatter_setup(alps_utils)
target_static_analysis_setup(alps_utils)
