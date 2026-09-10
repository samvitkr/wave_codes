set(HEADERS macros.h logging_fwd.h logging.h)
set(SOURCES logging.cpp)
add_library(alps_base STATIC ${HEADERS} ${SOURCES})
target_compile_host_only(alps_base)
target_link_libraries(
  alps_base
  PUBLIC spdlog::spdlog
  INTERFACE alps_src_root
  PRIVATE mpipp project_warnings
)
target_formatter_setup(alps_base)
