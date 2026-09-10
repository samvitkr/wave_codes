# Define the header and source files (and directories)
set(HEADERS config_table.h cli_options.h)

set(SOURCES config_table.cpp cli_options.cpp)

add_library(alps_program_options STATIC ${HEADERS} ${SOURCES})

target_link_libraries(
  alps_program_options
  PUBLIC tomlplusplus::tomlplusplus std::filesystem
  INTERFACE alps_src_root
  PRIVATE fmt::fmt cxxopts::cxxopts project_warnings
)
target_compile_host_only(alps_program_options)
target_formatter_setup(alps_program_options)
target_static_analysis_setup(alps_program_options)
