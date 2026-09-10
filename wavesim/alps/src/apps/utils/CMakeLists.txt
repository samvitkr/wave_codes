# Define the header and source files (and directories)
set(HEADERS std_cli_opts.h symlinks.h progress_tracker.h git_version.h)

set(SOURCES std_cli_opts.cpp symlinks.cpp progress_tracker.cpp git_version.cpp)

add_library(alps_apps_utils ${HEADERS} ${SOURCES})
target_link_libraries(
  alps_apps_utils
  PUBLIC alps_program_options std::filesystem
  INTERFACE alps_src_root
  PRIVATE cmake_git_version_tracking fmt::fmt project_warnings
)
target_compile_host_only(alps_apps_utils)
target_formatter_setup(alps_apps_utils)
target_static_analysis_setup(alps_apps_utils)
