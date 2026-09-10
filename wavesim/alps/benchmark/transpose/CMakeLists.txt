add_executable(bench_transpose main.cpp)
target_link_libraries(
  bench_transpose
  alps_transpose
  alps_program_options
  alps_main_wrapper
  project_warnings
)
target_formatter_setup(bench_transpose)
target_static_analysis_setup(bench_transpose)
