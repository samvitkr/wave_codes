add_executable(bench_fft main.cpp)
target_link_libraries(
  bench_fft
  alps_decomp
  alps_spectral
  alps_program_options
  alps_main_wrapper
  project_warnings
)
target_formatter_setup(bench_fft)
target_static_analysis_setup(bench_fft)
