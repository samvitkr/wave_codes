set(TARGET_NAME alps_wind_wave)

add_executable(${TARGET_NAME} wind_wave.cpp)
target_link_libraries(
  ${TARGET_NAME}
  PUBLIC alps_solvers_fow
         alps_diagnostics
         alps_apps_utils
         std::filesystem
         alps_main_wrapper
  PRIVATE project_warnings
)
set_target_properties(${TARGET_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})
target_formatter_setup(${TARGET_NAME})
target_static_analysis_setup(${TARGET_NAME})
