set(TARGET_NAME alps_hos)

add_executable(${TARGET_NAME} hos.cpp)
target_link_libraries(
  ${TARGET_NAME}
  PRIVATE alps_solvers_hos
          alps_common
          alps_apps_utils
          std::filesystem
          alps_main_wrapper
          mpipp
          project_warnings
)
set_target_properties(${TARGET_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})
target_formatter_setup(${TARGET_NAME})
target_static_analysis_setup(${TARGET_NAME})
