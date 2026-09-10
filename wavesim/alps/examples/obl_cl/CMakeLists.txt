set(TARGET_NAME alps_obl)

set(HEADERS vortex_force.h translate_velocity.h)
set(SRC vortex_force.cpp translate_velocity.cpp)

add_executable(${TARGET_NAME} obl.cpp ${HEADERS} ${SRC})
target_link_libraries(
  ${TARGET_NAME}
  PUBLIC alps_solvers_channel
         alps_diagnostics
         alps_apps_utils
         std::filesystem
         alps_main_wrapper
  PRIVATE project_warnings
)
set_target_properties(${TARGET_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})
target_formatter_setup(${TARGET_NAME})
target_static_analysis_setup(${TARGET_NAME})
