# This file defines functions to support adding static analyzers to targets

# Clang-tidy related setup
if(${PROJECT_NAME_UPPERCASE}_ENABLE_CLANG_TIDY)
  find_program(
    ${PROJECT_NAME_UPPERCASE}_CLANG_TIDY_BIN
    NAMES "clang-tidy"
    DOC "Path to clang-tidy executable"
  )
  if(NOT ${PROJECT_NAME_UPPERCASE}_CLANG_TIDY_BIN)
    message(SEND_ERROR "Clang-Tidy requested but executable not found.")
  endif()
endif()

function(target_clangtidy_setup target)
  if(${PROJECT_NAME_UPPERCASE}_ENABLE_CLANG_TIDY AND ${PROJECT_NAME_UPPERCASE}_CLANG_TIDY_BIN)
    set_target_properties(
      ${target} PROPERTIES CXX_CLANG_TIDY ${${PROJECT_NAME_UPPERCASE}_CLANG_TIDY_BIN}
    )
  else()
    set_property(TARGET ${target} PROPERTY CXX_CLANG_TIDY)
  endif()
endfunction()

# Cppcheck related setup
if(${PROJECT_NAME_UPPERCASE}_ENABLE_CPPCHECK)
  find_program(
    ${PROJECT_NAME_UPPERCASE}_CPPCHECK_BIN
    NAMES "cppcheck"
    DOC "Path to cppcheck executable"
  )
  if(NOT ${PROJECT_NAME_UPPERCASE}_CPPCHECK_BIN)
    message(SEND_ERROR "Cppcheck requested but executable not found.")
  endif()
endif()

function(target_cppcheck_setup target)
  if(${PROJECT_NAME_UPPERCASE}_ENABLE_CPPCHECK AND ${PROJECT_NAME_UPPERCASE}_CPPCHECK_BIN)
    set_target_properties(
      ${target} PROPERTIES CXX_CPPCHECK ${${PROJECT_NAME_UPPERCASE}_CPPCHECK_BIN}
    )
  else()
    set_property(TARGET ${target} PROPERTY CXX_CPPCHECK)
  endif()
endfunction()

function(target_static_analysis_setup target)
  target_clangtidy_setup(${target})
  target_cppcheck_setup(${target})
endfunction()
