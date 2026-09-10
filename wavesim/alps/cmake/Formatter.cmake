# This file defines support adding clang-format as custom targets

if(${PROJECT_NAME_UPPERCASE}_ENABLE_CLANG_FORMAT)
  find_program(
    ${PROJECT_NAME_UPPERCASE}_CLANG_FORMAT_BIN
    NAMES clang-format
    DOC "Path to clang-format executable"
  )
  if(NOT ${PROJECT_NAME_UPPERCASE}_CLANG_FORMAT_BIN)
    message(SEND_ERROR "Clang-format requested but executable not found")
  endif()
endif()

# Copyright Tomas Zeman 2019. Distributed under the Boost Software License, Version 1.0. (See
# accompanying file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

# Add a custom ${target_prefix}_clangformat to execute clang-format on all sources passed as
# arguments following target_prefix Also add ${target}_clangformat as the dependency of the
# aggregate 'clangformat' target
function(clangformat_setup target_prefix)
  unset(clangformat_sources)

  foreach(clangformat_source ${ARGN})
    if(clangformat_source MATCHES "^\\$<")
      continue()
    endif()

    get_filename_component(clangformat_source ${clangformat_source} ABSOLUTE)
    if(EXISTS ${clangformat_source} AND NOT IS_DIRECTORY ${clangformat_source})
      list(APPEND clangformat_sources ${clangformat_source})
    endif()
  endforeach()

  if(NOT clangformat_sources)
    return()
  endif()

  add_custom_target(
    ${target_prefix}_clangformat
    COMMAND ${${PROJECT_NAME_UPPERCASE}_CLANG_FORMAT_BIN} -style=file -i ${clangformat_sources}
    COMMENT "Formatting with clang-format..."
    VERBATIM
  )

  if(TARGET clangformat)
    add_dependencies(clangformat ${target_prefix}_clangformat)
  else()
    add_custom_target(clangformat DEPENDS ${target_prefix}_clangformat)
  endif()
endfunction()

function(target_clangformat_setup target)
  if(${PROJECT_NAME_UPPERCASE}_ENABLE_CLANG_FORMAT
     AND ${PROJECT_NAME_UPPERCASE}_CLANG_FORMAT_BIN
  )
    get_target_property(target_sources ${target} SOURCES)
    if(target_sources AND NOT target_sources STREQUAL "target_sources-NOTFOUND")
      clangformat_setup(${target} ${target_sources})
    endif()
  endif()
endfunction()

function(target_formatter_setup target)
  target_clangformat_setup(${target})
endfunction()
