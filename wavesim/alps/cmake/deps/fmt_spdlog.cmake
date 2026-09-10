CPMAddPackage(
  NAME fmt
  URL https://github.com/fmtlib/fmt/archive/refs/tags/12.1.0/fmt-12.1.0.tar.gz
  URL_HASH SHA256=ea7de4299689e12b6dddd392f9896f08fb0777ac7168897a244a6d6085043fea
)

set(patch_name "fmt_suppress_nvcc_non-constexpr_warning.patch")
set(patch_file "${CMAKE_CURRENT_LIST_DIR}/${patch_name}")
message(STATUS "Checking patch: ${patch_file}")
# Apply the patch
find_package(Patch REQUIRED)
execute_process(
  COMMAND "${Patch_EXECUTABLE}" --forward -p1 -i "${patch_file}" --reject-file=-
  WORKING_DIRECTORY ${fmt_SOURCE_DIR}
  RESULT_VARIABLE ret
  ERROR_VARIABLE err
  TIMEOUT 5
)
# if the patch was already applied, patch returns a non-zero exit code
if(NOT ret EQUAL 0)
  execute_process(
    COMMAND "${Patch_EXECUTABLE}" --dry-run -R -p1 -i "${patch_file}"
    WORKING_DIRECTORY ${fmt_SOURCE_DIR}
    RESULT_VARIABLE check_ret
    ERROR_VARIABLE check_err
    TIMEOUT 5
  )
  # if the reverse patch can be applied, it means the patch was already applied
  if(NOT check_ret EQUAL 0)
    message(FATAL_ERROR "Failed to apply patch ${patch_name}: ${err}")
  else()
    message(STATUS "Patch ${patch_name} was already applied.")
  endif()
endif()

target_compile_host_only(fmt)
target_compile_definitions(fmt PUBLIC -DFMT_USE_FLOAT128=0) # disable float128 because not all
                                                            # archs support it
# suppress nvcc warnings about unreachable loop in fmt 11.0.2
target_compile_options_nvcc_only(fmt INTERFACE "SHELL:-Xcudafe --diag_suppress=128")

CPMAddPackage(
  NAME spdlog
  URL https://github.com/gabime/spdlog/archive/refs/tags/v1.16.0/spdlog-v1.16.0.tar.gz
  URL_HASH SHA256=8741753e488a78dd0d0024c980e1fb5b5c85888447e309d9cb9d949bdb52aa3e
  OPTIONS "SPDLOG_FMT_EXTERNAL On"
)
target_compile_host_only(spdlog)
