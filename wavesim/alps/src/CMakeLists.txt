add_library(${PROJECT_NAME}_src_root INTERFACE)
target_include_directories(${PROJECT_NAME}_src_root INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})

add_subdirectory(mpipp)
add_subdirectory(common)
add_subdirectory(transpose)
add_subdirectory(decomp)
add_subdirectory(io)
add_subdirectory(fft)
add_subdirectory(spectral)
add_subdirectory(linear_algebra)

add_subdirectory(solvers)
add_subdirectory(diagnostics)

add_subdirectory(apps)
