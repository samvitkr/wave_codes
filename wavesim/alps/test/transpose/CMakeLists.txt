add_executable(single_transpose_test single_transpose_test.cpp)
target_link_libraries(single_transpose_test Catch2::Catch2 alps_transpose)
target_formatter_setup(single_transpose_test)
target_static_analysis_setup(single_transpose_test)
add_test(NAME single_transpose COMMAND $<TARGET_FILE:single_transpose_test>)

add_executable(dist_transpose_test dist_transpose_test.cpp)
target_link_libraries(dist_transpose_test catch_test alps_transpose)
target_formatter_setup(dist_transpose_test)
target_static_analysis_setup(dist_transpose_test)
add_test(NAME dist_transpose COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS}
                                     $<TARGET_FILE:dist_transpose_test>
)
