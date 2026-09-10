add_executable(mdcomm_test mdcomm_test.cpp)
target_link_libraries(mdcomm_test catch_test alps_decomp)
# target_include_directories(mdcomm_test PRIVATE ${TEST_ROOT_DIR})
target_formatter_setup(mdcomm_test)
target_static_analysis_setup(mdcomm_test)
add_test(NAME mdcomm_test COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS}
                                  $<TARGET_FILE:mdcomm_test>
)

add_executable(blockpartition_test blockpartition_test.cpp)
target_link_libraries(blockpartition_test catch_test alps_decomp)
# target_include_directories(mdcomm_test PRIVATE ${TEST_ROOT_DIR})
target_formatter_setup(blockpartition_test)
target_static_analysis_setup(blockpartition_test)
add_test(NAME blockpartition_test
         COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS}
                 $<TARGET_FILE:blockpartition_test>
)

add_executable(pencilplan_test pencilplan_test.cpp)
target_link_libraries(pencilplan_test catch_test alps_decomp)
# target_include_directories(dtranspose_test PRIVATE ${TEST_ROOT_DIR})
target_formatter_setup(pencilplan_test)
target_static_analysis_setup(pencilplan_test)
add_test(NAME pencilplan_test COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS}
                                      $<TARGET_FILE:pencilplan_test>
)
