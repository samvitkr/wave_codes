add_library(mpipp_test_main OBJECT main.cpp)
target_link_libraries(mpipp_test_main PUBLIC Catch2 MPI::MPI_CXX)
target_include_directories(mpipp_test_main PUBLIC ${TEST_ROOT_DIR})

add_executable(rma_window_test rma_window_test.cpp)
target_link_libraries(rma_window_test PRIVATE mpipp_test_main mpipp)
add_test(NAME rma_window_test
         COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS}
                 $<TARGET_FILE:rma_window_test> ${MPIEXEC_POSTFLAGS}
)

add_executable(shared_window_test shared_window_test.cpp)
target_link_libraries(shared_window_test PRIVATE mpipp_test_main mpipp)
add_test(NAME shared_window_test
         COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS}
                 $<TARGET_FILE:shared_window_test> ${MPIEXEC_POSTFLAGS}
)
