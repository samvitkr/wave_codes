add_executable(poisson_test poisson_test.cpp)
target_link_libraries(
  poisson_test
  project_warnings
  catch_test
  alps_decomp
  alps_solvers_mesh
  alps_solvers_poisson
)
target_formatter_setup(poisson_test)
target_static_analysis_setup(poisson_test)
add_test(NAME poisson_test COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS}
                                   $<TARGET_FILE:poisson_test>
)

add_executable(poisson_mpi_test poisson_mpi_test.cpp)
target_link_libraries(
  poisson_mpi_test
  project_warnings
  catch_test
  alps_decomp
  alps_solvers_mesh
  alps_solvers_poisson
)
target_formatter_setup(poisson_mpi_test)
target_static_analysis_setup(poisson_mpi_test)
add_test(NAME poisson_mpi_test COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS}
                                       $<TARGET_FILE:poisson_mpi_test>
)

add_executable(diffusion_cn_zeta_mpi_test diffusion_cn_zeta_mpi_test.cpp)
target_link_libraries(
  diffusion_cn_zeta_mpi_test
  project_warnings
  catch_test
  alps_decomp
  alps_spectral
  alps_solvers_mesh
  alps_solvers_curve_ns_common
)
target_formatter_setup(diffusion_cn_zeta_mpi_test)
target_static_analysis_setup(diffusion_cn_zeta_mpi_test)
add_test(
  NAME cn_zeta_1rank
  COMMAND
    ${CMAKE_COMMAND} -E env OMP_NUM_THREADS=4 OMP_PROC_BIND=false ${MPIEXEC}
    ${MPIEXEC_NUMPROC_FLAG} 1 ${MPIEXEC_PREFLAGS} $<TARGET_FILE:diffusion_cn_zeta_mpi_test>
    "[single-rank]"
)
foreach(i RANGE 2 6)
  add_test(
    NAME cn_zeta_${i}rank
    COMMAND
      ${CMAKE_COMMAND} -E env OMP_NUM_THREADS=4 OMP_PROC_BIND=false ${MPIEXEC}
      ${MPIEXEC_NUMPROC_FLAG} ${i} ${MPIEXEC_PREFLAGS}
      $<TARGET_FILE:diffusion_cn_zeta_mpi_test> "[multi-rank]"
  )
endforeach()
