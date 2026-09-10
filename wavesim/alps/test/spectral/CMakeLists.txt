add_executable(spectral_test spectral_test.cpp)
target_link_libraries(
  spectral_test
  catch_test
  alps_decomp
  alps_spectral
  Kokkos::kokkosalgorithms
)
target_formatter_setup(spectral_test)
target_static_analysis_setup(spectral_test)
add_test(NAME spectral_test COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS}
                                    $<TARGET_FILE:spectral_test>
)

if(ALPS_ENABLE_CUDA OR ALPS_ENABLE_HIP)
  add_executable(vkfft_batched_1d_conv_test vkfft_batched_1d_conv_test.cpp)
  target_link_libraries(
    vkfft_batched_1d_conv_test
    Catch2
    Kokkos::kokkoscore
    fmt::fmt
    VkFFT
    alps_spectral
  )
  if(ALPS_ENABLE_HIP)
    find_package(rocthrust REQUIRED CONFIG)
    target_link_libraries(vkfft_batched_1d_conv_test roc::rocthrust)
  endif()
  target_include_directories(vkfft_batched_1d_conv_test PRIVATE ${CMAKE_SOURCE_DIR}/src)
  target_compile_definitions(vkfft_batched_1d_conv_test PRIVATE ALPS_JIT_HOOK)
  target_formatter_setup(vkfft_batched_1d_conv_test)
  target_static_analysis_setup(vkfft_batched_1d_conv_test)
  add_test(NAME vkfft_batched_1d_conv_test COMMAND $<TARGET_FILE:vkfft_batched_1d_conv_test>)
endif()

if(ALPS_ENABLE_CUDA OR ALPS_ENABLE_HIP)
  add_executable(vkfft_tuning_test vkfft_tuning_test.cpp)
  target_link_libraries(
    vkfft_tuning_test
    Catch2
    Kokkos::kokkoscore
    fmt::fmt
    VkFFT
    alps_spectral
  )
  if(ALPS_ENABLE_HIP)
    find_package(rocthrust REQUIRED CONFIG)
    target_link_libraries(vkfft_tuning_test roc::rocthrust)
  endif()
  target_include_directories(vkfft_tuning_test PRIVATE ${CMAKE_SOURCE_DIR}/src)
  target_compile_definitions(vkfft_tuning_test PRIVATE ALPS_JIT_HOOK)
  target_formatter_setup(vkfft_tuning_test)
  add_test(NAME vkfft_tuning_test COMMAND $<TARGET_FILE:vkfft_tuning_test>)
endif()
