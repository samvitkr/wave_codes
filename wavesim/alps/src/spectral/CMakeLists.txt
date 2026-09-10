add_library(alps_spectral_fwd INTERFACE)
target_sources(alps_spectral_fwd PRIVATE spectral_fwd.h)
target_link_libraries(alps_spectral_fwd INTERFACE alps_fft_fwd alps_src_root)
target_formatter_setup(alps_spectral_fwd)

add_library(alps_spectral STATIC)
target_link_libraries(
  alps_spectral
  INTERFACE alps_src_root
  PRIVATE project_warnings
)

# public interface (base class)
target_sources(alps_spectral PRIVATE spectral.h spectral.cpp spectral_base.h spectral_base.cpp)
target_link_libraries(
  alps_spectral
  PUBLIC alps_spectral_fwd
         alps_fft_fwd
         alps_base
         alps_device
         alps_decomp
         Kokkos::kokkoscore
)

# backend implementations (derived classes)
target_sources(
  alps_spectral PRIVATE spectral_host.h spectral_host.cpp detail/spectral_host_tuning.cpp
)
target_link_libraries(
  alps_spectral
  PRIVATE alps_base
          alps_utils
          alps_kokkos_abstraction
          alps_decomp
          alps_fft
          Kokkos::kokkoscore
)

if(${PROJECT_NAME_UPPERCASE}_ENABLE_CUDA)
  target_sources(alps_spectral PRIVATE spectral_cuda.h spectral_cuda.cpp)
  target_link_libraries(
    alps_spectral
    PRIVATE alps_base
            alps_device
            alps_kokkos_abstraction
            alps_decomp
            alps_fft
            Kokkos::kokkoscore
  )
endif()
if(${PROJECT_NAME_UPPERCASE}_ENABLE_CUDA OR ${PROJECT_NAME_UPPERCASE}_ENABLE_HIP)
  target_sources(
    alps_spectral
    PRIVATE spectral_cuda_vkfft.h
            spectral_cuda_vkfft.cpp
            detail/vkfft_tuning.h
            detail/vkfft_tuning.cpp
            detail/vkfft_app_utils.h
            detail/vkfft_jit_hook.h
            detail/vkfft_jit_hook.cpp
  )
  target_link_libraries(
    alps_spectral
    PRIVATE alps_base
            alps_hash
            alps_utils
            alps_kokkos_abstraction
            alps_decomp
            alps_fft
            Kokkos::kokkoscore
  )
  target_compile_definitions(alps_spectral PRIVATE ALPS_JIT_HOOK)
endif()

target_formatter_setup(alps_spectral)
target_static_analysis_setup(alps_spectral)
