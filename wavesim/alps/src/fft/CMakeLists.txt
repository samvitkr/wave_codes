# Define the header and source files (and directories)
add_library(alps_fft_fwd INTERFACE)
target_sources(alps_fft_fwd PRIVATE fftplan_fwd.h)
target_link_libraries(alps_fft_fwd INTERFACE Kokkos::kokkoscore alps_src_root)
target_formatter_setup(alps_fft_fwd)

add_library(alps_fft STATIC)
target_link_libraries(
  alps_fft
  INTERFACE alps_fft_fwd alps_src_root
  PRIVATE project_warnings
)

# Wrapper header for FFTW, header-only interface
target_sources(alps_fft PRIVATE fftplan_fftw.h dummy.cpp)
target_link_libraries(alps_fft INTERFACE FFTW::double FFTW::float alps_src_root)

if(${PROJECT_NAME_UPPERCASE}_ENABLE_CUDA)
  target_sources(alps_fft PRIVATE fftplan_cufft.h fftplan_cufft.cpp vkfft_structs.h vkfft.h)
  target_link_libraries(
    alps_fft
    PUBLIC alps_base CUDA::cufft
    INTERFACE VkFFT
    PRIVATE alps_memory alps_runtime fmt::fmt
  )
endif()
if(${PROJECT_NAME_UPPERCASE}_ENABLE_HIP)
  target_link_libraries(alps_fft INTERFACE VkFFT)
endif()

target_formatter_setup(alps_fft)
target_static_analysis_setup(alps_fft)
