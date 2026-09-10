set(_VKFFT_SOURCE_DIR ${PROJECT_SOURCE_DIR}/thirdparty/VkFFT)

add_library(VkFFT INTERFACE)
target_include_directories(VkFFT INTERFACE ${_VKFFT_SOURCE_DIR}/vkFFT)
target_include_directories(VkFFT INTERFACE ${_VKFFT_SOURCE_DIR}/half_lib)

if(${PROJECT_NAME_UPPERCASE}_ENABLE_CUDA)
  target_compile_definitions(VkFFT INTERFACE -DVKFFT_BACKEND=1)
  target_link_libraries(VkFFT INTERFACE CUDA::nvrtc)
endif()

if(${PROJECT_NAME_UPPERCASE}_ENABLE_HIP)
  target_compile_definitions(VkFFT INTERFACE -DVKFFT_BACKEND=2)
  find_package(hiprtc REQUIRED CONFIG)
  target_link_libraries(VkFFT INTERFACE hiprtc::hiprtc)
endif()
