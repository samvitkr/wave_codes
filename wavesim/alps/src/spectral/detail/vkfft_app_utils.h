#pragma once

#include <fft/vkfft.h>
#include <fft/vkfft_structs.h>

#include <memory>
#include <stdexcept>
#include <string>
#include <utility>

namespace alps::spectral::detail {

using VkApp = std::unique_ptr<VkFFTApplication, void (*)(VkFFTApplication*)>;

inline VkApp make_vkfft_app(VkFFTConfiguration config)
{
  auto app = VkApp{new VkFFTApplication(), &deleteVkFFT};

  auto result = initializeVkFFT(app.get(), std::move(config));
  if (result != VKFFT_SUCCESS) {
    std::string error_msg{getVkFFTErrorString(result)};
    throw std::runtime_error("Failed to initialize VkFFT: " + error_msg);
  }

  return app;
}

} // namespace alps::spectral::detail
