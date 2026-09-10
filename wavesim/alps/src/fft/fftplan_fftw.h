#pragma once

#include <fftw3.h>

#include <memory>
#include <type_traits>

namespace alps::fft {
template<class T>
struct FFTWInterface;

template<>
struct FFTWInterface<float>
{
  using handle_t   = fftwf_plan;
  using unique_ptr = std::unique_ptr<std::remove_pointer_t<handle_t>,
                                     decltype(&fftwf_destroy_plan)>;

  using complex_t = fftwf_complex;

  // NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
  static constexpr auto& plan_r2c = fftwf_plan_dft_r2c;
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
  static constexpr auto& plan_c2r = fftwf_plan_dft_c2r;
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
  static constexpr auto& exec_r2c = fftwf_execute_dft_r2c;
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
  static constexpr auto& exec_c2r = fftwf_execute_dft_c2r;
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
  static constexpr auto& destroy_plan = fftwf_destroy_plan;
};

template<>
struct FFTWInterface<double>
{
  using handle_t   = fftw_plan;
  using unique_ptr = std::unique_ptr<std::remove_pointer_t<handle_t>,
                                     decltype(&fftw_destroy_plan)>;

  using complex_t = fftw_complex;

  // NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
  static constexpr auto& plan_r2c = fftw_plan_dft_r2c;
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
  static constexpr auto& plan_c2r = fftw_plan_dft_c2r;
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
  static constexpr auto& exec_r2c = fftw_execute_dft_r2c;
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
  static constexpr auto& exec_c2r = fftw_execute_dft_c2r;
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
  static constexpr auto& destroy_plan = fftw_destroy_plan;
};
} // namespace alps::fft
