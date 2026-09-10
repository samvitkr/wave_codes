#pragma once

#include <Kokkos_Macros.hpp>

namespace alps {
namespace fft {

/** Class FFTPlan wraps FFT handle and methods. */
template<class T, class TransformType, class Backend>
class FFTPlan;

/** Types of FFT */
class R2C
{};

class C2R
{};

class Mixed // Mainly for vkfft plans that can be used for both R2C and C2R
{};

/** Supported backends for FFT */
class FFTW
{};
template<class T, class TransformType>
using FFTWPlan = FFTPlan<T, TransformType, FFTW>;

#if defined(KOKKOS_ENABLE_CUDA)
class CUFFT
{};
template<class T, class TransformType>
using CUFFTPlan = FFTPlan<T, TransformType, CUFFT>;
#endif

#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
class VKFFT
{};
template<class T>
using VKFFTPlan = FFTPlan<T, Mixed, VKFFT>;
#endif

} // namespace fft
} // namespace alps
