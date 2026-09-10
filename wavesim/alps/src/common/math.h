//
// Created by xuananqing on 9/27/22.
// Define some math functions not available from Kokkos
//

#pragma once

#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE
#include <Kokkos_MathematicalFunctions.hpp>
#include <Kokkos_MathematicalSpecialFunctions.hpp>
#undef KOKKOS_IMPL_PUBLIC_INCLUDE
#endif

namespace alps {
/// @brief Compute the square of a scalar.
//  For low integer power, this is much faster than pow functions
template<typename T>
[[nodiscard]] KOKKOS_FORCEINLINE_FUNCTION constexpr auto square(T x)
{
  return x * x;
}

/// @brief Compute the cube of a scalar.
//  For low integer power, this is much faster than pow functions
template<typename T>
[[nodiscard]] KOKKOS_FORCEINLINE_FUNCTION constexpr auto cube(T x)
{
  return x * x * x;
}

/// @brief Calculate one over the square root of the sum of squares of two
/// coordinates.
[[nodiscard]] KOKKOS_FORCEINLINE_FUNCTION auto rhypot(float x, float y)
{
#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
  return ::rhypotf(x, y);
#else
  return 1 / std::hypot(x, y);
#endif
}

/// @brief Calculate one over the square root of the sum of squares of two
/// coordinates.
[[nodiscard]] KOKKOS_FORCEINLINE_FUNCTION auto rhypot(double x, double y)
{
#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
  return ::rhypot(x, y);
#else
  return 1 / std::hypot(x, y);
#endif
}

/// @brief Calculate one over the square root of the sum of squares of three
/// coordinates.
[[nodiscard]] KOKKOS_FORCEINLINE_FUNCTION auto
rnorm3d(float x, float y, float z)
{
#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
  return ::rnorm3df(x, y, z);
#else
  return 1 / std::hypot(x, y, z);
#endif
}

/// @brief Calculate one over the square root of the sum of squares of three
/// coordinates.
[[nodiscard]] KOKKOS_FORCEINLINE_FUNCTION auto
rnorm3d(double x, double y, double z)
{
#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
  return ::rnorm3d(x, y, z);
#else
  return 1 / std::hypot(x, y, z);
#endif
}

/// @brief Kahan compensated summation
/**
 * The algorithm significantly reduces the numerical error accumulation when
 * performing many small additions.
 *
 * @param sum Accumulated sum
 * @param x Value to be added to the sum
 * @param error Running error for compensation
 * @return New sum
 */
template<typename T1, typename T2>
[[nodiscard]] auto accumulate_compensated(T1 sum, T2 x, T1& error)
{
  auto corrected_x = (T1)x - error;
  // volatile to prevent compiler optimizations
  volatile auto tmp = sum + corrected_x;
  volatile auto z   = tmp - sum;
  error             = z - corrected_x;
  return tmp;
}
} // namespace alps
