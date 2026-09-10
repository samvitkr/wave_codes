#pragma once

#define KOKKOS_TR_LAMBDA [&]

// Workaround the false warnings of nvcc compiler about missing returns when "if
// constexpr" is used.
// CUDA 11.3 or later supports __builtin_unreachable(), otherwise a return
// statement is appended.
#if (defined(__NVCC__))
#if (__CUDACC_VER_MAJOR__ >= 11 && __CUDACC_VER_MINOR__ >= 3)
#define ALPS_UNREACHABLE(...) __builtin_unreachable()
#else
#define ALPS_UNREACHABLE(...) return __VA_ARGS__
#endif
#else
#define ALPS_UNREACHABLE(...)
#endif

#define ALPS_EXPORT __attribute__((visibility("default")))

#define ALPS_C_EXPORT extern "C" ALPS_EXPORT
