//
// Description: This file is used to suppress warnings from vkFFT library.
// Created by xuanx004 on 4/26/24.
//

#pragma once

#if !defined(pfLD)
#define ALPS_DEFINE_VKFFT_TYPES
#define VKFFT_MAX_FFT_DIMENSIONS 4
#define pfLD long double
#define pfUINT uint64_t
#define pfINT int64_t
#endif
#include <vkFFT/vkFFT_Structs/vkFFT_Structs.h>
#if defined(ALPS_DEFINE_VKFFT_TYPES)
#undef ALPS_DEFINE_VKFFT_TYPES
#undef VKFFT_MAX_FFT_DIMENSIONS
#undef pfLD
#undef pfUINT
#undef pfINT
#endif
