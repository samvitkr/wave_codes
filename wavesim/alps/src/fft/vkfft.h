//
// Description: This file is used to suppress warnings from vkFFT library.
// Created by xuanx004 on 4/26/24.
//

#pragma once

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wdouble-promotion"
#pragma GCC diagnostic ignored "-Wsign-compare"
#pragma GCC diagnostic ignored "-Wcomment"
#ifdef __NVCC_DIAG_PRAGMA_SUPPORT__
#pragma nv_diagnostic push
#pragma nv_diag_suppress 68
#pragma nv_diag_suppress 177
#pragma nv_diag_suppress 550
#else
#pragma diagnostic push
#pragma diag_suppress 68
#pragma diag_suppress 177
#pragma diag_suppress 550
#endif
#include <vkFFT.h>
#ifdef __NVCC_DIAG_PRAGMA_SUPPORT__
#pragma nv_diagnostic pop
#else
#pragma diagnostic pop
#endif
#pragma GCC diagnostic pop
