#pragma once

#include <cstddef>

#ifdef __cplusplus
extern "C"
{
#endif

  // C-linkage hooks called from vendored vkFFT_CompileKernel.h.
  // On a cache hit, *binary_out is set to a malloc'd buffer the caller must
  // free(), and *binary_size_out to its byte length; returns 1.
  // Returns 0 on miss.
  int alps_jit_cache_lookup(char const*        source,
                            char const* const* opts,
                            int                num_opts,
                            char**             binary_out,
                            size_t*            binary_size_out);

  void alps_jit_cache_store(char const*        source,
                            char const* const* opts,
                            int                num_opts,
                            char const*        binary,
                            size_t             binary_size);

#ifdef __cplusplus
} // extern "C"
#endif
