#include "vkfft_jit_hook.h"

#include <common/base/logging.h>
#include <common/utils/kernel_cache.h>

#include <string_view>
#include <vector>

namespace {

alps::utils::KernelCacheStore& jit_kernel_cache_store()
{
  static alps::utils::KernelCacheStore store;
  return store;
}

alps::Logger jit_cache_logger()
{
  static auto logger = alps::get_logger("fft");
  return logger;
}

} // namespace

extern "C"
{

  int alps_jit_cache_lookup(char const*        source,
                            char const* const* opts,
                            int                num_opts,
                            char**             binary_out,
                            size_t*            binary_size_out)
  {
    if (source == nullptr || binary_out == nullptr
        || binary_size_out == nullptr) {
      return 0;
    }

    std::vector<std::string_view> key_opts;
    key_opts.reserve(num_opts > 0 ? static_cast<std::size_t>(num_opts) : 0U);
    if (opts != nullptr) {
      for (int i = 0; i < num_opts; ++i) {
        if (opts[i] != nullptr) {
          key_opts.emplace_back(opts[i]);
        }
      }
    }

    auto cached_binary =
      jit_kernel_cache_store().lookup_binary(source, key_opts);
    if (cached_binary.has_value()) {
      *binary_size_out = cached_binary->size;
      *binary_out      = cached_binary->data.release();

      jit_cache_logger()->debug("VkFFT JIT cache hit (binary_size={} bytes)",
                                *binary_size_out);
      return 1;
    }

    return 0;
  }

  void alps_jit_cache_store(char const*        source,
                            char const* const* opts,
                            int                num_opts,
                            char const*        binary,
                            size_t             binary_size)
  {
    if (source == nullptr || binary == nullptr || binary_size == 0) {
      return;
    }

    std::vector<std::string_view> key_opts;
    key_opts.reserve(num_opts > 0 ? static_cast<std::size_t>(num_opts) : 0U);
    if (opts != nullptr) {
      for (int i = 0; i < num_opts; ++i) {
        if (opts[i] != nullptr) {
          key_opts.emplace_back(opts[i]);
        }
      }
    }

    std::vector<char> binary_copy(binary, binary + binary_size);
    jit_kernel_cache_store().store(source, key_opts, std::move(binary_copy));
  }

} // extern "C"
