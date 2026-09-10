#pragma once

#include <cstddef>
#include <cstdlib>
#include <memory>
#include <optional>
#include <string_view>
#include <vector>

namespace alps::utils {

class KernelCacheStoreImpl;

class KernelCacheStore
{
 public:
  explicit KernelCacheStore(std::size_t capacity = 128);

  ~KernelCacheStore();

  void resize(std::size_t capacity);

  [[nodiscard]] std::size_t capacity() const noexcept;

  struct Binary
  {
    struct FreeDeleter
    {
      void operator()(void* ptr) const noexcept
      {
        // NOLINTNEXTLINE(cppcoreguidelines-no-malloc)
        std::free(ptr);
      }
    };

    using DataPtr = std::unique_ptr<char, FreeDeleter>;

    DataPtr     data;
    std::size_t size;
  };

  // Lookup a cached binary and return a deep-copied owned result.
  // Returns std::nullopt on miss.
  [[nodiscard]] std::optional<Binary>
  lookup_binary(std::string_view                     source,
                std::vector<std::string_view> const& options);

  // Store: deep-copies binary into cache.
  void store(std::string_view                     source,
             std::vector<std::string_view> const& options,
             std::vector<char>                    binary);

  // Clear all entries (cache remains active).
  void clear();

  KernelCacheStore(KernelCacheStore const&)            = delete;
  KernelCacheStore& operator=(KernelCacheStore const&) = delete;

 private:
  std::unique_ptr<KernelCacheStoreImpl> impl_;
};

} // namespace alps::utils
