#pragma once

#include <unordered_map>

namespace alps::memory {

template<class BaseAllocator>
class AlignedAllocation
{
 public:
  explicit AlignedAllocation(std::size_t alignment) noexcept;

  /// Return an allocation of `size` bytes that is aligned on the configured
  /// alignment boundary.
  [[nodiscard]] void* allocate(std::size_t size);

  /// Deallocate previously aligned allocation
  void deallocate(void* ptr);

  /// Round up `size` bytes to multiples of alignment
  std::size_t round_up(std::size_t size) const noexcept
  {
    return (size + alignment_ - 1u) & -alignment_;
  }

 private:
  std::unordered_map<void*, void*> ptr_records;
  std::size_t                      alignment_;
};

} // namespace alps::memory
