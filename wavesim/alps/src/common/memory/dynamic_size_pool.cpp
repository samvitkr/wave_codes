//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2016-22, Lawrence Livermore National Security, LLC and Umpire
// project contributors. See the COPYRIGHT file for details.
//
// SPDX-License-Identifier: (MIT)
//////////////////////////////////////////////////////////////////////////////
#include "dynamic_size_pool.h"

#include <common/base/logging.h>

#include <Kokkos_Core_fwd.hpp>

#include <string>

namespace alps::memory {

using namespace std::string_literals;

template<class BaseSpace>
DynamicSizePool<BaseSpace>& DynamicSizePool<BaseSpace>::global_instance()
{
  constexpr std::size_t Alignment =
    512u; // an alignment consistent with recent Cuda behaviour

  static DynamicSizePool instance{Alignment};
  return instance;
}

template<class BaseSpace>
DynamicSizePool<BaseSpace>::~DynamicSizePool()
{
  if (!usedBlocks.empty()) {
    logger_->warn("Memory pool is not empty at destruction");
  }
  freeReleasedBlocks();
  logger_->trace("Memory pool destroyed");
}

template<class BaseSpace>
void DynamicSizePool<BaseSpace>::release()
{
  std::lock_guard lock(mutex_);
  freeReleasedBlocks();
}

template<class BaseSpace>
[[nodiscard]] void* DynamicSizePool<BaseSpace>::allocate(std::size_t bytes)
{
  if (bytes == 0u) {
    return &zero_byte_data;
  }

  logger_->trace("Try allocating {:}B", bytes);

  std::lock_guard lock(mutex_);

  return impl_allocate(bytes);
}

template<class BaseSpace>
[[nodiscard]] void* DynamicSizePool<BaseSpace>::impl_allocate(std::size_t bytes)
{
  const std::size_t rounded_bytes = allocator_.round_up(bytes);

  auto best_block = findUsableBlock(rounded_bytes);

  // Allocate a block if needed
  if (best_block == freeBlocks.cend()) {
    best_block = allocateBlock(rounded_bytes);
    if (best_block == freeBlocks.cend()) {
      // Block allocation failed
      return nullptr;
    }
  }

  // Split the needed bytes from the memory block and add it to the used
  // list
  current_size_ += rounded_bytes;
  usedBlocks.push_back(splitBlock(best_block, rounded_bytes));
  const auto& new_block = usedBlocks.back();

  // Return the new pointer
  logger_->trace("Spliced a block of size {:}B at {}",
                 rounded_bytes,
                 fmt::ptr(new_block.data));
  return new_block.data;
}

template<class BaseSpace>
void DynamicSizePool<BaseSpace>::deallocate(void* ptr)
{
  logger_->trace("deallocate ptr={}", fmt::ptr(ptr));
  if (ptr == nullptr || ptr == &zero_byte_data) return;

  std::lock_guard lock(mutex_);
  impl_deallocate(ptr);

  if (auto suggested_size = should_coalesce(); suggested_size != 0) {
    logger_->debug("Coalesce heuristic satisfied.");
    coalesce(suggested_size);
  }
}

template<class BaseSpace>
void DynamicSizePool<BaseSpace>::impl_deallocate(void* ptr)
{
  // Find the associated block (reverse search)
  auto r_it = usedBlocks.crbegin();
  while (r_it != usedBlocks.crend()) {
    if (r_it->data == ptr) break;
    ++r_it;
  }
  if (r_it == usedBlocks.crend()) {
    logger_->warn("Invalid ptr ({}) to deallocate", fmt::ptr(ptr));
    return;
  }

  r_it++;
  auto it = r_it.base();

  // Release it
  current_size_ -= it->size;
  releaseBlock(*it);
  usedBlocks.erase(it);
}

template<class BaseSpace>
std::size_t DynamicSizePool<BaseSpace>::getActualSize() const noexcept
{
  std::lock_guard lock(mutex_);
  return actual_bytes_;
}

template<class BaseSpace>
std::size_t DynamicSizePool<BaseSpace>::getCurrentSize() const noexcept
{
  std::lock_guard lock(mutex_);
  return current_size_;
}

template<class BaseSpace>
std::size_t DynamicSizePool<BaseSpace>::getActualHighwaterMark() const noexcept
{
  std::lock_guard lock(mutex_);
  return actual_highwatermark_;
}

template<class BaseSpace>
std::size_t
DynamicSizePool<BaseSpace>::getLargestAvailableBlock() const noexcept
{
  std::size_t largest_block{0};

  std::lock_guard lock(mutex_);
  for (const auto& block : freeBlocks) {
    if (block.size > largest_block) largest_block = block.size;
  }
  return largest_block;
}

template<class BaseSpace>
std::size_t DynamicSizePool<BaseSpace>::getReleasableSize() const noexcept
{
  std::size_t n_blocks = 0;
  std::size_t n_bytes  = 0;

  for (const auto& block : freeBlocks) {
    if (block.size == block.blockSize) {
      n_bytes += block.blockSize;
      n_blocks += 1;
    }
  }
  return n_blocks > 1 ? n_bytes : 0;
}

template<class BaseSpace>
std::size_t DynamicSizePool<BaseSpace>::getFreeBlocks() const noexcept
{
  std::size_t nb{0};

  for (const auto& block : freeBlocks) {
    if (block.size == block.blockSize) nb++;
  }
  return nb;
}

template<class BaseSpace>
void DynamicSizePool<BaseSpace>::coalesce(std::size_t suggested_size)
{
  if (getFreeBlocks() > 1) {
    freeReleasedBlocks();
    std::size_t size_post{actual_bytes_};

    if (size_post < suggested_size) {
      std::size_t alloc_size{suggested_size - size_post};
      logger_->debug("Coalescing {:}B (actual: {:}B, requested: {:}B)",
                     alloc_size,
                     size_post,
                     suggested_size);
      if (void* ptr = impl_allocate(alloc_size); ptr != nullptr) {
        impl_deallocate(ptr);
        if (auto new_size = should_coalesce(); new_size != 0) {
          coalesce(new_size);
        }
      } else {
        logger_->warn("Failed to coalesce to the requested size ({:}B)",
                      suggested_size);
      }
    }
  }
}

template<class BaseSpace>
void DynamicSizePool<BaseSpace>::coalesce()
{
  auto actual_size = getActualSize();

  std::lock_guard lock(mutex_);
  coalesce(actual_size);
}

template<class BaseSpace>
typename DynamicSizePool<BaseSpace>::BlockIter
DynamicSizePool<BaseSpace>::findUsableBlock(std::size_t size)
{
  auto        best = freeBlocks.end();
  std::size_t best_size{};
  for (auto it = freeBlocks.begin(); it != freeBlocks.end(); ++it) {
    if (it->size >= size
        && (best == freeBlocks.end() || it->size < best_size)) {
      best      = it;
      best_size = best->size;
      if (best_size == size) break; // exact match, stop searching
    }
  }
  return best;
}

template<class BaseSpace>
typename DynamicSizePool<BaseSpace>::BlockIter
DynamicSizePool<BaseSpace>::allocateBlock(std::size_t size)
{
  if (freeBlocks.empty() && usedBlocks.empty()) {
    if (min_pool_size_initial_ > size) size = min_pool_size_initial_;
  } else {
    if (min_pool_size_growth_ > size) size = min_pool_size_growth_;
  }

  logger_->trace("Allocating a new chunk of size {:}", size);

  using ptr_t   = decltype(Block::data);
  auto new_data = static_cast<ptr_t>(allocator_.allocate(size));
  if (new_data == nullptr) {
    logger_->debug("Allocation failed, release free chunks and retry...");
    freeReleasedBlocks();
    new_data = static_cast<ptr_t>(allocator_.allocate(size));
    if (new_data == nullptr) {
      logger_->warn("Allocation failed.");
      return freeBlocks.end();
    }
  }

  actual_bytes_ += size;
  actual_highwatermark_ = (actual_bytes_ > actual_highwatermark_)
                          ? actual_bytes_
                          : actual_highwatermark_;
  releasable_blocks_++;
  total_blocks_++;

  // Insert into free blocks
  // Search from back
  auto pos = freeBlocks.cend();
  while (pos != freeBlocks.cbegin() && new_data < (pos - 1)->data) {
    --pos;
  }

  logger_->trace(
    "Allocated a new chunk of size {:}B @ {}", size, fmt::ptr(new_data));

  return freeBlocks.insert(pos, {new_data, size, size});
}

template<class BaseSpace>
typename DynamicSizePool<BaseSpace>::Block
DynamicSizePool<BaseSpace>::splitBlock(BlockIter const block, std::size_t size)
{
  if (block->size == block->blockSize) releasable_blocks_--;

  Block sub_block = *block;
  if (block->size == size) {
    // Remove the entire block
    freeBlocks.erase(block);
  } else {
    // Split the block
    block->data += size;
    block->size -= size;
    block->blockSize = 0;
    sub_block.size   = size;
  }
  return sub_block;
}

template<class BaseSpace>
void DynamicSizePool<BaseSpace>::releaseBlock(Block block)
{
  // Search the location to insert in the freeBlocks list
  auto pos = freeBlocks.end();
  while (pos != freeBlocks.begin() && block.data < (pos - 1)->data) {
    --pos;
  }

  // Merge the next block (pointed by pos) if possible
  if (pos != freeBlocks.end() && pos->blockSize == 0 // must be a sub-block
      && block.data + block.size == pos->data) {
    block.size += pos->size;
    pos = freeBlocks.erase(pos);
  }

  // Check if the previous block can be merged
  if (pos != freeBlocks.begin()) {
    auto prev = pos - 1;
    if (block.blockSize == 0 // the block to be inserted must be sub-block
        && prev->data + prev->size == block.data) {
      // assign the previous block to merge the current block
      prev->size += block.size;
      if (prev->size == prev->blockSize) ++releasable_blocks_;
      return;
    }
  }

  freeBlocks.insert(pos, block);
  if (block.size == block.blockSize) ++releasable_blocks_;
}

template<class BaseSpace>
std::size_t DynamicSizePool<BaseSpace>::should_coalesce() const noexcept
{
  return std::visit(
    [this](auto&& heuristic) -> std::size_t {
      using T = std::decay_t<decltype(heuristic)>;
      if constexpr (std::is_same_v<T, percent_releasable>) {
        auto percentage = heuristic.percentage;
        if (percentage == 0) return 0u;
        if (percentage == 100) {
          return current_size_ == 0 ? actual_bytes_ : 0u;
        }
        auto f = (float)((float)percentage / 100.0f);
        // Calculate threshold in bytes from the percentage
        auto threshold = static_cast<std::size_t>(f * actual_bytes_);
        return getReleasableSize() >= threshold ? actual_bytes_ : 0;
      } else if constexpr (std::is_same_v<T, blocks_releasable>) {
        auto n_blocks = heuristic.n_blocks;
        return releasable_blocks_ > n_blocks ? actual_bytes_ : 0;
      }
      return 0;
    },
    coalesce_heuristic_);
}

template<class BaseSpace>
std::size_t DynamicSizePool<BaseSpace>::freeReleasedBlocks()
{
  std::size_t freed{0};

  auto it = freeBlocks.cbegin();
  while (it != freeBlocks.cend()) {
    if (it->size == it->blockSize) {
      logger_->trace(
        "Releasing chunk @ {} ({:}B)", fmt::ptr(it->data), it->size);

      allocator_.deallocate(it->data);

      actual_bytes_ -= it->size;
      releasable_blocks_--;
      total_blocks_--;

      freed += it->size;

      it = freeBlocks.erase(it);
    } else {
      ++it;
    }
  }
  return freed;
}

#define SPECIALIZE_POOL_CONSTRUCTOR(MEM_SPACE)          \
  template<>                                            \
  DynamicSizePool<Kokkos::MEM_SPACE>::DynamicSizePool(  \
    std::size_t           alignment,                    \
    std::size_t           minimum_pool_size_initial,    \
    std::size_t           minimum_pool_size_growth,     \
    PoolCoalesceHeuristic coalesce_heuristic)           \
    : allocator_(alignment)                             \
    , coalesce_heuristic_{coalesce_heuristic}           \
    , min_pool_size_initial_{minimum_pool_size_initial} \
    , min_pool_size_growth_{minimum_pool_size_growth}   \
    , logger_{alps::get_logger(#MEM_SPACE "_pool"s)}    \
  {                                                     \
    freeBlocks.reserve(1u << 6);                        \
    usedBlocks.reserve(1u << 6);                        \
    logger_->trace("Memory pool created");              \
  }

SPECIALIZE_POOL_CONSTRUCTOR(HostSpace)
template class DynamicSizePool<Kokkos::HostSpace>;
#if defined(KOKKOS_ENABLE_CUDA)
SPECIALIZE_POOL_CONSTRUCTOR(CudaSpace)
template class DynamicSizePool<Kokkos::CudaSpace>;
SPECIALIZE_POOL_CONSTRUCTOR(CudaUVMSpace)
template class DynamicSizePool<Kokkos::CudaUVMSpace>;
SPECIALIZE_POOL_CONSTRUCTOR(CudaHostPinnedSpace)
template class DynamicSizePool<Kokkos::CudaHostPinnedSpace>;
#elif defined(KOKKOS_ENABLE_HIP)
SPECIALIZE_POOL_CONSTRUCTOR(HIPSpace)
template class DynamicSizePool<Kokkos::HIPSpace>;
SPECIALIZE_POOL_CONSTRUCTOR(HIPManagedSpace)
template class DynamicSizePool<Kokkos::HIPManagedSpace>;
SPECIALIZE_POOL_CONSTRUCTOR(HIPHostPinnedSpace)
template class DynamicSizePool<Kokkos::HIPHostPinnedSpace>;
#endif

} // namespace alps::memory
