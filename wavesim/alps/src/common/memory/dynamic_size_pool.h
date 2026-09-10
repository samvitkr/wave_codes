//////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2016-22, Lawrence Livermore National Security, LLC and Umpire
// project contributors. See the COPYRIGHT file for details.
//
// SPDX-License-Identifier: (MIT)
//////////////////////////////////////////////////////////////////////////////
#pragma once

#include "aligned_allocation.h"

#include <common/base/logging_fwd.h>

#include <cstddef>
#include <mutex>
#include <variant>
#include <vector>

namespace alps::memory {
template<class BaseSpace>
class base_allocator;

struct percent_releasable
{
  int percentage;
};

struct blocks_releasable
{
  std::size_t n_blocks;
};

using PoolCoalesceHeuristic =
  std::variant<percent_releasable, blocks_releasable>;

template<class BaseSpace>
class DynamicSizePool
{
 public:
  static constexpr std::size_t DEFAULT_BLOCK_SIZE_INITIAL{
    static_cast<std::size_t>(4) * 1024 * 1024};
  static constexpr std::size_t DEFAULT_BLOCK_SIZE_GROWTH{
    static_cast<std::size_t>(1) * 1024 * 1024};

 private:
  struct Block
  {
    std::byte*  data;
    std::size_t size;
    std::size_t blockSize;
  };

  using BlockIter = typename std::vector<Block>::iterator;

 public:
  static DynamicSizePool& global_instance();

  DynamicSizePool(
    std::size_t alignment,
    std::size_t minimum_pool_size_initial    = DEFAULT_BLOCK_SIZE_INITIAL,
    std::size_t minimum_pool_size_growth     = DEFAULT_BLOCK_SIZE_GROWTH,
    PoolCoalesceHeuristic coalesce_heuristic = percent_releasable{100});

  DynamicSizePool(const DynamicSizePool&) = delete;

  ~DynamicSizePool();

  [[nodiscard]] void* allocate(std::size_t bytes);

  void deallocate(void* ptr);

  void release();

  std::size_t getActualSize() const noexcept;

  std::size_t getCurrentSize() const noexcept;

  std::size_t getActualHighwaterMark() const noexcept;

  std::size_t getLargestAvailableBlock() const noexcept;

  void coalesce();

 private:
  [[nodiscard]] void* impl_allocate(std::size_t bytes);

  void impl_deallocate(void* ptr);

  /// Search the list of free blocks and return an iterator to the best block
  BlockIter findUsableBlock(std::size_t size);

  /// Allocate a new block and add it to the list of free blocks
  [[nodiscard]] BlockIter allocateBlock(std::size_t size);

  /// Obtain a sub-block with the requested size from the specified block
  [[nodiscard]] Block splitBlock(BlockIter block, std::size_t size);

  std::size_t getReleasableSize() const noexcept;

  std::size_t getFreeBlocks() const noexcept;

  /// Release the specified block into freeBlocks
  void releaseBlock(Block block);

  void coalesce(std::size_t suggested_size);

  std::size_t should_coalesce() const noexcept;

  // Release all unsplitted blocks from freeBlocks
  std::size_t freeReleasedBlocks();

  mutable std::mutex mutex_;

  // Allocator for the underlying data
  AlignedAllocation<base_allocator<BaseSpace>> allocator_;

  // Start of the nodes of used and free block lists
  std::vector<Block> usedBlocks{};
  std::vector<Block> freeBlocks{};

  PoolCoalesceHeuristic coalesce_heuristic_;

  std::byte zero_byte_data{0};

  // Total size allocated (bytes)
  std::size_t actual_bytes_{0};
  std::size_t current_size_{0};
  std::size_t actual_highwatermark_{0};

  // Minimum size of initial block
  std::size_t min_pool_size_initial_;

  // Minimum size for growth block
  std::size_t min_pool_size_growth_;

  std::size_t releasable_blocks_{0};
  std::size_t total_blocks_{0};

  alps::Logger logger_;
};
} // namespace alps::memory
