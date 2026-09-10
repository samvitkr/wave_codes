#pragma once

#include <common/base/logging_fwd.h>

#include <Kokkos_Core_fwd.hpp>

#include <atomic>
#include <vector>

#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE
#if defined(KOKKOS_ENABLE_OPENMP)
#include <OpenMP/Kokkos_OpenMP.hpp> // needed for definition of Kokkos::OpenMP
#endif
#if defined(KOKKOS_ENABLE_CUDA)
#include <Cuda/Kokkos_Cuda.hpp> // needed for definition of Kokkos::Cuda
#endif
#if defined(KOKKOS_ENABLE_HIP)
#include <HIP/Kokkos_HIP.hpp> // needed for definition of Kokkos::HIP
#endif
#undef KOKKOS_IMPL_PUBLIC_INCLUDE
#endif

namespace alps {

/// @brief A pool of instances of Kokkos execution spaces.
/**
 * Provids access to a pool of asynchronous execution streams (or queues) for
 * supported backends, such as CUDA.
 *
 * The default implementation always return the same instance, i.e. multiple
 * execution streams are not supported.
 */
template<class ExecutionSpace>
class StreamPool
{
 public:
  using non_owning_space_t = ExecutionSpace;

  explicit StreamPool(int pool_size) noexcept
    : size{pool_size}
  {}

  StreamPool(StreamPool&&)      = delete;
  StreamPool(StreamPool const&) = delete;

  StreamPool& operator=(StreamPool&&)      = delete;
  StreamPool& operator=(StreamPool const&) = delete;

  non_owning_space_t get_next_stream() const noexcept { return {}; }

  non_owning_space_t get_stream(int /*stream_id*/) const { return {}; }

  static non_owning_space_t get_default_stream() noexcept { return {}; }

  non_owning_space_t operator[](int /*stream_id*/) const { return {}; }

  void fence() const {}

  int actual_size() const { return 0; }

  void clear() const {}

  int size;
};

#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)

/// @brief A pool of instances of Kokkos Cuda spaces.
/**
 * @note This class maintains a circular buffer of CUDA streams, i.e.
 * get_next_stream() returns streams in a round-robin fashion. Calls to
 * get_next_stream() are thread safe.
 */
template<>
class StreamPool<Kokkos::DefaultExecutionSpace>
{
 public:
  using non_owning_space_t = Kokkos::DefaultExecutionSpace;
  using space_t            = Kokkos::DefaultExecutionSpace;

  explicit StreamPool(int pool_size) noexcept;

  StreamPool(StreamPool&&)      = delete;
  StreamPool(StreamPool const&) = delete;

  StreamPool& operator=(StreamPool&&)      = delete;
  StreamPool& operator=(StreamPool const&) = delete;

  /// Return the next stream in the pool in loop. Thread safe.
  non_owning_space_t get_next_stream() const noexcept
  {
    if (streams_.empty()) initialize();
    next_stream_ = (next_stream_ + 1) % streams_.size();
    return streams_[next_stream_];
  }

  /// Get a stream in the pool by id. The stream 0 is the default stream.
  non_owning_space_t get_stream(int stream_id) const
  {
    if (stream_id == 0) return {};
    if (streams_.empty()) initialize();
    return streams_.at(stream_id - 1);
  }

  /// Get the default stream (stream 0).
  static non_owning_space_t get_default_stream() noexcept { return {}; }

  /// Get a stream in the pool by id. The stream 0 is the default stream.
  non_owning_space_t operator[](int stream_id) const
  {
    return get_stream(stream_id);
  }

  void fence() const
  {
    for (auto const& stream : streams_) {
      stream.fence();
    }
  }

  int actual_size() const { return static_cast<int>(streams_.size()); }

  void clear() const;

  ~StreamPool();

  int size;

 private:
  void initialize() const;

  mutable std::vector<space_t> streams_{};
  mutable std::atomic_size_t   next_stream_{};

  Logger logger_;
};
#endif

/// Obtain a stream from a stream pool
template<class ExecSpace>
auto get_next_stream(const StreamPool<ExecSpace>& stream_pool)
{
  return stream_pool.get_next_stream();
}

template<class ExecSpace>
void fence(std::vector<ExecSpace>&& streams)
{
  for (auto&& stream : streams) {
    stream.fence();
  }
}

} // namespace alps
