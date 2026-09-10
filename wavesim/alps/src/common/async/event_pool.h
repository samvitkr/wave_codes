#pragma once

#include "event_cuda.h"
#include "event_openmp.h"
#include <common/base/logging_fwd.h>

#include <mutex>
#include <vector>

namespace alps {
namespace async {
template<class ExecutionSpace>
class event_pool
{
 public:
  using event_t  = async::queue_event<ExecutionSpace>;
  using handle_t = typename event_t::handle_t; // underlying event type

 private:
  static constexpr std::size_t default_incremental_size = 6u;

 public:
  event_pool() noexcept;

  void resize(std::size_t new_size);

  void release(handle_t e) noexcept;

  [[nodiscard]] event_t acquire() noexcept;

  void clear() noexcept;

  ~event_pool();

 private:
  void impl_add(std::size_t n);

  std::vector<handle_t> events_;
  std::size_t           total_size_{0};
  alps::Logger          logger_;
  mutable std::mutex    mutex_;
};

} // namespace async
} // namespace alps
