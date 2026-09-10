#pragma once

#include "manager.h"
#include <common/async/event.h>
#include <common/async/streams.h>

namespace alps {
/// Obtain a stream from the default stream pool
[[nodiscard]] inline auto get_next_stream()
{
  return RuntimeManager::instance().stream_pool().get_next_stream();
}

inline auto& get_default_thread_pool()
{
  return RuntimeManager::instance().thread_pool();
}

[[nodiscard]] inline auto get_device_event()
{
  return RuntimeManager::instance().event_pool().acquire();
}
} // namespace alps
