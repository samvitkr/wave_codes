#pragma once

#include "config.h"

#include <spdlog/fwd.h>

#include <memory>
#include <string>

namespace mpipp {

enum class threading_modes
{
  single     = MPI_THREAD_SINGLE,
  funneled   = MPI_THREAD_FUNNELED,
  serialized = MPI_THREAD_SERIALIZED,
  multiple   = MPI_THREAD_MULTIPLE
};

bool initialized() noexcept;

bool finalized() noexcept;

threading_modes init(threading_modes thread_mode = threading_modes::single);

void finalize() noexcept;

int tag_up();

std::string processor_name();

threading_modes threading_mode();

bool is_thread_main();

bool wtime_is_global();

std::shared_ptr<spdlog::logger> logger();

} // namespace mpipp
