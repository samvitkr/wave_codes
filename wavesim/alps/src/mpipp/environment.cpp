#include "environment.h"

#include "config.h"

#include <enum.hpp/enum.hpp>
#include <spdlog/spdlog.h>

#include <stdexcept>

namespace mpipp {

ENUM_HPP_TRAITS_DECL(threading_modes, (single)(funneled)(serialized)(multiple))

void finalize() noexcept
{
  logger()->debug("Finalize MPI.");
  if (!::mpipp::finalized()) {
    if (auto result = MPI_Finalize(); result != MPI_SUCCESS) {
      logger()->error("MPI finalization failed with error code {}.", result);
    }
  }
}

bool initialized() noexcept
{
  int flag{};
  MPI_Initialized(&flag);
  return flag != 0;
}

bool finalized() noexcept
{
  int flag{};
  MPI_Finalized(&flag);
  return flag != 0;
}

threading_modes init(int* argc, char*** argv, threading_modes thread_mode)
{
  if (finalized()) {
    throw std::runtime_error("MPI has already been finalized");
  }
  int provided{};
  if (!initialized()) {
    logger()->debug("Initializing MPI with threading mode {}.",
                    threading_modes_traits::to_string_or_empty(thread_mode));
    if (MPI_SUCCESS
        != MPI_Init_thread(argc, argv, int(thread_mode), &provided)) {
      throw std::runtime_error("MPI initialization failed.");
    }
  } else {
    provided = threading_modes_traits::to_underlying(threading_mode());
    int rank{};
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
      logger()->warn("MPI already initialized.");
    }
  }

  if (provided < threading_modes_traits::to_underlying(thread_mode)) {
    int rank{};
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
      logger()->warn(
        "MPI threading level support ({}) is less than requested {}.",
        threading_modes_traits::to_string_or_empty(threading_modes(provided)),
        threading_modes_traits::to_string_or_empty(thread_mode));
    }
  }
  return static_cast<threading_modes>(provided);
}

threading_modes init(threading_modes thread_mode)
{
  return init(nullptr, nullptr, thread_mode);
}

int tag_up()
{
  void* p;
  int   flag;
  MPI_Comm_get_attr(MPI_COMM_WORLD, MPI_TAG_UB, &p, &flag);
  return *static_cast<int*>(p);
}

std::string processor_name()
{
  char name[MPI_MAX_PROCESSOR_NAME];
  int  len{};
  MPI_Get_processor_name(static_cast<char*>(name), &len);
  return name;
}

threading_modes threading_mode()
{
  int thread_mode_{};
  MPI_Query_thread(&thread_mode_);
  return static_cast<threading_modes>(thread_mode_);
}

bool is_thread_main()
{
  int res{};
  MPI_Is_thread_main(&res);
  return static_cast<bool>(res);
}

bool wtime_is_global()
{
  void* p{};
  int   flag{};
  MPI_Comm_get_attr(MPI_COMM_WORLD, MPI_WTIME_IS_GLOBAL, &p, &flag);
  return *static_cast<int*>(p) != 0;
}

std::shared_ptr<spdlog::logger> logger()
{
  static auto _logger = spdlog::default_logger()->clone("mpi");
  return _logger;
}

} // namespace mpipp
