#include "timeout_mpi_barrier.h"

#include <common/base/logging.h>

#include <chrono>
#include <future>
#include <thread>

namespace alps::utils {

bool timeout_mpi_barrier(const mpipp::communicator& comm, std::size_t timeout)
{
  default_logger()->debug("Timeout MPI barrier with {} s", timeout);

  // post a non-blocking barrier
  MPI_Request req{MPI_REQUEST_NULL};
  MPI_Ibarrier(comm.raw_handle(), &req);
  int flag{0};
  MPI_Test(&req, &flag, MPI_STATUS_IGNORE);
  if (flag) {
    return true;
  }

  // start a timer in a separate thread for the timeout
  for (size_t i = 0; i < timeout; ++i) {
    auto fut = std::async(std::launch::async, [] {
      std::this_thread::sleep_for(std::chrono::seconds(1));
    });
    fut.wait();

    MPI_Test(&req, &flag, MPI_STATUS_IGNORE);
    if (flag) {
      return true;
    }
  }
  default_logger()->warn("MPI barrier timeout");
  // do not free or cancel the request.
  return false;
}
} // namespace alps::utils
