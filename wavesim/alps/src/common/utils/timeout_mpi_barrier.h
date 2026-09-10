#pragma once

#include <mpipp/comm.h>

#include <cstdint>

namespace alps::utils {

/**
 * @brief Performs a timeout MPI barrier operation.
 *
 * This function blocks until all processes in the given communicator have
 * reached the barrier or the specified timeout has elapsed.
 *
 * @param comm The MPI communicator.
 * @param timeout The timeout value in milliseconds.
 * @return True if all processes reached the barrier within the timeout, false
 * otherwise.
 */
bool timeout_mpi_barrier(mpipp::communicator const& comm, std::size_t timeout);

} // namespace alps::utils
