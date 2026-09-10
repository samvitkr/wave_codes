#pragma once

#include <spdlog/fwd.h>

#include <memory>
#include <string>

namespace alps {

using Logger = std::shared_ptr<spdlog::logger>;

/// Returns the string pattern (incl. rank) for the MPI logger, only valid after
/// MPI has been initialized
[[nodiscard]] std::string get_logger_pattern_mpi();

/// Returns the string pattern for the non-MPI logger
[[nodiscard]] std::string get_logger_pattern_nompi();

[[nodiscard]] Logger default_logger();

/// Clone a new MPI logger with the given name
[[nodiscard]] Logger get_logger(std::string logger_name);

struct RotatingFileSinkConfig
{
  std::string filename;
  std::size_t max_size{(std::size_t)1024 * 1024 * 10}; // 10MB
  std::size_t max_files{2};

  explicit RotatingFileSinkConfig(std::string filename_);

  RotatingFileSinkConfig(std::string filename_,
                         std::size_t max_size_,
                         std::size_t max_files_);

  /// Set the maximum size of the log file in bytes
  RotatingFileSinkConfig& set_max_size(std::size_t max_size_);

  RotatingFileSinkConfig& set_max_files(std::size_t max_files_);
};

/// Create a new rotating file logger with the given name and configuration
[[nodiscard]] Logger create_logger(std::string            logger_name,
                                   RotatingFileSinkConfig config);
} // namespace alps
