#include "logging.h"

#include <mpipp/comm.h>
#include <mpipp/environment.h>
#include <spdlog/cfg/env.h>
#include <spdlog/sinks/rotating_file_sink.h>

#include <mutex>
#include <string>

namespace alps {

namespace {
Logger initialize_default_logger(std::string logger_name)
{
  spdlog::cfg::load_env_levels();

  if (mpipp::initialized()) {
    spdlog::default_logger()->set_pattern(get_logger_pattern_mpi());
  } else {
    spdlog::default_logger()->set_pattern(get_logger_pattern_nompi());
  }
  spdlog::set_default_logger(spdlog::default_logger()->clone(logger_name));

  return spdlog::default_logger();
}
} // anonymous namespace

[[nodiscard]] std::string get_logger_pattern_mpi()
{
  if (!mpipp::initialized()) {
    throw std::runtime_error(
      "MPI logger can only be setup after MPI has been initialized");
  }
  auto const rank_str = std::to_string(mpipp::COMM_WORLD().rank());
  // format: [time][level][rank][logger_name] message
  return "[%T.%e][%L][rk" + rank_str + "][%n] %v";
}

[[nodiscard]] std::string get_logger_pattern_nompi()
{
  // format: [time][level][process][logger_name] message
  return "[%T.%e][%L][%P][%n] %v";
}

[[nodiscard]] Logger default_logger()
{
  static std::once_flag _default_logger_flag;
  std::call_once(_default_logger_flag,
                 []() { initialize_default_logger("main"); });
  return spdlog::default_logger();
}

[[nodiscard]] Logger get_logger(std::string logger_name)
{
  return default_logger()->clone(std::move(logger_name));
}

RotatingFileSinkConfig::RotatingFileSinkConfig(std::string filename_)
  : filename(std::move(filename_))
{
  if (this->filename.empty()) {
    throw std::invalid_argument("Filename cannot be empty");
  }
}

RotatingFileSinkConfig::RotatingFileSinkConfig(std::string filename_,
                                               std::size_t max_size_,
                                               std::size_t max_files_)
  : filename(std::move(filename_))
  , max_size(max_size_)
  , max_files(max_files_)
{
  if (this->filename.empty()) {
    throw std::invalid_argument("Filename cannot be empty");
  }
}

RotatingFileSinkConfig&
RotatingFileSinkConfig::set_max_size(std::size_t max_size_)
{
  this->max_size = max_size_;
  return *this;
}

RotatingFileSinkConfig&
RotatingFileSinkConfig::set_max_files(std::size_t max_files_)
{
  this->max_files = max_files_;
  return *this;
}

[[nodiscard]] Logger create_logger(std::string            logger_name,
                                   RotatingFileSinkConfig config)
{
  auto new_logger = spdlog::rotating_logger_mt(
    logger_name, config.filename, config.max_size, config.max_files);
  return new_logger;
}

} // namespace alps
