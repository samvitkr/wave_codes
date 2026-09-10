//
// Created by xuananqing on 4/1/23.
//

#include "diagnostics.h"

#include "field_statistics.h"
#include "mesh_statistics.h"
#include <common/base/logging.h>
#include <solvers/mesh/curvilinear_mesh.h>
#include <solvers/mesh/mesh.h>

#include <spdlog/fmt/ostr.h>

#include <filesystem>
#include <fstream>
#include <string_view>

namespace alps::diagnostics {

namespace fs = std::filesystem;

namespace {
/// Return the index sizes of the Tecplot header given a shape
std::string get_tecplot_header_index_header(std::vector<int> shape)
{
  if (shape.size() > 3) {
    throw std::invalid_argument("Tecplot supports max 3 dimensions");
  }

  std::string header;
  for (size_t i = 0; i < shape.size(); ++i) {
    header += fmt::format("{}={}", static_cast<char>(73 + (char)i), shape[i]);
    if (i != shape.size() - 1) header += ", ";
  }
  return header;
}
} // namespace

StatisticalDiagnosticsManager::StatisticalDiagnosticsManager(
  alps::solver::Mesh const& mesh)
  : logger{get_logger("stats")}
{
  field_stats.push_back(std::make_unique<MeshStatistics>(mesh));
}

StatisticalDiagnosticsManager::StatisticalDiagnosticsManager(
  alps::solver::BottomWaveMesh const& mesh)
  : logger{get_logger("stats")}
{
  field_stats.push_back(
    std::make_unique<CurvilinearMeshStatistics<solver::BottomWaveMesh>>(mesh));
}

StatisticalDiagnosticsManager::StatisticalDiagnosticsManager(
  alps::solver::TopWaveMesh const& mesh)
  : logger{get_logger("stats")}
{
  field_stats.push_back(
    std::make_unique<CurvilinearMeshStatistics<solver::TopWaveMesh>>(mesh));
}

StatisticalDiagnosticsManager& StatisticalDiagnosticsManager::attach_statistics(
  std::unique_ptr<FieldStatistics>&& stat)
{
  field_stats.push_back(std::move(stat));
  return *this;
}

StatisticalDiagnosticsManager& StatisticalDiagnosticsManager::attach_statistics(
  std::vector<std::unique_ptr<FieldStatistics>> stat)
{
  field_stats.insert(field_stats.end(),
                     std::make_move_iterator(stat.begin()),
                     std::make_move_iterator(stat.end()));
  return *this;
}

void StatisticalDiagnosticsManager::calculate()
{
  for (auto& field_stat : field_stats) {
    field_stat->calculate();
  }
}

namespace {
/* Create a directory if it does not exist, return the actual path to the
 * directory. If the provide path exists but is not a directory, return the
 * current working directory.
 */
fs::path get_directory(fs::path const& directory, Logger logger)
{
  auto actual_directory =
    !directory.empty() ? fs::absolute(directory) : fs::current_path();

  if (fs::exists(actual_directory)) {
    // test if the path points to a directory
    if (fs::is_directory(fs::canonical(actual_directory))) {
      return actual_directory;
    }

    logger->warn("{} already exists and is not a directory. Will write to the "
                 "working directory {}",
                 actual_directory.string(),
                 fs::current_path().string());
    return fs::current_path();
  }

  try {
    fs::create_directories(actual_directory);
  } catch (fs::filesystem_error const& e) {
    logger->warn("Error when creating directory {}: {}. Will write to the "
                 "working directory {}",
                 directory.string(),
                 e.what(),
                 fs::current_path().string());
    actual_directory = fs::current_path();
  }

  return actual_directory;
}

void write_stats_to_tecplot(
  std::vector<std::unique_ptr<FieldStatistics>> const& field_stats,
  FieldStatisticsType                                  stat_type,
  fs::path                                             directory,
  std::string_view                                     zone_name,
  Logger                                               logger)
{
  std::string const filename =
    fmt::format("{}_{}.dat",
                FieldStatisticsType_traits::to_string_or_empty(stat_type),
                zone_name);

  // Get all the variable labels
  std::vector<std::string> variable_labels{};
  for (auto const& field_stat : field_stats) {
    auto const labels = field_stat->get_variable_labels(stat_type);
    variable_labels.insert(
      variable_labels.end(), labels.cbegin(), labels.cend());
  }

  auto const output_shape = field_stats.front()->get_shape(stat_type);

  try {
    // Open the file; on processor that do not write, open /dev/null
    // This avoids opening the same file on different processes while
    // keeping the same stream object
    auto file = !output_shape.empty()
                ? std::ofstream{get_directory(directory, logger) / filename,
                                std::ios::out | std::ios::trunc}
                : std::ofstream{"/dev/null", std::ios::out | std::ios::app};

    // Write the header
    fmt::print(file, "TITLE = {}\n", "Statistical diagnostics");
    fmt::print(file, "VARIABLES =");
    for (auto const& label : variable_labels) {
      fmt::print(file, " \"{}\"", label);
    }
    fmt::print(file, "\n");
    fmt::print(file,
               "ZONE T = \"{}\", {}, F=BLOCK\n",
               zone_name,
               get_tecplot_header_index_header(output_shape));

    for (auto const& field_stat : field_stats) {
      field_stat->write_to_tecplot(file, stat_type);
    }
  } catch (std::exception const& e) {
    logger->error(
      "Error when opening and writing to {}: {}", filename, e.what());
  }
}
} // namespace

void StatisticalDiagnosticsManager::write_to_tecplot(
  std::string           zone_name,
  std::filesystem::path directory) const
{
  decltype(FieldStatistics::statistic_types) statistic_types_available;

  // Collect all types of statistics
  for (auto const& field_stat : field_stats) {
    statistic_types_available.insert(field_stat->statistic_types.cbegin(),
                                     field_stat->statistic_types.cend());
  }

  // For each type of statistics, write to a separate file
  for (auto stat_type : statistic_types_available) {
    write_stats_to_tecplot(
      field_stats, stat_type, directory, zone_name, logger);
  }
}

void StatisticalDiagnosticsManager::cleanup()
{
  for (auto& field_stat : field_stats) {
    field_stat->cleanup();
  }
}

StatisticalDiagnosticsManager::~StatisticalDiagnosticsManager() = default;

} // namespace alps::diagnostics
