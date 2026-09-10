//
// Created by xuananqing on 4/1/23.
//

#pragma once

#include <common/base/logging_fwd.h>
#include <solvers/mesh/mesh_fwd.h>

#include <filesystem>
#include <memory>
#include <string>
#include <vector>

namespace alps::diagnostics {

class FieldStatistics;

/// A class for managing different statistics
/** To use this class, you need to provide a mesh to construct the manager.
 * Then, attach statistics objects to the manager. When the statistics are
 * needed, call calculate() and corresponding write functions. When statistics
 * are no longer needed, call cleanup() to release the storage of the cached
 * statistics.
 */
class StatisticalDiagnosticsManager
{
 public:
  explicit StatisticalDiagnosticsManager(solver::Mesh const& mesh);

  explicit StatisticalDiagnosticsManager(solver::BottomWaveMesh const& mesh);

  explicit StatisticalDiagnosticsManager(solver::TopWaveMesh const& mesh);

  StatisticalDiagnosticsManager&
  attach_statistics(std::unique_ptr<FieldStatistics>&& stat);

  StatisticalDiagnosticsManager&
  attach_statistics(std::vector<std::unique_ptr<FieldStatistics>> stat);

  void calculate();

  void write_to_tecplot(
    std::string           zone_name,
    std::filesystem::path directory = std::filesystem::current_path()) const;

  void cleanup();

  ~StatisticalDiagnosticsManager();

  std::vector<std::unique_ptr<FieldStatistics>> field_stats;

 private:
  Logger logger;
};

} // namespace alps::diagnostics
