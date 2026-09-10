//
// Created by xuananqing on 4/1/23.
//

#pragma once

#include <enum.hpp/enum.hpp>

#include <memory>
#include <ostream>
#include <string>
#include <unordered_set>
#include <vector>

namespace alps::solver {
class FlowField;
} // namespace alps::solver

namespace alps::diagnostics {

ENUM_HPP_CLASS_DECL(FieldStatisticsType, int, (XY = 0)(Y))

class FieldStatistics;

/// @brief Create a velocity statistics object from a given flow field
std::unique_ptr<FieldStatistics> create_velocity_statistics(
  alps::solver::FlowField const&          flow_field,
  std::unordered_set<FieldStatisticsType> field_statistics_types);

/// @brief Create statistics objects for each scalar in a given flow field
std::vector<std::unique_ptr<FieldStatistics>> create_scalar_statistics(
  alps::solver::FlowField const&          flow_field,
  std::unordered_set<FieldStatisticsType> field_statistics_types);

/// @brief Base class for computing and writing statistics of a field variable
class FieldStatistics
{
 public:
  /// @brief The types of the field statistics to be calculated
  std::unordered_set<FieldStatisticsType> statistic_types;

  /// @brief Calculate the statistics
  virtual void calculate() = 0;

  /// @brief Cleanup (release the stored statistics)
  virtual void cleanup();

  /// @brief Query the shape of the statistics
  /**
   * @note On processes that do not write statistics, return an empty vector.
   */
  virtual std::vector<int> get_shape(FieldStatisticsType type) const;

  /// @brief Query the variable labels
  virtual std::vector<std::string>
  get_variable_labels(FieldStatisticsType type) const;

  virtual void write_to_tecplot(std::ostream&       out,
                                FieldStatisticsType type) const;

  virtual ~FieldStatistics();

 protected:
  explicit FieldStatistics(
    std::unordered_set<FieldStatisticsType> statistic_types_);
};

} // namespace alps::diagnostics
