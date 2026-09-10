#include "field_statistics.h"

namespace alps::diagnostics {

FieldStatistics::FieldStatistics(
  std::unordered_set<FieldStatisticsType> statistic_types_)
  : statistic_types{std::move(statistic_types_)}
{}

FieldStatistics::~FieldStatistics() = default;

void FieldStatistics::cleanup() {}

std::vector<std::string>
FieldStatistics::get_variable_labels(FieldStatisticsType /*type*/) const
{
  return {};
}

std::vector<int> FieldStatistics::get_shape(FieldStatisticsType /*type*/) const
{
  return {};
}

void FieldStatistics::write_to_tecplot(std::ostream& /*out*/,
                                       FieldStatisticsType /*type*/) const
{}

} // namespace alps::diagnostics
