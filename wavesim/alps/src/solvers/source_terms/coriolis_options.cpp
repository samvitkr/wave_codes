//
// Created by xuanx004 on 7/16/24.
//

#include "coriolis_options.h"

#include <common/real_type.h>

#include <fmt/format.h>

#include <limits>

namespace alps::solver {

CoriolisOptions CoriolisOptions::parse_from(ConfigTable const& config)
{
  if (!config.contains(config_key)) {
    return {};
  }
  auto const table = config.extract_table(config_key);

  CoriolisOptions options;

  options.fz = table.get_value_or<double>("f", 0.0);

  if (Real const abs_fz = std::abs((Real)options.fz);
      abs_fz < std::numeric_limits<Real>::epsilon() && abs_fz > 0) {
    throw std::runtime_error(
      fmt::format("Coriolis parameter is too small: {}", options.fz));
  }

  if (options.fz != 0) {
    options.enabled = true;
  }

  return options;
}

} // namespace alps::solver
