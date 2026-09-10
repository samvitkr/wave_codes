//
// Created by xuanx004 on 7/16/24.
//

#include "pressure_grad_options.h"

#include <common/real_type.h>

#include <fmt/format.h>

#include <limits>

namespace alps::solver {

PressureGradOptions PressureGradOptions::parse_from(ConfigTable const& config)
{
  if (!config.contains(config_key)) {
    return {};
  }
  auto const table = config.extract_table(config_key);

  if (auto const pg_type = table.get_value_or<std::string>("type", "constant");
      pg_type != "constant") {
    throw std::runtime_error("Only constant pressure gradient is supported");
  }

  PressureGradOptions options;
  options.gx = table.get_value_or<double>("values[0]", 0);
  options.gy = table.get_value_or<double>("values[1]", 0);

  auto const near_zero = [](double x) {
    auto const abs_x = std::abs(x);
    return abs_x < (double)std::numeric_limits<Real>::epsilon() && abs_x > 0;
  };
  if (near_zero(options.gx) && near_zero(options.gy)) {
    throw std::runtime_error(fmt::format(
      "Pressure gradient is too small: ({}, {})", options.gx, options.gy));
  }

  if (options.gx != 0 || options.gy != 0) {
    options.enabled = true;
  }

  return options;
}

} // namespace alps::solver
