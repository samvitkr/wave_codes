//
// Created by xuanx004 on 7/16/24.
//

#include "rayleigh_damp_options.h"

#include <common/real_type.h>

#include <fmt/format.h>

#include <limits>

namespace alps::solver {

RayleighDampOptionsCommon
RayleighDampOptionsCommon::parse_from(ConfigTable const& config)
{
  RayleighDampOptionsCommon options;
  options.factor = config.get_value_or<double>("factor", 1.0);
  options.z_zero = config.get_value<double>("range[0]");
  options.z_one  = config.get_value<double>("range[1]");

  auto const near_zero = [](double x) {
    auto const abs_x = std::abs(x);
    return abs_x < (double)std::numeric_limits<Real>::epsilon() && abs_x > 0;
  };
  auto const valid_range = [](double z0, double z1) {
    auto z_min = std::min(z0, z1);
    auto z_max = std::max(z0, z1);
    return z_min <= 1 && z_max >= 0;
  };
  if (near_zero(options.factor)) {
    throw std::runtime_error(
      fmt::format("Rayleigh damping factor is too small: {}", options.factor));
  }
  if (!valid_range(options.z_zero, options.z_one)) {
    throw std::runtime_error(
      fmt::format("Rayleigh damping range ({}, {}) does not overlap 0..1",
                  options.z_zero,
                  options.z_one));
  }

  if (options.factor != 0) {
    options.enabled = true;
  }

  return options;
}

RayleighDampOptions RayleighDampOptions::parse_from(ConfigTable const& config)
{
  if (!config.contains(config_key)) {
    return RayleighDampOptions();
  }
  auto const table = config.extract_table(config_key);

  RayleighDampOptions options = RayleighDampOptionsCommon::parse_from(table);
  auto                b       = table.get_value_or<double>("b", 0.0);
  auto                a       = table.get_value_or<double>("a", 0.0);
  options.u_b                 = table.get_value_or<double>("u_b", b);
  options.u_a                 = table.get_value_or<double>("u_a", a);
  options.v_b                 = table.get_value_or<double>("v_b", b);
  options.v_a                 = table.get_value_or<double>("v_a", a);

  return options;
}

ScalarRayleighDampOptions
ScalarRayleighDampOptions::parse_from(ConfigTable const& config)
{
  if (!config.contains(config_key)) {
    return ScalarRayleighDampOptions();
  }
  auto const table = config.extract_table(config_key);

  ScalarRayleighDampOptions options =
    RayleighDampOptionsCommon::parse_from(table);
  options.b = table.get_value_or<double>("b", 0.0);
  options.a = table.get_value_or<double>("a", 0.0);

  return options;
}
} // namespace alps::solver
