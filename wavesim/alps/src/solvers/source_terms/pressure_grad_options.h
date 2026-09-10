//
// Created by xuanx004 on 7/16/24.
//

#pragma once

#include <common/program_options/config_table.h>

#include <string_view>

namespace alps::solver {

struct PressureGradOptions
{
  bool   enabled{false};
  double gx;
  double gy;

  static constexpr std::string_view config_key{"BodyForce.PressureGradient"};

  [[nodiscard]] static PressureGradOptions
  parse_from(ConfigTable const& config);
};

} // namespace alps::solver
