//
// Created by xuanx004 on 7/16/24.
//

#pragma once

#include <common/program_options/config_table.h>

#include <string_view>

namespace alps::solver {

struct CoriolisOptions
{
  bool   enabled{false};
  double fz;

  static constexpr std::string_view config_key{"BodyForce.Coriolis"};

  [[nodiscard]] static CoriolisOptions parse_from(ConfigTable const& config);
};

} // namespace alps::solver
