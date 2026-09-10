//
// Created by xuanx004 on 1/21/24.
//

#pragma once

#include <common/program_options/config_table.h>

namespace alps::solver {

struct MonochromaticWave
{
  int    n_waves{0};
  double ak{0};
  double phase_c{0};

  MonochromaticWave() = default;

  static MonochromaticWave parse_from(ConfigTable const& config);
};

} // namespace alps::solver
