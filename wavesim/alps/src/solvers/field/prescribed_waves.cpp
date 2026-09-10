//
// Created by xuanx004 on 1/21/24.
//

#include "prescribed_waves.h"

#include <common/program_options/config_table.h>

namespace alps::solver {

MonochromaticWave MonochromaticWave::parse_from(ConfigTable const& config)
{
  MonochromaticWave wave;

  wave.n_waves = config.get_value_or("wave_n", 0);
  wave.ak      = config.get_value_or("wave_ak", 0.0);
  wave.phase_c = config.get_value_or("wave_c", 0.0);

  return wave;
}

} // namespace alps::solver
