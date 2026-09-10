//
// Created by xuananqing on 7/4/23.
//

#pragma once

#include <common/program_options/config_table.h>
#include <solvers/ns/channel_options.h>

#include <limits>

namespace alps::solver {

struct FreeSurfaceSolverOptions : ChannelSolverOptions
{
  double Fr2{std::numeric_limits<double>::infinity()};
  double WeR{0};

  double PressureEqnAbsTol{1e-6};
  double PressureEqnRelTol{1e-12};
  int    PressureEqnIterations{7};

  static FreeSurfaceSolverOptions parse_from(ConfigTable const& config);

 private:
  using base_t = ChannelSolverOptions;
};
} // namespace alps::solver
