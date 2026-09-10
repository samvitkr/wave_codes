#pragma once

#include <common/program_options/config_table.h>
#include <solvers/ns/channel_options.h>

#include <string>

namespace alps::solver {
struct FlowOverWaveSolverOptions : ChannelSolverOptions
{
  std::string integrator{"ab2"};

  double PressureEqnAbsTol{1e-6};
  double PressureEqnRelTol{1e-12};
  int    PressureEqnIterations{7};

  std::string DiffusionCNZetaAlgorithm{"tdma"};
  int         DiffusionCNZetaJacobiMaxIters{100};
  double      DiffusionCNZetaJacobiAbsTol{1e-8};
  double      DiffusionCNZetaJacobiRelTol{1e-12};

  static FlowOverWaveSolverOptions parse_from(ConfigTable const& config);

 private:
  using base_t = ChannelSolverOptions;
};
} // namespace alps::solver
