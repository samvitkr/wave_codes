//
// Created by xuananqing on 7/4/23.
//

#include "options.h"

namespace alps::solver {
FreeSurfaceSolverOptions
FreeSurfaceSolverOptions::parse_from(ConfigTable const& config)
{
  FreeSurfaceSolverOptions options{base_t::parse_from(config)};

  options.Fr2 = config.get_value<double>("Fr2");
  options.WeR = config.get_value<double>("WeR");

  options.PressureEqnIterations = config.get_value_or("Iter_MaxN", 7);
  if (options.PressureEqnIterations <= 0) {
    throw std::invalid_argument("Iter_MaxN must be positive.");
  }

  options.PressureEqnAbsTol = config.get_value_or("Iter_eps", 1e-6);
  if (options.PressureEqnAbsTol <= 0) {
    throw std::invalid_argument("Iter_eps must be positive.");
  }

  options.PressureEqnRelTol = config.get_value_or("Iter_rel_eps", 1e-12);
  if (options.PressureEqnRelTol <= 0) {
    throw std::invalid_argument("Iter_rel_eps must be positive.");
  }

  return options;
}
} // namespace alps::solver
