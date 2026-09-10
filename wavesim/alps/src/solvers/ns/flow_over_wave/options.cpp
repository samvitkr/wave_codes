#include "options.h"

#include <common/utils/to_lower_case.h>

#include <stdexcept>

namespace alps::solver {
FlowOverWaveSolverOptions
FlowOverWaveSolverOptions::parse_from(ConfigTable const& config)
{
  FlowOverWaveSolverOptions options{base_t::parse_from(config)};

  options.integrator =
    to_lower_case(config.get_value_or<std::string>("integrator", "ab2"));
  if (options.integrator != "ab2" && options.integrator != "ab2cn") {
    throw std::invalid_argument("Unknown integrator type: "
                                + options.integrator);
  }

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

  options.DiffusionCNZetaAlgorithm = to_lower_case(
    config.get_value_or<std::string>("DiffusionCNZetaAlgorithm", "tdma"));
  if (options.DiffusionCNZetaAlgorithm != "tdma"
      && options.DiffusionCNZetaAlgorithm != "pcr"
      && options.DiffusionCNZetaAlgorithm != "jacobi") {
    throw std::invalid_argument("Unknown DiffusionCNZetaAlgorithm: "
                                + options.DiffusionCNZetaAlgorithm);
  }

  options.DiffusionCNZetaJacobiMaxIters =
    config.get_value_or("DiffusionCNZetaJacobiMaxIters", 100);
  if (options.DiffusionCNZetaJacobiMaxIters <= 0) {
    throw std::invalid_argument(
      "DiffusionCNZetaJacobiMaxIters must be positive.");
  }

  options.DiffusionCNZetaJacobiAbsTol =
    config.get_value_or("DiffusionCNZetaJacobiAbsTol", 1e-8);
  if (options.DiffusionCNZetaJacobiAbsTol <= 0) {
    throw std::invalid_argument(
      "DiffusionCNZetaJacobiAbsTol must be positive.");
  }

  options.DiffusionCNZetaJacobiRelTol =
    config.get_value_or("DiffusionCNZetaJacobiRelTol", 1e-12);
  if (options.DiffusionCNZetaJacobiRelTol <= 0) {
    throw std::invalid_argument(
      "DiffusionCNZetaJacobiRelTol must be positive.");
  }

  return options;
}
} // namespace alps::solver
