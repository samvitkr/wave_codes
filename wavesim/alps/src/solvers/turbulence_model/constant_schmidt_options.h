//
// Created by xuanx004 on 7/14/24.
//

#pragma once

#include <common/program_options/config_table.h>

#include <array>
#include <string_view>

namespace alps::solver {

struct ConstantTurbulentSc;

/// @brief Options for the constant turbulent Schmidt/Prandtl number model
/**
 * The parameters for this model may be specified using a table in the config
 * file as:
 * ```toml
 * model = "constantsct" # or "constantprt"
 * Sct = 1.0 # or Prt = 1.0
 * ```
 *
 * The "model" is handled by the solver, and the other parameters are parsed
 * by the @ref parse_from function.
 */
struct ConstantTurbulentScOptions
{
  double Sct{1.0}; // constant turbulent Schmidt/Prandtl number

  using ModelType = ConstantTurbulentSc;
  static constexpr std::array<std::string_view, 4> match_names{"constantsct",
                                                               "constantprt",
                                                               "constsct",
                                                               "constprt"};

  [[nodiscard]] static ConstantTurbulentScOptions
  parse_from(ConfigTable const& config);
};

} // namespace alps::solver
