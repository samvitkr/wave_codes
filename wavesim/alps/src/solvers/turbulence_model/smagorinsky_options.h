//
// Created by xuanx004 on 1/23/24.
//

#pragma once

#include <common/program_options/config_table.h>

#include <array>
#include <string_view>

namespace alps::solver {

struct ConstantSmagorinsky;

/// @brief Options for the constant Smagorinsky model
/**
 * The parameters for this model may be specified using a table in the config
 * file as:
 * ```toml
 * model = "constantsmagorinsky" # or "CS"
 * C0 = 0.14
 * kappa = 0.4
 * wall_z0 = 0.0
 * wall_damp_exp = 2
 * horizontal_filter_width_scale = 1.0
 * debug_C0 = false
 * ```
 *
 * All parameters are optional, and the defaults are listed in the source file.
 *
 * The "model" is handled by the solver, and the other parameters are parsed
 * by the @ref parse_from function.
 */
struct ConstantSmagorinskyOptions
{
  double C0{0.14};         // constant Smagorinsky coefficient
  double kappa{0.4};       // von Karman coefficient
  double wall_z0{0};       // wall roughness
  int    wall_damp_exp{2}; // wall damping exponent

  double horizontal_filter_width_scale{1.0};

  bool enable_debug_C0{false};

  using ModelType = ConstantSmagorinsky;
  // values in the config file to match
  static constexpr std::array<std::string_view, 2> match_names{
    "constantsmagorinsky",
    "cs"};

  [[nodiscard]] static ConstantSmagorinskyOptions
  parse_from(ConfigTable const& config);
};

} // namespace alps::solver
