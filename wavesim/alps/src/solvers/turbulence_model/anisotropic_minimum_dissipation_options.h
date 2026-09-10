//
// Created by xuanx004 on 1/23/24.
//

#pragma once

#include <common/program_options/config_table.h>

#include <array>
#include <string_view>

namespace alps::solver {

struct AnisotropicMinimumDissipation;
struct AnisotropicMinimumDissipationScalar;

/// @brief Options for the anisotropic minimum dissipation model
/**
 * The parameters for this model may be specified using a table in the config
 * file as:
 * ```toml
 * model = "anisotropicminimumdissipation" # or "amd"
 * PoincareCxy = 0.083333333
 * PoincareCz = 0.3333
 * horizontal_filter_width_scale = 1.0
 * debug_C0 = false
 * ```
 *
 * All parameters are optional, and the defaults are listed in the source file.
 *
 * The "model" is handled by the solver, and the other parameters are parsed
 * by the @ref parse_from function.
 */
struct AnisotropicMinimumDissipationOptions
{
  double C0{1 / 12.0}; /// Poincaré constant.

  double Cz{1 / 3.0}; /// Poincaré constant for central difference

  double horizontal_filter_width_scale{1.0};

  bool enable_debug_C0{false};

  using ModelType = AnisotropicMinimumDissipation;
  // values in the config file to match
  static constexpr std::array<std::string_view, 2> match_names{
    "anisotropicminimumdissipation",
    "amd"};

  [[nodiscard]] static AnisotropicMinimumDissipationOptions
  parse_from(ConfigTable const& config);
};

/// @brief Options for the anisotropic minimum dissipation model
/**
 * The parameters for this model may be specified using a table in the config
 * file as:
 * ```toml
 * model = "anisotropicminimumdissipation" # or "amd"
 * PoincareCxy = 0.083333333
 * PoincareCz = 0.3333
 * horizontal_filter_width_scale = 1.0
 * debug_C0 = false
 * ```
 *
 * All parameters are optional, and the defaults are listed in the source file.
 *
 * The "model" is handled by the solver, and the other parameters are parsed
 * by the @ref parse_from function.
 */
struct AnisotropicMinimumDissipationScalarOptions
{
  using ModelType = AnisotropicMinimumDissipationScalar;
  // values in the config file to match
  static constexpr std::array<std::string_view, 2> match_names{
    "anisotropicminimumdissipation",
    "amd"};

  [[nodiscard]] static AnisotropicMinimumDissipationScalarOptions
  parse_from(ConfigTable const& config);
};

} // namespace alps::solver
