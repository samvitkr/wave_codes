//
// Created by xuanx004 on 1/23/24.
//

#pragma once

#include <common/program_options/config_table.h>

#include <array>
#include <string_view>

namespace alps::solver {

struct DynamicSmagorinsky;
struct DynamicSmagorinskyScalar;

/** @brief Options for the dynamic Smagorinsky model
 * The parameters for this model may be specified using a table in the config
 * file as:
 * ```toml
 * model = "dynamicsmagorinsky" # or "DS"
 * debug_C0 = false
 * ```
 *
 * All parameters are optional, and the defaults are listed in the source file.
 *
 * The "model" is handled by the solver, and the other parameters are parsed
 * by the @ref parse_from function.
 */
struct DynamicSmagorinskyOptions
{
  double horizontal_filter_width_scale{1.0};

  bool enable_debug_C0{false};

  using ModelType = DynamicSmagorinsky;
  // values in the config file to match
  static constexpr std::array<std::string_view, 2> match_names{
    "dynamicsmagorinsky",
    "ds"};

  [[nodiscard]] static DynamicSmagorinskyOptions
  parse_from(ConfigTable const& config);
};

/** @brief Options for the dynamic Smagorinsky model for scalars
 * The parameters for this model may be specified using a table in the config
 * file as:
 * ```toml
 * model = "dynamicsmagorinsky" # or "DS"
 * ```
 *
 * All parameters are optional, and the defaults are listed in the source file.
 *
 * The "model" is handled by the solver, and the other parameters are parsed
 * by the @ref parse_from function.
 */
struct DynamicSmagorinskyScalarOptions
{
  using ModelType = DynamicSmagorinskyScalar;
  // values in the config file to match
  static constexpr std::array<std::string_view, 2> match_names{
    "dynamicsmagorinsky",
    "ds"};

  [[nodiscard]] static DynamicSmagorinskyScalarOptions
  parse_from(ConfigTable const& config);
};

} // namespace alps::solver
