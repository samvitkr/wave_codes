//
// Created by xuanx004 on 5/26/23.
//

#pragma once

#include <common/program_options/config_table.h>

#include <array>
#include <string_view>
#include <variant>

namespace alps::solver {

// forward declaration
class LogLawWallModel;

/// @brief Options for the log law wall model
/**
 * The parameters for this model may be specified using a table in the config
 * file as:
 * ```
 * { type = "loglaw", kappa = 0.4, z0 = 0.001 }
 * ```
 *
 * z0 is mandatory, and the default kappa are listed in the source file.
 *
 * The "type" part is handled by the solver, and the other parameters are parsed
 * by the @ref parse_from function.
 */
struct LogLawWallModelOptions
{
  double kappa{0.4}; /// von Kármán constant

  double z0; /// roughness length

  using ModelType = LogLawWallModel;

  static constexpr std::array<std::string_view, 1> match_names{"loglaw"};

  [[nodiscard]] static LogLawWallModelOptions
  parse_from(ConfigTable const& config);
};

using WallLayerModelOptionVariants =
  std::variant<std::monostate, LogLawWallModelOptions>;

class WallLayerModel
{
 public:
  virtual std::string info() const = 0;

  virtual ~WallLayerModel() = default;

  static WallLayerModelOptionVariants parse_options(ConfigTable const& config);
};

} // namespace alps::solver
