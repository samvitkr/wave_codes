//
// Created by xuanx004 on 5/26/23.
//

#include "wall_models.h"

#include <common/program_options/config_table.h>
#include <common/utils/to_lower_case.h>

#include <fmt/format.h>

namespace alps::solver {

LogLawWallModelOptions
LogLawWallModelOptions::parse_from(ConfigTable const& config)
{
  LogLawWallModelOptions options;
  options.kappa = config.get_value_or<double>("kappa", 0.4);
  options.z0    = config.get_value<double>("z0");

  // Validate the wall model parameters
  if (options.z0 <= 0) {
    throw std::runtime_error(
      fmt::format("Invalid wall model z0={}", options.z0));
  }
  if (options.kappa <= 0) {
    throw std::runtime_error(
      fmt::format("Invalid wall model kappa={}", options.kappa));
  }

  return options;
}

namespace {
template<typename... Ts>
std::variant<std::monostate, Ts...>
MatchOptionVariant(std::variant<std::monostate, Ts...> const& /*unused*/,
                   std::string const& name,
                   ConfigTable const& config)
{
  std::variant<std::monostate, Ts...> result = std::monostate{};

  auto try_match = [&](auto&& type) {
    if (!std::holds_alternative<std::monostate>(result)) return;

    using Type          = std::decay_t<decltype(type)>;
    auto constexpr list = Type::match_names;
    auto itr            = list.cbegin();
    for (; itr != list.cend(); ++itr) {
      if (*itr == name) {
        break;
      }
    }
    if (itr != list.cend()) {
      result = Type::parse_from(config);
    }
  };

  // fold expression: test each type in the variant until a match
  (try_match(Ts{}), ...);
  return result;
}
} // namespace

WallLayerModelOptionVariants
WallLayerModel::parse_options(ConfigTable const& config)
{
  auto const input_model_name =
    config.get_value_or<std::string>("type", "none");
  auto const model_name = to_lower_case(input_model_name);

  if (model_name == "none") {
    return std::monostate{};
  }

  // loop through all types and match the one that is in the config
  auto const option =
    MatchOptionVariant(WallLayerModelOptionVariants{}, model_name, config);
  if (std::holds_alternative<std::monostate>(option)) {
    // throw an exception if no match is found
    throw std::runtime_error("Unknown wall model: " + input_model_name);
  }
  return option;
}

} // namespace alps::solver
