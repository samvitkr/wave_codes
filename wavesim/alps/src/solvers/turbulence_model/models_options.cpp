//
// Created by xuananqing on 7/14/24.
//

#include "models_options.h"

#include <common/utils/to_lower_case.h>

#include <string>

namespace alps::solver {

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

SGSModelOptionVariants parse_sgs_model_options(ConfigTable const& config)
{
  auto const input_model_name = config.get_value<std::string>("model");
  auto const model_name       = to_lower_case(input_model_name);

  if (model_name == "none") {
    return std::monostate{};
  }

  // loop through all types and match the one that is in the config
  auto const option =
    MatchOptionVariant(SGSModelOptionVariants{}, model_name, config);
  if (std::holds_alternative<std::monostate>(option)) {
    // throw an exception if no match is found
    throw std::runtime_error("Unknown LES SGS model: " + input_model_name);
  }
  return option;
}

ScalarSGSModelOptionVariants
parse_scalar_sgs_model_options(ConfigTable const& config)
{
  auto const input_model_name = config.get_value<std::string>("model");
  auto const model_name       = to_lower_case(input_model_name);

  if (model_name == "none") {
    return std::monostate{};
  }

  auto const option =
    MatchOptionVariant(ScalarSGSModelOptionVariants{}, model_name, config);
  if (std::holds_alternative<std::monostate>(option)) {
    throw std::runtime_error("Unknown scalar SGS model: " + input_model_name);
  }
  return option;
}

} // namespace alps::solver
