#pragma once

#include <common/program_options/cli_options.h>

#include <string_view>

namespace alps::apps {

inline constexpr std::string_view DEFAULT_CONFIG_FILENAME = "param.toml";

/// Add standard CLI options
void add_std_cli_options(alps::CLIOptions& options);

} // namespace alps::apps
