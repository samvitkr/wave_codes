#include "std_cli_opts.h"

#include <cstdint>

namespace alps::apps {

void add_std_cli_options(alps::CLIOptions& options)
{
  options.add_option("h,help", "Print usage");
  options.add_option<std::int64_t>(
    "restart-from", "Specify a restart file by ID", "ID");
  options.add_option_with_default<std::string>(
    "c,config",
    "Specify a configuration file",
    std::string{DEFAULT_CONFIG_FILENAME},
    "filename");
}

} // namespace alps::apps
