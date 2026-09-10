//
// Created by xuanx004 on 2/17/24.
//

#pragma once

#include <memory>
#include <string>

namespace cxxopts {
class Options;
class ParseResult;
} // namespace cxxopts

namespace alps {

class CLIOptions
{
 public:
  explicit CLIOptions(std::string program_name, std::string help_string = "");

  /// @brief Add an option to the CLI parser
  /// @tparam T The type of the option
  /// @param opts The option string
  /// @param desc The description of the option
  /// @param arg_help Alternative string for the option argument in help message
  /// The option Options must have a long form, and may have a short form. To
  /// specify a command line option, the string "opts" has the form: "l,long" or
  /// "long"
  template<typename T = bool>
  void
  add_option(std::string opts, std::string desc, std::string arg_help = "");

  /// @brief Add an option to the CLI parser with a default value
  /// @tparam T The type of the option
  /// @param opts The option string
  /// @param desc The description of the option
  /// @param default_value The default string of the option
  /// @param arg_help Alternative string for the option argument in help message
  /// The option Options must have a long form, and may have a short form. To
  /// specify a command line option, the string "opts" has the form: "l,long" or
  /// "long"
  template<typename T>
  void add_option_with_default(std::string opts,
                               std::string desc,
                               std::string default_value,
                               std::string arg_help = "");

  /// @brief Parse the command line
  cxxopts::ParseResult const& parse(int argc, char const* const* argv);

  /// @brief Get the value of an option
  template<typename T>
  [[nodiscard]] T get(std::string name) const;

  /// @brief Query whether an option is present
  [[nodiscard]] bool has(std::string name) const;

  /// @brief Count the number of times an option is present
  [[nodiscard]] int count(std::string name) const;

  /// @brief Get the help message
  [[nodiscard]] std::string help() const;

  ~CLIOptions();

  CLIOptions() = delete;

 private:
  std::unique_ptr<cxxopts::Options>     options_;
  std::unique_ptr<cxxopts::ParseResult> result_;
};
} // namespace alps
