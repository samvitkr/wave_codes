//
// Created by xuanx004 on 2/17/24.
//

#include "cli_options.h"

#define CXXOPTS_NO_REGEX true
#include <cxxopts.hpp>
#undef CXXOPTS_NO_REGEX

#include <cstdint>

namespace alps {

CLIOptions::CLIOptions(std::string program_name, std::string help_string)
  : options_{std::make_unique<cxxopts::Options>(program_name, help_string)}
  , result_{nullptr}
{}

template<typename T>
void CLIOptions::add_option(std::string opts,
                            std::string desc,
                            std::string arg_help)
{
  options_->add_option("", {opts, desc, cxxopts::value<T>(), arg_help});
}

template<typename T>
void CLIOptions::add_option_with_default(std::string opts,
                                         std::string desc,
                                         std::string default_value,
                                         std::string arg_help)
{
  options_->add_option(
    "",
    {opts, desc, cxxopts::value<T>()->default_value(default_value), arg_help});
}

cxxopts::ParseResult const& CLIOptions::parse(int argc, const char* const* argv)
{
  // place the return value of parse in the unique_ptr result_
  result_ = std::make_unique<cxxopts::ParseResult>(options_->parse(argc, argv));
  return *result_;
}

bool CLIOptions::has(std::string option) const
{
  if (!result_) {
    throw std::runtime_error("CLIOptions::parse() must be called before has()");
  }
  return result_->count(option) > 0;
}

int CLIOptions::count(std::string option) const
{
  if (!result_) {
    throw std::runtime_error(
      "CLIOptions::parse() must be called before count()");
  }
  return result_->count(option);
}

template<typename T>
T CLIOptions::get(std::string option) const
{
  if (!result_) {
    throw std::runtime_error("CLIOptions::parse() must be called before get()");
  }
  return result_->operator[](option).as<T>();
}

std::string CLIOptions::help() const
{
  return options_->help();
}

CLIOptions::~CLIOptions() = default;

template void CLIOptions::add_option<bool>(std::string opts,
                                           std::string desc,
                                           std::string arg_help);
template void CLIOptions::add_option<std::string>(std::string opts,
                                                  std::string desc,
                                                  std::string arg_help);
template void CLIOptions::add_option<int>(std::string opts,
                                          std::string desc,
                                          std::string arg_help);
template void CLIOptions::add_option<std::int64_t>(std::string opts,
                                                   std::string desc,
                                                   std::string arg_help);

template void
CLIOptions::add_option_with_default<std::string>(std::string opts,
                                                 std::string desc,
                                                 std::string default_value,
                                                 std::string arg_help);
template void
CLIOptions::add_option_with_default<int>(std::string opts,
                                         std::string desc,
                                         std::string default_value,
                                         std::string arg_help);
// template void
// CLIOptions::add_option_with_default<std::int64_t>(std::string opts,
//                                                   std::string desc,
//                                                   std::string default_value,
//                                                   std::string arg_help);

template std::string  CLIOptions::get<std::string>(std::string option) const;
template int          CLIOptions::get<int>(std::string option) const;
template std::int64_t CLIOptions::get<std::int64_t>(std::string option) const;

} // namespace alps
