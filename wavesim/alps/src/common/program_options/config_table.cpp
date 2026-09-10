#include "config_table.h"

#include <fmt/core.h>
#include <fmt/ostream.h>
#include <toml++/toml.h>

namespace alps {

namespace {
std::string error_msg_not_found(std::string_view path)
{
  return fmt::format("key '{}' not found in table", path);
}

std::string error_msg_incorrect_type(std::string_view path,
                                     std::string_view type)
{
  return fmt::format("key '{}' is not {}", path, type);
}
} // namespace

template<typename Container>
struct is_vector : std::false_type
{};

template<typename... Ts>
struct is_vector<std::vector<Ts...>> : std::true_type
{};

ConfigTable::ConfigTable(toml::table const& table)
  : table_{std::make_unique<toml::table>(table)}
{}

ConfigTable::ConfigTable(ConfigTable const& other)
  : table_{std::make_unique<toml::table>(*other.table_)}
{}

ConfigTable::ConfigTable(std::unique_ptr<toml::table> table)
  : table_{std::move(table)}
{}

ConfigTable::ConfigTable(toml::table&& table)
  : table_{std::make_unique<toml::table>(std::move(table))}
{}

ConfigTable& ConfigTable::operator=(ConfigTable const& other)
{
  if (this != &other) {
    table_ = std::make_unique<toml::table>(*other.table_);
  }
  return *this;
}

ConfigTable& ConfigTable::operator=(ConfigTable&& other) noexcept
{
  if (this != &other) {
    table_ = std::move(other.table_);
  }
  return *this;
}

ConfigTable::ConfigTable()                       = default;
ConfigTable::ConfigTable(ConfigTable&&) noexcept = default;
ConfigTable::~ConfigTable()                      = default;

bool ConfigTable::contains(std::string_view path) const
{
  return table_->at_path(path).operator bool();
}

ConfigTable ConfigTable::extract_table(std::string_view path) const
{
  const auto node = table_->at_path(path);
  if (!node) {
    throw std::runtime_error(error_msg_not_found(path));
  }
  if (!node.is<toml::table>()) {
    throw std::runtime_error(error_msg_incorrect_type(path, "a table"));
  }
  return ConfigTable(*node.as_table());
}

ConfigTable ConfigTable::parse_from_file(std::filesystem::path filename)
{
  ConfigTable result;
  try {
    auto config = std::make_unique<toml::table>();
    *config     = toml::parse_file(filename.c_str());
    result.table_.swap(config);
  } catch (const toml::parse_error& err) {
    auto msg = fmt::format("{} ({}, {})",
                           err.description(),
                           fmt::streamed(*err.source().path),
                           fmt::streamed(err.source().begin));
    throw std::runtime_error(msg);
  }
  return result;
}

ConfigTable ConfigTable::parse_from_string(std::string_view str)
{
  ConfigTable result;
  try {
    auto config = std::make_unique<toml::table>();
    *config     = toml::parse(str);
    result.table_.swap(config);
  } catch (const toml::parse_error& err) {
    auto msg = fmt::format("{} ({}, {})",
                           err.description(),
                           fmt::streamed(*err.source().path),
                           fmt::streamed(err.source().begin));
    throw std::runtime_error(msg);
  }
  return result;
}

std::string ConfigTable::to_string() const
{
  return fmt::format("{}", fmt::streamed(*table_));
}

template<class T>
T ConfigTable::get_value(std::string_view path) const
{
  if constexpr (is_vector<T>::value) {
    using element_t = typename T::value_type;
    static_assert(toml::impl::is_native<element_t>
                    || toml::impl::can_represent_native<element_t>
                    || toml::impl::can_partially_represent_native<element_t>,
                  "element_t not supported");

    const auto value = table_->at_path(path);
    if (!value) {
      throw std::runtime_error(error_msg_not_found(path));
    }
    const auto* array = value.as_array();
    if (array == nullptr) {
      throw std::runtime_error(error_msg_incorrect_type(path, "an array"));
    }
    T vec;
    for (const auto& element : *array) {
      const auto array_val = element.template value<element_t>();
      if (array_val) {
        vec.emplace_back(*array_val);
      } else {
        auto msg = fmt::format("error parsing values in '{}'", path);
        throw std::runtime_error(msg);
      }
    }
    return vec;
  } else if constexpr ((toml::impl::is_native<T>
                        || toml::impl::can_represent_native<T>
                        || toml::impl::can_partially_represent_native<T>)
                       && !toml::impl::is_cvref<T>) {
    const auto value = table_->at_path(path);
    if (!value) {
      throw std::runtime_error(error_msg_not_found(path));
    }
    const auto parsed_val = value.value<T>();
    if (!parsed_val) {
      throw std::runtime_error(
        error_msg_incorrect_type(path, "the requested type"));
    }
    return *parsed_val;
  } else {
    static_assert(!std::is_same_v<T, T>, "type not supported");
  }
}

template<class T>
T ConfigTable::get_value_or(std::string_view path, T default_value) const
{
  return table_->at_path(path).value_or<T>(std::move(default_value));
}

// Explicit instantiations
template float  ConfigTable::get_value<float>(std::string_view) const;
template double ConfigTable::get_value<double>(std::string_view) const;
template int    ConfigTable::get_value<int>(std::string_view) const;
template std::int64_t
  ConfigTable::get_value<std::int64_t>(std::string_view) const;
template std::string
              ConfigTable::get_value<std::string>(std::string_view) const;
template bool ConfigTable::get_value<bool>(std::string_view) const;
template std::vector<float>
  ConfigTable::get_value<std::vector<float>>(std::string_view) const;
template std::vector<double>
  ConfigTable::get_value<std::vector<double>>(std::string_view) const;
template std::vector<int>
  ConfigTable::get_value<std::vector<int>>(std::string_view) const;
template std::vector<std::int64_t>
  ConfigTable::get_value<std::vector<std::int64_t>>(std::string_view) const;

template float  ConfigTable::get_value_or<float>(std::string_view, float) const;
template double ConfigTable::get_value_or<double>(std::string_view,
                                                  double) const;
template int    ConfigTable::get_value_or<int>(std::string_view, int) const;
template std::string ConfigTable::get_value_or<std::string>(std::string_view,
                                                            std::string) const;
template bool ConfigTable::get_value_or<bool>(std::string_view, bool) const;

std::vector<ConfigTable>
ConfigTable::extract_array_of_tables(std::string_view path) const
{
  const auto node = table_->at_path(path);
  if (!node) {
    throw std::runtime_error(error_msg_not_found(path));
  }
  const auto* array = node.as_array();
  if (array == nullptr || !array->is_array_of_tables()) {
    throw std::runtime_error(
      error_msg_incorrect_type(path, "an array of tables"));
  }

  std::vector<ConfigTable> results;
  array->for_each(
    [&results](toml::table const& el) { results.emplace_back(el); });
  return results;
}

} // namespace alps
