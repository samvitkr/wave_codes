#pragma once

#include <toml++/impl/version.hpp>

#include <filesystem>
#include <memory>
#include <string_view>
#include <vector>

#if !defined(ALPS_CONCAT)
#define ALPS_CONCAT_1(x, y) x##y
#define ALPS_CONCAT(x, y)   ALPS_CONCAT_1(x, y)
#define ALPS_INLINE_CONCAT_MACRO
#endif

namespace toml {
inline namespace ALPS_CONCAT(v, TOML_LIB_MAJOR) {
class table;
}
} // namespace toml

#if defined(ALPS_INLINE_CONCAT_MACRO)
#undef ALPS_CONCAT
#undef ALPS_CONCAT_1
#endif

namespace alps {

/// @brief A table of configuration values
/// @note This is a plmpl wrapper around toml::table and provides access methods
/// to interact with TOML++ library. The function templates are not defined in
/// the header file to workaround nvcc bug with TOML++.
class ConfigTable
{
 public:
  static ConfigTable parse_from_file(std::filesystem::path filename);

  static ConfigTable parse_from_string(std::string_view str);

  ConfigTable(ConfigTable const&);

  ConfigTable(ConfigTable&&) noexcept;

  ConfigTable& operator=(ConfigTable const&);

  ConfigTable& operator=(ConfigTable&&) noexcept;

  explicit ConfigTable(std::unique_ptr<toml::table> table);

  explicit ConfigTable(toml::table const& table);

  explicit ConfigTable(toml::table&& table);

  explicit operator toml::table const&() const
  {
    if (table_ == nullptr) {
      throw std::runtime_error("ConfigTable is empty");
    }
    return *table_;
  }

  /// @brief Returns true if the table contains the given path
  bool contains(std::string_view path) const;

  /// @brief Extract a subtable from the current table
  ConfigTable extract_table(std::string_view path) const;

  /// Extract an array of tables
  std::vector<ConfigTable> extract_array_of_tables(std::string_view path) const;

  /// @brief Get a value from the table (supports single value type T and
  /// vector<T>)
  template<class T>
  T get_value(std::string_view path) const;

  /// @brief Get a value from the table with a fallback default
  template<class T>
  T get_value_or(std::string_view path, T default_value) const;

  /// @brief Print the table to a string
  std::string to_string() const;

  ~ConfigTable();

 private:
  ConfigTable();

  std::unique_ptr<toml::table> table_;
};

} // namespace alps
