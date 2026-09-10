#pragma once

#include <filesystem>

namespace alps::apps {
void force_create_symlink(std::filesystem::path target,
                          std::filesystem::path link);
} // namespace alps::apps
