#include "symlinks.h"

namespace alps::apps {
void force_create_symlink(std::filesystem::path const target,
                          std::filesystem::path const link)
{
  namespace fs = std::filesystem;
  auto status  = fs::symlink_status(link);
  if (fs::is_regular_file(status) || fs::is_symlink(status)) {
    fs::remove(link);
  }
  fs::create_symlink(target, link);
}
} // namespace alps::apps
