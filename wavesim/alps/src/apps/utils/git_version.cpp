//
// Created by xuananqing on 5/23/23.
//

#include "git_version.h"

#include <git.h>

namespace alps::apps {
std::string git_version(std::string prefix)
{
  if (!git_IsPopulated()) {
    return prefix.empty() ? "unknown_version" : prefix + "_unknown_version";
  }

  const std::string version{git_Describe()};
  return prefix.empty() ? version : prefix + "_" + version;
}
} // namespace alps::apps
