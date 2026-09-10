//
// Created by xuananqing on 5/23/23.
//

#pragma once

#include <string>

namespace alps::apps {

/// Returns a version string based on git commit hash.
//* @param prefix A prefix to prepend to the version string. */
std::string git_version(std::string prefix = "");

} // namespace alps::apps
