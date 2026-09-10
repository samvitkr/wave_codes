//
// Created by xuanx004 on 7/15/24.
//

#include "to_lower_case.h"

#include <algorithm>

namespace alps {

std::string to_lower_case(std::string const& str)
{
  std::string result = str;
  std::transform(result.begin(),
                 result.end(),
                 result.begin(),
                 [](unsigned char c) { return std::tolower(c); });
  return result;
}

} // namespace alps
