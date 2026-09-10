//
// Created by xuananqing on 6/14/24.
//

#pragma once

#include <type_traits>

namespace alps::solver {

template<class FieldType>
struct is_curvilinear_field : std::false_type
{};

} // namespace alps::solver
