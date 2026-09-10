//
// Created by xuananqing on 3/29/23.
//

#pragma once

#include <common/base/logging_fwd.h>
#include <common/container/view_types.h>
#include <common/real_type.h>

#include <string>

namespace alps::solver {
class Mesh;

namespace detail {
void report_nu_t_stats(MDView<Real const***> const& nu_t,
                       double                       time,
                       Mesh const&                  mesh,
                       std::string                  logger_name,
                       RotatingFileSinkConfig       file_config);

} // namespace detail
} // namespace alps::solver
