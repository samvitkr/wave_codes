#pragma once

#include <common/base/logging_fwd.h>
#include <common/container/view_types.h>
#include <common/real_type.h>

namespace alps::solver {
class Mesh;
struct ConstantSmagorinsky;
struct DynamicSmagorinsky;
struct DynamicSmagorinskyScalar;
namespace detail {

void report_smagorinsky_C0(ConstantSmagorinsky const&              sgs_model,
                           MDView<Real*, Kokkos::HostSpace> const& scaled_delta,
                           double                                  time,
                           Mesh const&                             mesh,
                           std::string const&                      logger_name,
                           RotatingFileSinkConfig                  file_config);

void report_dynamic_smagorinsky_C0(DynamicSmagorinsky const& sgs_model,
                                   double                    time,
                                   Mesh const&               mesh,
                                   std::string const&        logger_name,
                                   RotatingFileSinkConfig    file_config);

void report_dynamic_smagorinsky_C0(DynamicSmagorinskyScalar const& sgs_model,
                                   double                          time,
                                   Mesh const&                     mesh,
                                   std::string const&              logger_name,
                                   RotatingFileSinkConfig          file_config);

} // namespace detail
} // namespace alps::solver
