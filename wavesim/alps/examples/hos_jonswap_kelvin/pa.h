#pragma once

#include <common/container/view_types.h>
#include <common/real_type.h>
#include <spectral/spectral_fwd.h>

namespace alps {

void radial_moving_pressure(MDView<Real**>                       pa,
                            double                               time,
                            double                               Pmax,
                            double                               R,
                            std::pair<double, double>            x0_initial,
                            double                               U,
                            Grid const&                          grid,
                            Kokkos::DefaultExecutionSpace const& space);
} // namespace alps
