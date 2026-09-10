#pragma once

#include <common/container/view_types.h>
#include <common/real_type.h>
#include <spectral/spectral_fwd.h>

namespace alps::solver::hos {

void calc_evolution_rhs(MDView<Real** [2]>                   F_t,
                        MDView<const Real** [2]>             F,
                        MDView<const Real**>                 ws,
                        MDView<const Real**>                 pa,
                        Real                                 Fr2,
                        Real                                 We,
                        Grid const&                          grid,
                        Kokkos::DefaultExecutionSpace const& space);

void calc_surface_tension(MDView<Real**>                       p_st,
                          MDView<const Real** [1]>             eta_x,
                          MDView<const Real** [1]>             eta_y,
                          Real                                 We,
                          Grid const&                          grid,
                          Kokkos::DefaultExecutionSpace const& space);
} // namespace alps::solver::hos
