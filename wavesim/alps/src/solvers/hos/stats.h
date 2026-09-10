//
// Created by xuanx004 on 5/25/24.
//

#pragma once

#include <common/container/view_types.h>
#include <common/real_type.h>
#include <solvers/hos/field.h>
#include <spectral/spectral_fwd.h>

namespace alps::solver::hos {

bool check_validity(MDView<Real const**> const& eta, Grid const& grid);

void report_wave_stats(HOSField const& field, double time);

Real compute_eta_rms(MDView<Real const**> const& eta, Grid const& grid);

MDView<Real**, Kokkos::HostSpace>
compute_eta_density_2d(MDView<Real const**> const& eta, Grid const& grid);

} // namespace alps::solver::hos
