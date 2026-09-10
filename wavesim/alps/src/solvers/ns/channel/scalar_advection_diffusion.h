#pragma once

#include <common/container/matrix_field.h>
#include <solvers/field/flow_field.h>

namespace alps::solver {
void calculate_advection_fluxes(const Vector3Field<Real***>&   fluxes,
                                const HaloView<Real const***>& f,
                                const FlowField&               flow);

void add_diffusion_fluxes(Vector3Field<Real***> const&   fluxes,
                          const HaloView<Real const***>& f,
                          const FlowField&               flow,
                          Real                           D,
                          Real                           gamma = 1);
} // namespace alps::solver
