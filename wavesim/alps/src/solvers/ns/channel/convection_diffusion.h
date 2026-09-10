#pragma once

#include <common/container/matrix_field.h>
#include <solvers/field/flow_field.h>

namespace alps::solver {

void calculate_convection_fluxes(const Tensor33Field<Real***>& fluxes,
                                 const FlowField&              flow);

void add_viscous_fluxes_from_u(Tensor33Field<Real***> const& fluxes,
                               const FlowField&              flow,
                               Real                          nu,
                               Real                          gamma = 1);

} // namespace alps::solver
