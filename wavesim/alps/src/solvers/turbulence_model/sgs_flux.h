//
// Created by xuananqing on 3/11/23.
//

#pragma once

#include <common/container/matrix_field.h>
#include <common/container/view_types.h>
#include <solvers/field/flow_field.h>

namespace alps::solver {
void add_sgs_momentum_fluxes_from_gradu(Tensor33Field<Real***> const& fluxes,
                                        FlowField const&              flow,
                                        HaloView<Real***> const&      nu_t,
                                        Real                          nu,
                                        Vector3Field<Real***> const&  grad_u,
                                        Vector3Field<Real***> const&  grad_v,
                                        Vector3Field<Real***> const&  grad_w,
                                        Real gamma = 1);

void add_sgs_momentum_fluxes_from_gradu(Tensor33Field<Real***> const& fluxes,
                                        FlowOverWaveField const&      flow,
                                        HaloView<Real***> const&      nu_t,
                                        Vector3Field<Real***> const&  grad_u,
                                        Vector3Field<Real***> const&  grad_v,
                                        Vector3Field<Real***> const&  grad_w);

void add_sgs_momentum_fluxes_from_gradu(Tensor33Field<Real***> const& fluxes,
                                        FreeSurfaceFlowField const&   flow,
                                        HaloView<Real***> const&      nu_t,
                                        Vector3Field<Real***> const&  grad_u,
                                        Vector3Field<Real***> const&  grad_v,
                                        Vector3Field<Real***> const&  grad_w);

void add_sgs_momentum_fluxes_from_Sij(Tensor33Field<Real***> const&     fluxes,
                                      FlowField const&                  flow,
                                      HaloView<Real***> const&          nu_t,
                                      SymmTensor33Field<Real***> const& Sij);

void add_sgs_momentum_fluxes_from_Sij(Tensor33Field<Real***> const&     fluxes,
                                      FlowOverWaveField const&          flow,
                                      HaloView<Real***> const&          nu_t,
                                      SymmTensor33Field<Real***> const& Sij);

void add_sgs_momentum_fluxes_from_Sij(Tensor33Field<Real***> const&     fluxes,
                                      FreeSurfaceFlowField const&       flow,
                                      HaloView<Real***> const&          nu_t,
                                      SymmTensor33Field<Real***> const& Sij);
} // namespace alps::solver
