//
// Created by xuanx004 on 7/14/24.
//

#pragma once

#include <common/container/vector_field.h>
#include <solvers/field/flow_field.h>

namespace alps::solver {
void add_SGS_and_molecular_diffusion_fluxes(Vector3Field<Real***> const& fluxes,
                                            HaloView<Real const***> const& f,
                                            FlowField const&               flow,
                                            HaloView<Real***> const&       nu_D,
                                            Real                           D,
                                            Real gamma = 1);

void add_SGS_and_molecular_diffusion_fluxes(Vector3Field<Real***> const& fluxes,
                                            HaloView<Real const***> const& f,
                                            FlowOverWaveField const&       flow,
                                            HaloView<Real***> const&       nu_D,
                                            Real                           D,
                                            Real gamma = 1);

void add_SGS_and_molecular_diffusion_fluxes(Vector3Field<Real***> const& fluxes,
                                            HaloView<Real const***> const& f,
                                            FreeSurfaceFlowField const&    flow,
                                            HaloView<Real***> const&       nu_D,
                                            Real                           D,
                                            Real gamma = 1);

void add_sgs_and_molecular_diffusion_fluxes_from_gradf(
  Vector3Field<Real***> const&       fluxes,
  Vector3Field<Real const***> const& grad_f,
  FlowField const&                   flow,
  HaloView<Real***> const&           nuD,
  Real                               D,
  Real                               gamma = 1);

void add_sgs_and_molecular_diffusion_fluxes_from_gradf(
  Vector3Field<Real***> const& fluxes,
  Vector3Field<Real***> const& grad_f,
  FlowOverWaveField const&     flow,
  HaloView<Real***> const&     nuD,
  Real                         D,
  Real                         gamma = 1);

void add_sgs_and_molecular_diffusion_fluxes_from_gradf(
  Vector3Field<Real***> const& fluxes,
  Vector3Field<Real***> const& grad_f,
  FreeSurfaceFlowField const&  flow,
  HaloView<Real***> const&     nuD,
  Real                         D,
  Real                         gamma = 1);

} // namespace alps::solver
