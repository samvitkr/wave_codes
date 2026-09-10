//
// Created by xuanx004 on 3/8/23.
//

#pragma once

#include "anisotropic_minimum_dissipation_options.h"

#include <common/container/matrix_field.h>
#include <solvers/field/flow_field.h>
#include <solvers/source_terms/boussinesq.h>

#include <string>
#include <string_view>

namespace alps::solver {

class ChannelFlowSolverAB2;
class FlowOverWaveSolver;
class FreeSurfaceSolver;

/// @brief Class describing an anisotropic minimum dissipation model
/** see @ref AnisotropicMinimumDissipationOptions for the parameters */
struct AnisotropicMinimumDissipation
{
  explicit AnisotropicMinimumDissipation(
    AnisotropicMinimumDissipationOptions options);

  std::string show() const;

  static constexpr std::string_view name = "anisotropic minimum dissipation";

  AnisotropicMinimumDissipationOptions options_{};

  mutable Tensor33Field<Real***> grad_u;

  bool keep_grad{false};

  void cleanup_grad() const;
};

struct AnisotropicMinimumDissipationScalar
{
  explicit AnisotropicMinimumDissipationScalar(
    AnisotropicMinimumDissipation& sgs_model);

  std::string show() const;

  static constexpr std::string_view name = "anisotropic minimum dissipation";

  AnisotropicMinimumDissipation const* sgs_model_;
};

void calculate_eddy_viscosity(
  HaloView<Real***> const&                     nu,
  AnisotropicMinimumDissipation const&         sgs_model,
  Tensor33Field<Real***> const&                grad_u,
  FlowField const&                             flow,
  BoussinesqForce<ChannelFlowSolverAB2> const* buoyancy);

void calculate_eddy_viscosity(
  HaloView<Real***> const&                   nu,
  AnisotropicMinimumDissipation const&       sgs_model,
  Tensor33Field<Real***> const&              grad_u,
  FlowOverWaveField const&                   flow,
  BoussinesqForce<FlowOverWaveSolver> const* buoyancy);

void calculate_eddy_viscosity(
  HaloView<Real***> const&                  nu,
  AnisotropicMinimumDissipation const&      sgs_model,
  Tensor33Field<Real***> const&             grad_u,
  FreeSurfaceFlowField const&               flow,
  BoussinesqForce<FreeSurfaceSolver> const* buoyancy);

namespace detail {
template<bool has_buoyancy>
void calculate_eddy_viscosity_impl(
  HaloView<Real***> const&             nu,
  AnisotropicMinimumDissipation const& sgs_model,
  Tensor33Field<Real***> const&        grad_u,
  Mesh const&                          mesh,
  double                               time,
  Vector3Field<Real***> const*         grad_b);
} // namespace detail

/*
 * @note halo cells of grad_u should be updated before calling this function
 */
void calculate_eddy_diffusivity(
  HaloView<Real***> const&                   De,
  AnisotropicMinimumDissipationScalar const& scalar_sgs_model,
  ScalarField const&                         f,
  Vector3Field<Real***> const&               grad_f,
  FlowField const&                           flow);

void calculate_eddy_diffusivity(
  HaloView<Real***> const&                   De,
  AnisotropicMinimumDissipationScalar const& scalar_sgs_model,
  ScalarField const&                         f,
  Vector3Field<Real***> const&               grad_f,
  FlowOverWaveField const&                   flow);

void calculate_eddy_diffusivity(
  HaloView<Real***> const&                   De,
  AnisotropicMinimumDissipationScalar const& scalar_sgs_model,
  ScalarField const&                         f,
  Vector3Field<Real***> const&               grad_f,
  FreeSurfaceFlowField const&                flow);
} // namespace alps::solver
