//
// Created by xuananqing on 3/7/23.
//

#pragma once

#include "dynamic_smagorinsky_options.h"
#include <common/container/matrix_field.h>
#include <common/container/view_types.h>
#include <solvers/field/flow_field.h>

#include <string>
#include <string_view>

namespace alps::solver {

/// @brief Class describing a dynamic Smagorinsky model
/** see @ref DynamicSmagorinskyOptions for the parameters */
struct DynamicSmagorinsky
{
#if (ALPS_SGS_USE_LOWER_PRECISION)
  using LM_floating_t = float;
#else
  using LM_floating_t = Real;
#endif

  explicit DynamicSmagorinsky(DynamicSmagorinskyOptions options);

  std::string show() const;

  static constexpr std::string_view name = "dynamic Smagorinsky";

  DynamicSmagorinskyOptions options_{};

  mutable MDView<Real*> C0_Delta2;

  mutable HaloView<Real***>          Smag;
  mutable HaloView<LM_floating_t***> Smag_test;

  bool keep_Smag{false};

  void cleanup_Smag() const;
};

struct DynamicSmagorinskyScalar
{
  explicit DynamicSmagorinskyScalar(DynamicSmagorinsky& sgs_model);

  std::string show() const;

  static constexpr std::string_view name = "dynamic Smagorinsky";

  mutable MDView<Real*> C0_Delta2;

  DynamicSmagorinsky const* sgs_model_;
};

void calculate_eddy_viscosity(HaloView<Real***> const&          nu,
                              DynamicSmagorinsky const&         sgs_model,
                              SymmTensor33Field<Real***> const& Sij,
                              FlowField const&                  flow,
                              bool                              skip_update_C0);

void calculate_eddy_viscosity(HaloView<Real***> const&          nu,
                              DynamicSmagorinsky const&         sgs_model,
                              SymmTensor33Field<Real***> const& Sij,
                              FlowOverWaveField const&          flow,
                              bool                              skip_update_C0);

void calculate_eddy_viscosity(HaloView<Real***> const&          nu,
                              DynamicSmagorinsky const&         sgs_model,
                              SymmTensor33Field<Real***> const& Sij,
                              FreeSurfaceFlowField const&       flow,
                              bool                              skip_update_C0);

void calculate_eddy_diffusivity(
  HaloView<Real***> const&        nuD,
  DynamicSmagorinskyScalar const& scalar_sgs_model,
  ScalarField const&              scalar,
  Vector3Field<Real***> const&    grad_f,
  FlowField const&                flow,
  bool                            skip_update_C0);

void calculate_eddy_diffusivity(
  HaloView<Real***> const&        nuD,
  DynamicSmagorinskyScalar const& scalar_sgs_model,
  ScalarField const&              scalar,
  Vector3Field<Real***> const&    grad_f,
  FlowOverWaveField const&        flow,
  bool                            skip_update_C0);

void calculate_eddy_diffusivity(
  HaloView<Real***> const&        nuD,
  DynamicSmagorinskyScalar const& scalar_sgs_model,
  ScalarField const&              scalar,
  Vector3Field<Real***> const&    grad_f,
  FreeSurfaceFlowField const&     flow,
  bool                            skip_update_C0);

} // namespace alps::solver
