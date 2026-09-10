//
// Created by xuanx004 on 7/14/24.
//

#pragma once

#include "constant_schmidt_options.h"

#include <solvers/field/flow_field.h>

#include <string>
#include <string_view>

namespace alps::solver {
struct ConstantTurbulentSc
{
  explicit ConstantTurbulentSc(ConstantTurbulentScOptions options);

  std::string show() const;

  static constexpr std::string_view name = "constant Sc_t/Pr_t";

  ConstantTurbulentScOptions options_{};
};

void calculate_eddy_diffusivity(HaloView<Real***> const&   nuD,
                                ConstantTurbulentSc const& sgs_model,
                                FlowField const&           flow);

void calculate_eddy_diffusivity(HaloView<Real***> const&   nuD,
                                ConstantTurbulentSc const& sgs_model,
                                FlowOverWaveField const&   flow);

void calculate_eddy_diffusivity(HaloView<Real***> const&    nuD,
                                ConstantTurbulentSc const&  sgs_model,
                                FreeSurfaceFlowField const& flow);
} // namespace alps::solver
