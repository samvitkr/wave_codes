//
// Created by xuananqing on 3/7/23.
//

#pragma once

#include "smagorinsky_options.h"

#include <common/container/matrix_field.h>
#include <solvers/field/flow_field.h>

#include <string>
#include <string_view>

namespace alps::solver {

/// @brief Class describing a constant Smagorinsky model
/** see @ref ConstantSmagorinskyOptions for the parameters */
struct ConstantSmagorinsky
{
  explicit ConstantSmagorinsky(ConstantSmagorinskyOptions options);

  std::string show() const;

  static constexpr std::string_view name = "constant Smagorinsky";

  ConstantSmagorinskyOptions options_{};
};

void calculate_eddy_viscosity(HaloView<Real***> const&          nu,
                              ConstantSmagorinsky const&        sgs_model,
                              SymmTensor33Field<Real***> const& Sij,
                              FlowField const&                  flow);

void calculate_eddy_viscosity(HaloView<Real***> const&          nu,
                              ConstantSmagorinsky const&        sgs_model,
                              SymmTensor33Field<Real***> const& Sij,
                              FlowOverWaveField const&          flow);

void calculate_eddy_viscosity(HaloView<Real***> const&          nu,
                              ConstantSmagorinsky const&        sgs_model,
                              SymmTensor33Field<Real***> const& Sij,
                              FreeSurfaceFlowField const&       flow);

} // namespace alps::solver
