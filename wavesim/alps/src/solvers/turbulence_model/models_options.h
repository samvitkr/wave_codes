//
// This file contains classes for options of turbulence models.
// Created by xuanx004 on 1/23/24.
//

#pragma once

#include <solvers/turbulence_model/anisotropic_minimum_dissipation_options.h>
#include <solvers/turbulence_model/dynamic_smagorinsky_options.h>
#include <solvers/turbulence_model/smagorinsky_options.h>

#include <solvers/turbulence_model/constant_schmidt_options.h>

#include <variant>

namespace alps::solver {

// Variant type for all possible SGS model options, used to store the options
// for the SGS model
using SGSModelOptionVariants =
  std::variant<std::monostate,
               ConstantSmagorinskyOptions,
               DynamicSmagorinskyOptions,
               AnisotropicMinimumDissipationOptions>;

using ScalarSGSModelOptionVariants =
  std::variant<std::monostate,
               ConstantTurbulentScOptions,
               DynamicSmagorinskyScalarOptions,
               AnisotropicMinimumDissipationScalarOptions>;

SGSModelOptionVariants parse_sgs_model_options(ConfigTable const& config);

ScalarSGSModelOptionVariants
parse_scalar_sgs_model_options(ConfigTable const& config);

} // namespace alps::solver
