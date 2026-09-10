//
// Created by xuanx004 on 7/12/24.
//

#include "scalar_bc_types.h"

#include <fmt/format.h>

#include <string>

namespace alps::solver {

std::unique_ptr<ScalarBC> ScalarBC::get_updated_bc(double /*time*/)
{
  return {};
}

void ScalarBC::load_from(const ConfigTable& /*config*/) {}

std::string ConstantDirichletBC::info() const
{
  return fmt::format("ConstantDirichletBC(Cb={})", value_);
}

void ConstantDirichletBC::load_from(const ConfigTable& config)
{
  value_ = config.get_value_or<Real>("Cb", 0.0);
}

std::string ConstantGradientBC::info() const
{
  return fmt::format("ConstantGradientBC(grad={})", grad_);
}

void ConstantGradientBC::load_from(const ConfigTable& config)
{
  grad_ = config.get_value_or<Real>("grad", 0.0);
}

std::string ConstantFluxBC::info() const
{
  return fmt::format("ConstantFluxBC(flux={})", flux_);
}

void ConstantFluxBC::load_from(const ConfigTable& config)
{
  flux_ = config.get_value_or<Real>("flux", 0.0);
}

} // namespace alps::solver
