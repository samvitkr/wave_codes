//
// Created by xuananqing on 10/25/22.
//

#include "bc_types.h"

#include <common/kokkos_abstraction/pool_space.h>

#include <fmt/format.h>

#include <string>

namespace alps::solver {

std::unique_ptr<VelocityBC> VelocityBC::get_updated_bc(double /*time*/)
{
  return {};
}

void VelocityBC::load_from(const ConfigTable& /*config*/) {}

std::string NoSlipWall::info() const
{
  return fmt::format("NoSlipWall(u={}, v={}, w={})", u_, v_, w_);
}

void NoSlipWall::load_from(const ConfigTable& config)
{
  using bc_value_t = decltype(u_);
  u_ = static_cast<bc_value_t>(config.get_value_or<double>("u", 0.0));
  v_ = static_cast<bc_value_t>(config.get_value_or<double>("v", 0.0));
  w_ = static_cast<bc_value_t>(config.get_value_or<double>("w", 0.0));
}

NoSlipWallVarying::NoSlipWallVarying(int nx, int ny, bool is_time_dependent_)
  : VelocityBC(is_time_dependent_)
  , u_(Kokkos::View<Real**, Kokkos::LayoutLeft, default_memory_pool>("ub",
                                                                     nx,
                                                                     ny))
  , v_(Kokkos::View<Real**, Kokkos::LayoutLeft, default_memory_pool>(
      "vb",
      u_.layout()))
  , w_(Kokkos::View<Real**, Kokkos::LayoutLeft, default_memory_pool>(
      "wb",
      u_.layout()))
  , eta_(Kokkos::View<Real**, Kokkos::LayoutLeft, default_memory_pool>(
      "eta",
      u_.layout()))
  , eta_t_(Kokkos::View<Real**, Kokkos::LayoutLeft, default_memory_pool>(
      "eta_t",
      u_.layout()))
{}

std::string NoSlipWallVarying::info() const
{
  return "NoSlipWallVarying";
}

std::string GradientWall::info() const
{
  return fmt::format("GradientWall(grad_1={}, grad_2={})", grad_1, grad_2);
}

void GradientWall::load_from(const ConfigTable& config)
{
  using bc_value_t = decltype(grad_1);
  grad_1 = static_cast<bc_value_t>(config.get_value_or<double>("grad_1", 0.0));
  grad_2 = static_cast<bc_value_t>(config.get_value_or<double>("grad_2", 0.0));
}

std::string TangentialStressWall::info() const
{
  return fmt::format("TangentialStressWall(tau_1={}, tau_2={})", tau_1, tau_2);
}

void TangentialStressWall::load_from(const ConfigTable& config)
{
  using bc_value_t = decltype(tau_1);
  tau_1 = static_cast<bc_value_t>(config.get_value_or<double>("tau_1", 0.0));
  tau_2 = static_cast<bc_value_t>(config.get_value_or<double>("tau_2", 0.0));
}

} // namespace alps::solver
