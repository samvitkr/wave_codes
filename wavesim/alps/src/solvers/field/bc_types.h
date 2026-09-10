/*
 * @Author: Anqing Xuan
 * @Date: 2020-07-30 02:02:10
 */
#pragma once

#include <common/program_options/config_table.h>
#include <common/real_type.h>

#include <Kokkos_Core.hpp>

#include <memory>
#include <string>

namespace alps {
namespace solver {

/// Base class for boundary conditions of velocities
/** @note This class is not meant to be used directly, but rather to be
 * inherited by specific boundary conditions.
 *  @note The derived classes should contain the parameters and state variables
 * of the BC.
 */
struct VelocityBC
{
  /// Get the information of the boundary condition
  virtual std::string info() const = 0;

  /// Return the boundary conditions at a given time
  virtual std::unique_ptr<VelocityBC> get_updated_bc(double time);

  /// Load configurations from a TOML table
  virtual void load_from(const ConfigTable& config);

  virtual ~VelocityBC() = default;

  bool is_time_dependent{false};

 protected:
  VelocityBC() = default;

  explicit VelocityBC(bool is_time_dependent_)
    : is_time_dependent(is_time_dependent_)
  {}
};

/// @brief A Dirichlet boundary condition where the velocities are scalar values
/**
 * This BC is specified in the config file as a table:
 * ```toml
 * BC.top = { type = "noslipwall", u = 0.0, v = 0.0, w = 0.0 }
 * BC.bottom = { type = "noslipwall" }
 * ```
 * The velocity components `u`, `v`, and `w` are optional and default to 0.
 */
struct NoSlipWall : public VelocityBC
{
  std::string info() const override;

  void load_from(const ConfigTable& config) override;

  Real u_{0};
  Real v_{0};
  Real w_{0};
};

/// @brief A Dirichlet boundary condition where the velocities are 2D fields
struct NoSlipWallVarying : public VelocityBC
{
  [[maybe_unused]] NoSlipWallVarying(int  nx,
                                     int  ny,
                                     bool is_time_dependent_ = false);

  std::string info() const override;

  Kokkos::View<Real**, Kokkos::LayoutLeft> u_;
  Kokkos::View<Real**, Kokkos::LayoutLeft> v_;
  Kokkos::View<Real**, Kokkos::LayoutLeft> w_;

  Kokkos::View<Real**, Kokkos::LayoutLeft> eta_;
  Kokkos::View<Real**, Kokkos::LayoutLeft> eta_t_;
};

/**
 * @brief A Neumann boundary condition where the vertical gradients in two
 * horizontal directions are scalar values
 *
 * This BC is specified in the config file as a table:
 * ```toml
 * BC.top = { type = "gradientwall", grad_1 = 1.0, grad_2 = 0.0 }
 * BC.bottom = { type = "gradientwall" }
 * ```
 * The gradients `grad_1` and `grad_2` are optional and default to 0.
 */
struct GradientWall : public VelocityBC
{
  std::string info() const override;

  void load_from(const ConfigTable& config) override;

  Real grad_1{0};
  Real grad_2{0};
};

/**
 * @brief A Neumann boundary condition where the tangential stresses are scalar
 * values
 *
 * This BC is specified in the config file as a table:
 * ```toml
 * BC.top = { type = "tangentialstresswall", tau_1 = 1.0, tau_2 = 0.0 }
 * ```
 * The tangential stresses `tau_1` and `tau_2` are optional and default to 0.
 * Currently, this BC is only implemented for the top boundary.
 */
struct TangentialStressWall : public VelocityBC
{
  std::string info() const override;

  void load_from(const ConfigTable& config) override;

  Real tau_1{0};
  Real tau_2{0};
};

/// An enum to specify which boundary is applied to
enum class WhichBoundary
{
  Both,
  TopBC,
  BottomBC
};

} // namespace solver
} // namespace alps
