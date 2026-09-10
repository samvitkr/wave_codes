#pragma once

#include "body_force_types.h"
#include "pressure_grad_options.h"
#include <common/container/vector_field.h>
#include <solvers/field/flow_field.h>

#include <fmt/format.h>

#include <array>

namespace alps {
namespace solver {

/// @brief Class describing a constant pressure gradient force
/** The constant pressure gradient force is given by
 * \f$\mathbf{F}_G = G_x \hat{x} + G_y \hat{y}\f$.
 *
 * The constant pressure gradient force is specified in the config file as a
 * table under BodyForce, e.g.
 * ```toml
 * [BodyForce]
 * PressureGradient = { type = "constant", values = [ 1.0, 0.0 ] }
 * ```
 * See @ref ChannelSolverOptions::parse_body_forces for config parsing.
 */
template<typename SolverType>
class ConstantPressureGradient : public BodyForce
{
 public:
  /** @brief Constructor with a given flow field and given gradients
   */
  ConstantPressureGradient(const SolverType&   flow_solver,
                           std::array<Real, 2> gradients)
    : solver_{flow_solver}
    , gradients_{gradients}
  {}

  ConstantPressureGradient(const SolverType&          flow_solver,
                           PressureGradOptions const& options)
    : solver_{flow_solver}
    , gradients_{(Real)options.gx, (Real)options.gy}
  {}

  void add_forces(const Vector3Field<Real***>&         fu,
                  const Kokkos::DefaultExecutionSpace& space) const override;

  std::string info() const override
  {
    return fmt::format(
      "Constant pressure gradient: G = ({}, {})", gradients_[0], gradients_[1]);
  }

 private:
  SolverType const&   solver_;
  std::array<Real, 2> gradients_;
};

} // namespace solver
} // namespace alps
