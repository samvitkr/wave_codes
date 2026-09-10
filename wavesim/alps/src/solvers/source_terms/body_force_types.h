//
// Created by xuanx004 on 12/27/23.
//

#pragma once

#include <common/container/vector_field.h>
#include <common/real_type.h>

#include <Kokkos_Core.hpp>

#include <string>

namespace alps::solver {

/// @brief Parent class for describing a body force
/**
 * @note The body forces are added to a vector of BodyForce, which is then
 * invoked by the solver to add the forces to the momentum equation
 * @note Child classes should implement the @ref add_forces and @ref info
 * function; see @ref ConstantPressureGradient and @ref CoriolisForce for
 * examples.
 */
class BodyForce
{
 public:
  /// @brief Return a string describing the body force
  virtual std::string info() const = 0;

  virtual ~BodyForce() = default;

  /// @brief Add the body force to the momentum equation
  /**
   * @param fu Right hand side of the momentum equation
   * @param space Execution space
   */
  virtual void add_forces(const Vector3Field<Real***>&         fu,
                          const Kokkos::DefaultExecutionSpace& space) const = 0;

 protected:
  BodyForce() = default;
};

} // namespace alps::solver
