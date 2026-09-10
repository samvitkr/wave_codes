//
// Created by xuanx004 on 12/27/23.
//

#pragma once

#include <common/real_type.h>

#include <Kokkos_Core.hpp>

#include <string>

namespace alps::solver {

/// @brief Parent class for describing a scalar source term
/**
 * @note Child classes should implement the @ref add_source and @ref info
 * function
 */
class ScalarSource
{
 public:
  /// @brief Return a string describing the body force
  virtual std::string info() const = 0;

  virtual ~ScalarSource() = default;

  /// @brief Add the source term to the scalar equation
  /**
   * @param fc Right hand side of the scalar equation
   * @param space Execution space
   */
  virtual void add_source(const Kokkos::View<Real***, Kokkos::LayoutLeft>& fc,
                          const Kokkos::DefaultExecutionSpace& space) const = 0;

 protected:
  ScalarSource() = default;
};

} // namespace alps::solver
