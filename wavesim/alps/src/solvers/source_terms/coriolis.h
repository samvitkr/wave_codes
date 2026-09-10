//
// Created by xuanx004 on 12/30/23.
//

#pragma once

#include "body_force_types.h"
#include "coriolis_options.h"

#include <fmt/format.h>

namespace alps::solver {

/// @brief Class describing a Coriolis force
/** The Coriolis force is given by \f$\mathbf{F}_C = f_z\hat{\mathbf{z}} \times
 * \mathbf{u}\f$, where \f$f_z\f$ is the Coriolis parameter.
 *
 * The Coriolis force is specified in the config file as a table under
 * BodyForce, e.g.
 * ```toml
 * [BodyForce]
 * Coriolis = { f = 1e-5 }
 * ```
 * See @ref CoriolisOptions::parse for config parsing.
 */
template<typename SolverType>
class CoriolisForce : public BodyForce
{
 public:
  /** @brief Constructor with a given flow field and Coriolis parameter f_z
   */
  CoriolisForce(const SolverType& flow_solver, double f_z)
    : solver{flow_solver}
    , fz{(Real)f_z}
  {}

  CoriolisForce(const SolverType& flow_solver, CoriolisOptions const& options)
    : solver{flow_solver}
    , fz{(Real)options.fz}
  {}

  void add_forces(const Vector3Field<Real***>&         fu,
                  const Kokkos::DefaultExecutionSpace& space) const override;

  std::string info() const override
  {
    return fmt::format("Coriolis force: f = {}", fz);
  }

 private:
  SolverType const& solver;
  Real              fz;
};

} // namespace alps::solver
