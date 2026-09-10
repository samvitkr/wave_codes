//
// Created by xuanx004 on 7/22/24.
//

#pragma once

#include "body_force_types.h"
#include "boussinesq_options.h"
#include <solvers/field/scalar_fields.h>

#include <fmt/format.h>

namespace alps::solver {

class ScalarField;
class Mesh;

/// @brief Class describing a Boussinesq buoyancy force
/** The buoyancy force is given by \f$\mathbf{F} = -Ri (theta - theta_0)
 * \hat{e}_z\f$.
 *
 * ```toml
 * [BodyForce]
 * Boussinesq = { scalar_id = 1, Ri = 980.0, theta0 = 300.0 }
 * ```
 * The scalar theta can specified in the config file either by scalar_id or
 * scalar_label. If theta0 is omitted or cannot be parsed as a number, the
 * reference value theta_0 would be the plane average of the scalar.
 */
template<typename SolverType>
class BoussinesqForce : public BodyForce
{
 public:
  BoussinesqForce(const SolverType&        flow_solver,
                  BoussinesqOptions const& options)
    : solver{flow_solver}
    , scalar{[&scalars = flow_solver.flow_field.scalars, &options] {
      if (!options.scalar_label.empty()) {
        for (auto const& scalar_ : scalars) {
          if (scalar_.label() == options.scalar_label) {
            return &scalar_;
          }
        }
      }
      if (options.scalar_id >= 0) {
        if (std::size_t(options.scalar_id) >= scalars.size()) {
          throw std::runtime_error(fmt::format(
            "Scalar field id specified in Boussinesq force (id={}) is out of "
            "bounds (max id: {})",
            options.scalar_id,
            scalars.size() - 1));
        }
      } else {
        throw std::runtime_error("Matching scalar field id or label for "
                                 "Boussinesq force cannot be found");
      }
      return &scalars.at(options.scalar_id);
    }()}
    , Ri{(Real)options.Ri}
    , theta0{options.ref_scalar}
  {}

  void add_forces(const Vector3Field<Real***>&         fu,
                  const Kokkos::DefaultExecutionSpace& space) const override;

  std::string info() const override
  {
    if (theta0 == std::numeric_limits<double>::max()) {
      return fmt::format("Boussinesq buoyancy force: beta = {}, theta = "
                         "\"{}\", theta0 = <average>",
                         Ri,
                         scalar->label());
    }
    return fmt::format(
      "Boussinesq buoyancy force: beta = {}, theta = \"{}\", theta0 = {}",
      Ri,
      scalar->label(),
      theta0);
  }

  auto const* get_scalar() const { return scalar; }

  auto get_Ri() const { return Ri; }

 private:
  SolverType const&  solver;
  ScalarField const* scalar;
  Real               Ri;
  double             theta0{std::numeric_limits<double>::max()};
};

void add_buoyancy_force_with_ref_scalar(
  HaloView<Real***> const&             fz,
  HaloView<Real const***> const&       theta,
  Real                                 theta0,
  Real                                 beta,
  Mesh const&                          mesh,
  Kokkos::DefaultExecutionSpace const& stream);

void add_buoyancy_force_without_ref_scalar(
  HaloView<Real***> const&             fz,
  HaloView<Real const***> const&       theta,
  Real                                 beta,
  Mesh const&                          mesh,
  Kokkos::DefaultExecutionSpace const& stream);

} // namespace alps::solver
