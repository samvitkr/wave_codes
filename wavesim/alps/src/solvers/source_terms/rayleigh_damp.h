#pragma once

#include "body_force_types.h"
#include "rayleigh_damp_options.h"
#include "scalar_source_types.h"
#include <common/container/vector_field.h>
#include <solvers/field/flow_field.h>

#include <fmt/format.h>

#include <array>

namespace alps {
namespace solver {

/**
 * @brief Rayleigh damping body force
 * Adds a body force to the momentum equation to damp the flow:
 *
 * \f[
 * f_x = -r(\zeta) * (u - u_0), \quad f_y = -r(\zeta) * (v - v_0)
 * \f]
 *
 * where \f$r(\zeta)\f$ varies from 0 to 1 between z0 and z1. \f$u0\f$ and
 *\f$v0\f$ are background states, which can be specified as linear functions of
 *\f$zeta\$:
 *
 *\f[
 * u_0 = u_b + u_a * \zeta, \quad v0 = v_b + v_a * \zeta
 *\f]
 *
 * The damping parameter ss specified in the config file as a table under
 * ```toml
 * [BodyForce]
 * RayleighDamping = { factor = 1.0, z_zero = 0.1, z_one = 0.0, u_b = 1.0, u_a =
 *0.0, v_b = -1.0, v_a = 0.0 }
 * ```
 */
template<typename SolverType>
class RayleighDamp : public BodyForce
{
 public:
  /** @brief Constructor with a given flow field and given gradients
   */
  RayleighDamp(const SolverType& flow_solver,
               double            factor_,
               double            z_zero_,
               double            z_one_,
               double            u_b_,
               double            u_a_,
               double            v_b_,
               double            v_a_)
    : solver_{flow_solver}
    , factor{factor_}
    , z_zero(static_cast<Real>(z_zero_))
    , z_one(static_cast<Real>(z_one_))
    , u_b(static_cast<Real>(u_b_))
    , u_a(static_cast<Real>(u_a_))
    , v_b(static_cast<Real>(v_b_))
    , v_a(static_cast<Real>(v_a_))
  {}

  RayleighDamp(const SolverType&          flow_solver,
               RayleighDampOptions const& options)
    : RayleighDamp(flow_solver,
                   options.factor,
                   options.z_zero,
                   options.z_one,
                   options.u_b,
                   options.u_a,
                   options.v_b,
                   options.v_a)
  {}

  void add_forces(const Vector3Field<Real***>&         fu,
                  const Kokkos::DefaultExecutionSpace& space) const override;

  std::string info() const override
  {
    return fmt::format("Rayleigh damping: factor = {}, zeta range = ({}, {}), "
                       "background: u = {}*z{:+}, v = {}*z{:+}",
                       factor,
                       z_zero,
                       z_one,
                       u_a,
                       u_b,
                       v_a,
                       v_b);
  }

 private:
  SolverType const& solver_;
  double            factor; /// damping strength is multiplied by this factor
  Real              z_zero; /// zeta coordinate where damping starts
  Real              z_one; /// zeta coordinate where damping reaches its maximum
  Real              u_b;   /// background state u = b + a * z
  Real              u_a;   /// background state u = b + a * z
  Real              v_b;   /// background state v = b + a * z
  Real              v_a;   /// background state v = b + a * z
};

/**
 * @brief Rayleigh damping scalar source term
 * Adds a source term to the scalar equation to damp the scalar field:
 *
 * \f[
 * f = -r(\zeta) * (c - c_0)
 * \f]
 *
 * where \f$r(\zeta)\f$ varies from 0 to 1 between z0 and z1. \f$c0\f$ is a
 * background state, which can be specified as a linear function of \f$zeta\$:
 * \f[
 * c_0 = c_b + c_a * \zeta
 * \f]
 */
template<typename SolverType>
class RayleighDampScalar : public ScalarSource
{
 public:
  /** @brief Constructor with a given flow field and given gradients
   */
  RayleighDampScalar(const SolverType& flow_solver,
                     int               scalar_index_,
                     double            factor_,
                     double            z_zero_,
                     double            z_one_,
                     double            background_b_ = 0.0,
                     double            background_a_ = 0.0)
    : solver_{flow_solver}
    , scalar_index(scalar_index_)
    , factor{factor_}
    , z_zero(static_cast<Real>(z_zero_))
    , z_one(static_cast<Real>(z_one_))
    , background_b(static_cast<Real>(background_b_))
    , background_a(static_cast<Real>(background_a_))
  {
    if (scalar_index < 0
        || scalar_index >= (int)solver_.flow_field.scalars.size()) {
      throw std::runtime_error(
        fmt::format("Invalid scalar index: {}", scalar_index));
    }
  }

  RayleighDampScalar(const SolverType&                flow_solver,
                     int                              scalar_index_,
                     ScalarRayleighDampOptions const& options)
    : RayleighDampScalar(flow_solver,
                         scalar_index_,
                         options.factor,
                         options.z_zero,
                         options.z_one,
                         options.b,
                         options.a)
  {}

  void add_source(const Kokkos::View<Real***, Kokkos::LayoutLeft>& fc,
                  const Kokkos::DefaultExecutionSpace& space) const override;

  std::string info() const override
  {
    return fmt::format("Rayleigh damping: factor = {}, zeta range = ({}, {}), "
                       "background c_b = {}*z{:+}",
                       factor,
                       z_zero,
                       z_one,
                       background_a,
                       background_b);
  }

 private:
  SolverType const& solver_;
  int               scalar_index;
  double            factor; /// damping strength is multiplied by this factor
  Real              z_zero; /// zeta coordinate where damping starts
  Real              z_one; /// zeta coordinate where damping reaches its maximum
  Real              background_b; /// background state y = b + a * z
  Real              background_a; /// background state y = b + a * z
};
} // namespace solver
} // namespace alps
