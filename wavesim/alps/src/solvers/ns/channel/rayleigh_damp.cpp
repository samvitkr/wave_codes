#include <solvers/source_terms/rayleigh_damp.h>

#include <common/container/view_utils.h>
#include <common/kokkos_abstraction/exec_policy.h>
#include <solvers/ns/channel/solver_ab2.h>
#include <solvers/ns/channel/solver_ab2cn.h>

#include <utility>

namespace alps::solver {

namespace {
template<typename FieldT>
void add_rayleigh_damping(Vector3Field<Real***> const&         fu,
                          FieldT const&                        field,
                          Real                                 z0,
                          Real                                 z1,
                          Real                                 rm,
                          Real                                 a_u,
                          Real                                 b_u,
                          Real                                 a_v,
                          Real                                 b_v,
                          Kokkos::DefaultExecutionSpace const& space)
{
  auto constexpr pi_v = Kokkos::numbers::pi_v<Real>;
  auto const& fux     = fu.x;
  auto const& fuy     = fu.y;
  auto const& fuz     = fu.z;
  auto const& u       = field.u.x;
  auto const& v       = field.u.y;
  auto const& w       = field.u.z;
  auto const& zz      = field.mesh.zz;
  auto const& zw      = field.mesh.zw;

  Kokkos::parallel_for(
    "add Rayleigh damping",
    LoopPolicy<3>(space, local_begins(fux), local_ends(fux)),
    KOKKOS_LAMBDA(int i, int j, int k) {
      // from z0 to z1, coeff varies from 0 to 1
      if ((z0 < z1 && zz(k) > z0) || (z0 >= z1 && zz(k) < z0)) {
        Real coeff = 1;
        if ((z0 < z1 && zz(k) < z1) || (z0 >= z1 && zz(k) > z1)) {
          coeff = (Kokkos::cos((zz(k) - z1) / (z1 - z0) * pi_v) + 1) / 2;
          coeff = Kokkos::clamp(coeff, (Real)0, (Real)1);
        }
        coeff *= rm;
        fux(i, j, k) -= coeff * (u(i, j, k) - b_u - a_u * zz(k));
        fuy(i, j, k) -= coeff * (v(i, j, k) - b_v - a_v * zz(k));
      }
      if ((z0 < z1 && zw(k) >= z0) || (z0 > z1 && zw(k) <= z0)) {
        Real coeff = 1;
        if ((z0 < z1 && zw(k) <= z1) || (z0 > z1 && zw(k) >= z1)) {
          coeff = (Kokkos::cos((zw(k) - z1) / (z1 - z0) * pi_v) + 1) / 2;
          coeff = Kokkos::clamp(coeff, (Real)0, (Real)1);
        }
        fuz(i, j, k) -= rm * coeff * w(i, j, k);
      }
    });
}
} // anonymous namespace

template<>
void RayleighDamp<ChannelFlowSolverAB2>::add_forces(
  const Vector3Field<Real***>&         fu,
  const Kokkos::DefaultExecutionSpace& space) const
{
  add_rayleigh_damping(fu,
                       this->solver_.flow_field,
                       z_zero,
                       z_one,
                       (Real)factor,
                       (Real)u_a,
                       (Real)u_b,
                       (Real)v_a,
                       (Real)v_b,
                       space);
}

template<>
void RayleighDamp<ChannelFlowSolverAB2CN>::add_forces(
  const Vector3Field<Real***>&         fu,
  const Kokkos::DefaultExecutionSpace& space) const
{
  add_rayleigh_damping(fu,
                       this->solver_.flow_field,
                       z_zero,
                       z_one,
                       (Real)factor,
                       (Real)u_a,
                       (Real)u_b,
                       (Real)v_a,
                       (Real)v_b,
                       space);
}

template<>
void RayleighDampScalar<ChannelFlowSolverAB2>::add_source(
  const Kokkos::View<Real***, Kokkos::LayoutLeft>& fc,
  const Kokkos::DefaultExecutionSpace&             space) const
{
  auto constexpr pi_v = Kokkos::numbers::pi_v<Real>;
  auto const& f       = solver_.flow_field.scalars.at(scalar_index).array;
  auto const& zz      = solver_.flow_field.mesh.zz;

  auto z0 = this->z_zero;
  auto z1 = this->z_one;
  auto rm = static_cast<Real>(this->factor);
  auto a  = static_cast<Real>(this->background_a);
  auto b  = static_cast<Real>(this->background_b);
  Kokkos::parallel_for(
    "add scalar Rayleigh damping",
    LoopPolicy<3>(space, local_begins(f), local_ends(f)),
    KOKKOS_LAMBDA(int i, int j, int k) {
      if ((z0 < z1 && zz(k) > z0) || (z0 >= z1 && zz(k) < z0)) {
        Real coeff = 1;
        if ((z0 < z1 && zz(k) < z1) || (z0 >= z1 && zz(k) > z1)) {
          coeff = (Kokkos::cos((zz(k) - z1) / (z1 - z0) * pi_v) + 1) / 2;
          coeff = Kokkos::clamp(coeff, (Real)0, (Real)1);
        }
        coeff *= rm;
        fc(i, j, k) -= coeff * (f(i, j, k) - b - a * zz(k));
      }
    });
}
} // namespace alps::solver
