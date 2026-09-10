#pragma once

#include <common/container/vector_field.h>
#include <common/container/view_utils.h>
#include <common/kokkos_abstraction/exec_policy.h>
#include <common/real_type.h>
#include <solvers/mesh/curvilinear_mesh.h>
#include <solvers/source_terms/rayleigh_damp.h>

#include <type_traits>

namespace alps::solver {
template<typename SolverType>
void RayleighDamp<SolverType>::add_forces(
  const Vector3Field<Real***>&         fu,
  const Kokkos::DefaultExecutionSpace& space) const
{
  auto const& mesh = solver_.flow_field.mesh;
  static_assert(
    std::is_base_of_v<CurvilinearMesh, std::decay_t<decltype(mesh)>>);

  auto constexpr pi_v = Kokkos::numbers::pi_v<Real>;
  auto const& fux     = fu.x;
  auto const& fuy     = fu.y;
  auto const& fuz     = fu.z;
  auto const& u       = solver_.flow_field.u.x;
  auto const& v       = solver_.flow_field.u.y;
  auto const& w       = solver_.flow_field.u.z;
  auto const& zz      = mesh.zz;
  auto const& zw      = mesh.zw;
  auto const& invJ    = mesh.invJ;

  auto z0  = this->z_zero;
  auto z1  = this->z_one;
  auto rm  = static_cast<Real>(this->factor);
  auto a_u = static_cast<Real>(this->u_a);
  auto b_u = static_cast<Real>(this->u_b);
  auto a_v = static_cast<Real>(this->v_a);
  auto b_v = static_cast<Real>(this->v_b);
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
        fux(i, j, k) -= coeff * (u(i, j, k) - b_u - a_u * zz(k)) * invJ(i, j);
        fuy(i, j, k) -= coeff * (v(i, j, k) - b_v - a_v * zz(k)) * invJ(i, j);
      }
      if ((z0 < z1 && zw(k) >= z0) || (z0 > z1 && zw(k) <= z0)) {
        Real coeff = 1;
        if ((z0 < z1 && zw(k) <= z1) || (z0 > z1 && zw(k) >= z1)) {
          coeff = (Kokkos::cos((zw(k) - z1) / (z1 - z0) * pi_v) + 1) / 2;
          coeff = Kokkos::clamp(coeff, (Real)0, (Real)1);
        }
        fuz(i, j, k) -= rm * coeff * w(i, j, k) * invJ(i, j);
      }
    });
}

template<typename SolverType>
void RayleighDampScalar<SolverType>::add_source(
  const Kokkos::View<Real***, Kokkos::LayoutLeft>& fc,
  const Kokkos::DefaultExecutionSpace&             space) const
{
  auto const& mesh = solver_.flow_field.mesh;
  static_assert(
    std::is_base_of_v<CurvilinearMesh, std::decay_t<decltype(mesh)>>);

  auto constexpr pi_v = Kokkos::numbers::pi_v<Real>;
  auto const& f       = solver_.flow_field.scalars.at(scalar_index).array;
  auto const& zz      = mesh.zz;
  auto const& invJ    = mesh.invJ;

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
        fc(i, j, k) -= coeff * (f(i, j, k) - b - a * zz(k)) * invJ(i, j);
      }
    });
}
} // namespace alps::solver
