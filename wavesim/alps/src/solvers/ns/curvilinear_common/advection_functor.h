//
// Created by xuananqing on 7/12/24.
//

#pragma once

#include <common/container/vector_field.h>
#include <common/container/view_utils.h>
#include <common/real_type.h>
#include <solvers/mesh/mesh.h>

#include <Kokkos_Core.hpp>

#include <type_traits>

namespace alps::solver {

template<typename mesh_t>
struct CurveAdvectionFunctor
{
  static_assert(std::is_base_of_v<CurvilinearMesh, mesh_t>);

  struct Bottom_0
  {}; // k == 0
  struct Top_NzMinus2
  {}; // k == nz - 2
  struct Top_NzMinus1
  {}; // k == nz - 1
  struct General
  {};

  template<typename T, typename S, typename U>
  KOKKOS_FORCEINLINE_FUNCTION static auto
  invJContravW(T u, T v, T w, S eta_x, S eta_y, U z)
  {
    return mesh_t::invJ_zeta_x(z, eta_x) * u + mesh_t::invJ_zeta_y(z, eta_y) * v
         + mesh_t::invJ_zeta_z() * w;
  }

  CurveAdvectionFunctor(const Vector3Field<Real***>&   fluxes,
                        const mesh_t&                  mesh,
                        const Vector3Field<Real***>&   vec_u,
                        const HaloView<Real const***>& scalar)
    : fx{fluxes.x}
    , fy{fluxes.y}
    , fz{fluxes.z}
    , u{vec_u.x}
    , v{vec_u.y}
    , w{vec_u.z}
    , f{scalar}
    , ex{mesh.ex}
    , ey{mesh.ey}
    , invJ{mesh.invJ}
    , et{mesh.et}
    , zeta{mesh.zz}
    , zeta_w{mesh.zw}
    , dzw{mesh.dzw}
    , nz{mesh.extent(2)}
  {
    if (!check_same_layout_and_offset(fx, fy, fz, u, v, w, f)) {
      throw std::runtime_error("Mismatched extents in vec_u, f and fluxes");
    }
  }

  using view_t       = HaloView<Real***>;
  using const_view_t = view_t::const_type;
  view_t       fx, fy, fz;
  const_view_t u, v, w;
  const_view_t f;

  MDView<Real const**>  ex, ey;
  MDView<Real const**>  invJ;
  MDView<Real const**>  et;
  HaloView<Real const*> zeta;
  HaloView<Real const*> zeta_w;
  HaloView<Real const*> dzw;
  int                   nz;

  KOKKOS_FUNCTION void operator()(const Bottom_0& /*tag*/, int i, int j) const
  {
    fx(i, j, 0) = -invJ(i, j) * u(i, j, 0) * f(i, j, 0);
    fy(i, j, 0) = -invJ(i, j) * v(i, j, 0) * f(i, j, 0);

    auto W =
      invJContravW(
        u(i, j, 0), v(i, j, 0), w(i, j, 0), ex(i, j), ey(i, j), zeta_w(0))
      - mesh_t::invJWg(zeta_w(0), et(i, j));
    fz(i, j, 0) = -W * f(i, j, 0);
  }

  KOKKOS_FUNCTION void
  operator()(const Top_NzMinus2& /*tag*/, int i, int j) const
  {
    const auto k = nz - 2;
    fx(i, j, k)  = -invJ(i, j) * u(i, j, k) * f(i, j, k);
    fy(i, j, k)  = -invJ(i, j) * v(i, j, k) * f(i, j, k);

    auto W = invJContravW(u(i, j, k + 1),
                          v(i, j, k + 1),
                          w(i, j, k),
                          ex(i, j),
                          ey(i, j),
                          zeta_w(k))
           - mesh_t::invJWg(zeta_w(k), et(i, j));
    fz(i, j, k) = -W * f(i, j, k + 1);
  }

  KOKKOS_FUNCTION void
  operator()(const Top_NzMinus1& /*tag*/, int i, int j) const
  {
    const auto k = nz - 1;
    fx(i, j, k)  = -invJ(i, j) * u(i, j, k) * f(i, j, k);
    fy(i, j, k)  = -invJ(i, j) * v(i, j, k) * f(i, j, k);

    auto W =
      invJContravW(
        u(i, j, k), v(i, j, k), w(i, j, k - 1), ex(i, j), ey(i, j), zeta(k))
      - mesh_t::invJWg(zeta(k), et(i, j));
    fz(i, j, k) = -W * f(i, j, k - 1);
  }

  KOKKOS_FUNCTION void
  operator()(const General& /*tag*/, int i, int j, int k) const
  {
    auto offset    = &f(i, j, k) - f.data();
    auto offset1   = offset + f.stride(2);
    auto offset_ij = &invJ(i, j) - invJ.data();

    fx(i, j, k) = -invJ.data()[offset_ij] * u.data()[offset] * f.data()[offset];
    fy(i, j, k) = -invJ.data()[offset_ij] * v.data()[offset] * f.data()[offset];

    auto beta = dzw(k - 1) / (dzw(k - 1) + dzw(k));
    auto u_w  = itp2node(u.data()[offset], u.data()[offset1], beta);
    auto v_w  = itp2node(v.data()[offset], v.data()[offset1], beta);
    auto f_w  = itp2node(f.data()[offset], f.data()[offset1], beta);

    auto W = invJContravW(u_w,
                          v_w,
                          w.data()[offset],
                          ex.data()[offset_ij],
                          ey.data()[offset_ij],
                          zeta_w(k))
           - mesh_t::invJWg(zeta_w(k), et.data()[offset_ij]);
    fz(i, j, k) = -W * f_w;
  }
};

} // namespace alps::solver
