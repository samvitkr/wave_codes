#pragma once

#include <common/container/matrix_field.h>
#include <common/container/view_utils.h>
#include <solvers/mesh/curvilinear_mesh.h>

#include <Kokkos_Core.hpp>

namespace alps::solver {

template<typename mesh_t>
struct CurveConvectionFunctor
{
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

  CurveConvectionFunctor(const Tensor33Field<Real***>& fluxes,
                         const mesh_t&                 mesh,
                         const Vector3Field<Real***>&  vec_u)
    : fxx{fluxes.xx}
    , fxy{fluxes.xy}
    , fxz{fluxes.xz}
    , fyx{fluxes.yx}
    , fyy{fluxes.yy}
    , fyz{fluxes.yz}
    , fzx{fluxes.zx}
    , fzy{fluxes.zy}
    , fzz{fluxes.zz}
    , u{vec_u.x}
    , v{vec_u.y}
    , w{vec_u.z}
    , ex{mesh.ex}
    , ey{mesh.ey}
    , invJ{mesh.invJ}
    , et{mesh.et}
    , zeta{mesh.zz}
    , zeta_w{mesh.zw}
    , dzw{mesh.dzw}
    , nz{local_extent(this->u, 2)}
  {
    if (!check_same_layout_and_offset(
          fxx, fxy, fxz, fyx, fyy, fyz, fzx, fzy, fzz, u, v, w)) {
      throw std::runtime_error("Mismatched extents in vec_u and fluxes");
    }
  }

  using view_t       = HaloView<Real***>;
  using const_view_t = view_t::const_type;
  view_t       fxx, fxy, fxz;
  view_t       fyx, fyy, fyz;
  view_t       fzx, fzy, fzz;
  const_view_t u, v, w;

  MDView<Real const**>  ex, ey;
  MDView<Real const**>  invJ;
  MDView<Real const**>  et;
  HaloView<Real const*> zeta;
  HaloView<Real const*> zeta_w;
  HaloView<Real const*> dzw;
  int                   nz;

  KOKKOS_FUNCTION void operator()(const Bottom_0& /*tag*/, int i, int j) const
  {
    fxx(i, j, 0) = -invJ(i, j) * u(i, j, 0) * u(i, j, 0);
    auto r       = -invJ(i, j) * u(i, j, 0) * v(i, j, 0);
    fxy(i, j, 0) = r;
    fyx(i, j, 0) = r;
    fyy(i, j, 0) = -invJ(i, j) * v(i, j, 0) * v(i, j, 0);

    auto W =
      invJContravW(
        u(i, j, 0), v(i, j, 0), w(i, j, 0), ex(i, j), ey(i, j), zeta_w(0))
      - mesh_t::invJWg(zeta_w(0), et(i, j));
    fxz(i, j, 0) = -W * u(i, j, 0);
    fzx(i, j, 0) = -invJ(i, j) * u(i, j, 0) * w(i, j, 0);
    fyz(i, j, 0) = -W * v(i, j, 0);
    fzy(i, j, 0) = -invJ(i, j) * v(i, j, 0) * w(i, j, 0);
    fzz(i, j, 0) = -W * w(i, j, 0);
  }

  KOKKOS_FUNCTION void
  operator()(const Top_NzMinus2& /*tag*/, int i, int j) const
  {
    const auto k = nz - 2;
    fxx(i, j, k) = -invJ(i, j) * u(i, j, k) * u(i, j, k);
    auto r       = -invJ(i, j) * u(i, j, k) * v(i, j, k);
    fxy(i, j, k) = r;
    fyx(i, j, k) = r;
    fyy(i, j, k) = -invJ(i, j) * v(i, j, k) * v(i, j, k);

    auto W = invJContravW(u(i, j, k + 1),
                          v(i, j, k + 1),
                          w(i, j, k),
                          ex(i, j),
                          ey(i, j),
                          zeta_w(k))
           - mesh_t::invJWg(zeta_w(k), et(i, j));
    fxz(i, j, k) = -W * u(i, j, k + 1);
    fzx(i, j, k) = -invJ(i, j) * u(i, j, k + 1) * w(i, j, k);
    fyz(i, j, k) = -W * v(i, j, k + 1);
    fzy(i, j, k) = -invJ(i, j) * v(i, j, k + 1) * w(i, j, k);

    r = (w(i, j, k) + w(i, j, k - 1)) / 2;
    W = invJContravW(u(i, j, k), v(i, j, k), r, ex(i, j), ey(i, j), zeta(k))
      - mesh_t::invJWg(zeta(k), et(i, j));
    fzz(i, j, k) = -W * r;
  }

  KOKKOS_FUNCTION void
  operator()(const Top_NzMinus1& /*tag*/, int i, int j) const
  {
    const auto k = nz - 1;
    fxx(i, j, k) = -invJ(i, j) * u(i, j, k) * u(i, j, k);
    auto r       = -invJ(i, j) * u(i, j, k) * v(i, j, k);
    fxy(i, j, k) = r;
    fyx(i, j, k) = r;
    fyy(i, j, k) = -invJ(i, j) * v(i, j, k) * v(i, j, k);

    auto W =
      invJContravW(
        u(i, j, k), v(i, j, k), w(i, j, k - 1), ex(i, j), ey(i, j), zeta(k))
      - mesh_t::invJWg(zeta(k), et(i, j));
    fzz(i, j, k) = -W * w(i, j, k - 1);
  }

  KOKKOS_FUNCTION void
  operator()(const General& /*tag*/, int i, int j, int k) const
  {
    auto offset = &u(i, j, k) - u.data(); // (i, j, k)
    auto uijk   = u.data()[offset];
    auto vijk   = v.data()[offset];
    auto wijk   = w.data()[offset];

    fxx.data()[offset] = -invJ(i, j) * uijk * uijk;
    auto r             = -invJ(i, j) * uijk * vijk;
    fxy.data()[offset] = r;
    fyx.data()[offset] = r;
    fyy.data()[offset] = -invJ(i, j) * vijk * vijk;

    auto beta = dzw(k - 1) / (dzw(k - 1) + dzw(k));
    auto u_w  = itp2node(uijk, u.data()[offset + u.stride(2)], beta);
    auto v_w  = itp2node(vijk, v.data()[offset + v.stride(2)], beta);

    auto W = invJContravW(u_w, v_w, wijk, ex(i, j), ey(i, j), zeta_w(k))
           - mesh_t::invJWg(zeta_w(k), et(i, j));
    fxz.data()[offset] = -W * u_w;
    fyz.data()[offset] = -W * v_w;
    fzx.data()[offset] = -invJ(i, j) * u_w * wijk;
    fzy.data()[offset] = -invJ(i, j) * v_w * wijk;

    r = itp2center(wijk, w.data()[offset - w.stride(2)]);
    W = invJContravW(uijk, vijk, r, ex(i, j), ey(i, j), zeta(k))
      - mesh_t::invJWg(zeta(k), et(i, j));
    fzz.data()[offset] = -W * r;
  }
};

} // namespace alps::solver
