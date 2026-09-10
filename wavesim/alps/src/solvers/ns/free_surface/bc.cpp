//
// Created by xuananqing on 7/4/23.
//

#include "bc.h"

#include "solver.h"
#include <common/base/logging.h>
#include <common/base/macros.h>
#include <common/container/view_utils.h>
#include <common/kokkos_abstraction/exec_policy.h>
#include <common/kokkos_abstraction/pool_space.h>
#include <common/math.h>
#include <common/runtime/async_utils.h>
#include <solvers/mesh/curvilinear_mesh.h>
#include <solvers/ns/curvilinear_common/apply_bc.h>
#include <solvers/operators/grad.h>
#include <spectral/spectral.h>

#include <mpipp/collectives.h>

namespace alps::solver {

std::string FreeSurfaceBC::info() const
{
  return fmt::format("FreeSurfaceBC, hyperviscosity (n={}, c={})",
                     options.hyper_viscosity_n,
                     options.hyper_viscosity_c);
}

FreeSurfaceBC::FreeSurfaceBC(TopWaveMesh const&          mesh,
                             FreeSurfaceBCOptions const& options_)
  : VelocityBC{false}
  , options(options_)
  , tau_1(MDView<Real**, default_memory_pool>("tau_1",
                                              mesh.extent(0),
                                              mesh.extent(1)))
  , tau_2(MDView<Real**, default_memory_pool>("tau_2",
                                              mesh.extent(0),
                                              mesh.extent(1)))
  , ue_(MDView<Real** [3], default_memory_pool>("ue",
                                                mesh.extent(0),
                                                mesh.extent(1)))
  , ve_(MDView<Real** [3], default_memory_pool>("ve",
                                                mesh.extent(0),
                                                mesh.extent(1)))
  , we_(MDView<Real** [2], default_memory_pool>("we",
                                                mesh.extent(0),
                                                mesh.extent(1)))
  , p_top_0(MDView<Real**, default_memory_pool>(
      "p_top_0", // only allocated at the top boundary
      mesh.comm().is_last(2) ? mesh.extent(0) : 1,
      mesh.comm().is_last(2) ? mesh.extent(1) : 1))
{}

void FreeSurfaceBC::update_eta_t(TopWaveMesh const&          mesh,
                                 MDView<Real const**> const& u_s,
                                 MDView<Real const**> const& v_s,
                                 MDView<Real const**> const& w_s) const
{
  using Kokkos::parallel_for;

  auto const  is_top = mesh.comm().is_last(2);
  auto const& et     = mesh.et;
  auto const& ex     = mesh.ex;
  auto const& ey     = mesh.ey;

  if (is_top) {
    auto stream = get_next_stream();
    parallel_for(
      "calc eta_t",
      LoopPolicy<2>(stream, {0, 0}, {et.extent(0), et.extent(1)}),
      KOKKOS_LAMBDA(int i, int j) {
        et(i, j) = w_s(i, j) - (u_s(i, j) * ex(i, j) + v_s(i, j) * ey(i, j));
      });
    stream.fence();
  }
  mpipp::bcast(nonstd::span(et.data(), et.size()),
               mesh.comm().dims[2] - 1,
               mesh.comm().axis_comm[2]);
}

namespace {
template<int N, typename T>
KOKKOS_FORCEINLINE_FUNCTION constexpr auto power(T v)
{
  static_assert(N >= 1, "power<N> requires N >= 1");
  if constexpr (N == 1) {
    return v;
  } else {
    return v * power<N - 1>(v);
  }
  ALPS_UNREACHABLE(T(1));
}
} // namespace

template<int E>
void apply_hyper_viscosity_eta(MDView<Real**> const& eta,
                               Real                  nu_hyper,
                               std::integral_constant<int, E> /*nu_power*/,
                               Real                                 dt,
                               Grid const&                          grid,
                               Kokkos::DefaultExecutionSpace const& stream)
{
  auto const r2c_layout =
    grid.get_r2c_xy_output_layout<Real, Kokkos::DefaultExecutionSpace>();

  MDView<Real** [1], default_memory_pool> const eta_r2c(
    Kokkos::view_alloc("eta_r2c", Kokkos::WithoutInitializing),
    r2c_layout.dimension[0],
    r2c_layout.dimension[1]);
  MDView<Real** [1]> eta1(eta.data(), eta.extent(0), eta.extent(1));

  spectral::fft_r2c_xy(eta_r2c, eta1, grid, stream);

  auto const kx0      = grid.pex;
  auto const ky0      = grid.pey;
  auto const offset_x = grid.offset(1, Pencil::Y);

  Kokkos::parallel_for(
    "eta hyper viscosity",
    LoopPolicy<2>(
      stream, {0, 0}, {grid.extent(0, Pencil::Y), grid.extent(1, Pencil::Y)}),
    KOKKOS_LAMBDA(int l, int m) {
      auto ax = Real(kx0 * int((m + offset_x) / 2));
      auto ay = Real(ky0 * int(l / 2));

      auto c = 1 + (power<E>(ax) + power<E>(ay)) * nu_hyper * dt;
      eta_r2c(l, m, 0) /= c;
    });

  spectral::fft_c2r_xy(eta1, eta_r2c, false, grid, stream);

  stream.fence(); // fence to ensure deallocation of eta_r2c
}

void FreeSurfaceBC::get_updated_eta(MDView<Real**> const&       eta,
                                    MDView<Real const**> const& u_s,
                                    MDView<Real const**> const& v_s,
                                    MDView<Real const**> const& w_s,
                                    double                      nu,
                                    Real                        dt,
                                    int                         rk_stage,
                                    TopWaveMesh const&          mesh) const
{
  using Kokkos::parallel_for;

  auto const  is_top = mesh.comm().is_last(2);
  auto const& et     = mesh.et;
  auto const& ex     = mesh.ex;
  auto const& ey     = mesh.ey;
  auto const  nx     = mesh.extent(0);
  auto const  ny     = mesh.extent(1);

  auto stream1 = get_next_stream();

  if (rk_stage == 1) {
    parallel_for(
      "calc eta",
      LoopPolicy<2>(stream1, {0, 0}, {nx, ny}),
      KOKKOS_LAMBDA(int i, int j) { eta(i, j) += dt * et(i, j); });
    dealias(MDView<Real** [1]>(eta.data(), nx, ny), mesh.grid, stream1);
    stream1.fence();
  }
  if (rk_stage == 2) {
    if (is_top) {
      auto const& et_old = mesh.et;
      parallel_for(
        "calc eta",
        LoopPolicy<2>(stream1, {0, 0}, {nx, ny}),
        KOKKOS_LAMBDA(int i, int j) {
          auto et_  = w_s(i, j) - (u_s(i, j) * ex(i, j) + v_s(i, j) * ey(i, j));
          eta(i, j) = Kokkos::fma(dt / 2, et_ - et_old(i, j), eta(i, j));
        });
      dealias(MDView<Real** [1]>(eta.data(), nx, ny), mesh.grid, stream1);
      stream1.fence();
    }
    mpipp::bcast(nonstd::span(eta.data(), eta.size()),
                 mesh.comm().dims[2] - 1,
                 mesh.comm().axis_comm[2]);
  }

  if (auto hyper_nu = Real(nu / options.hyper_viscosity_c);
      hyper_nu >= std::numeric_limits<Real>::min()) {
    if (options.hyper_viscosity_n == 1) {
      apply_hyper_viscosity_eta(eta,
                                hyper_nu,
                                std::integral_constant<int, 1>{},
                                dt,
                                mesh.grid,
                                stream1);
    }
    if (options.hyper_viscosity_n == 2) {
      apply_hyper_viscosity_eta(eta,
                                hyper_nu,
                                std::integral_constant<int, 2>{},
                                dt,
                                mesh.grid,
                                stream1);
    }
    if (options.hyper_viscosity_n == 3) {
      apply_hyper_viscosity_eta(eta,
                                hyper_nu,
                                std::integral_constant<int, 3>{},
                                dt,
                                mesh.grid,
                                stream1);
    }
    if (options.hyper_viscosity_n == 4) {
      apply_hyper_viscosity_eta(eta,
                                hyper_nu,
                                std::integral_constant<int, 4>{},
                                dt,
                                mesh.grid,
                                stream1);
    }
  }
}

// vec_u is (u)
// vec_invJuz is (d(J^{-1} u)/dzeta)
void FreeSurfaceBC::get_surface_uz(MDView<Real** [3]> const&    vec_invJuz,
                                   Vector3Field<Real***> const& vec_u,
                                   Real                         nu,
                                   const TopWaveMesh&           mesh) const
{
  using Kokkos::ALL;
  using Kokkos::parallel_for;
  using boundary_functor_t =
    FDBoundaryFunctor<Real, Kokkos::DefaultExecutionSpace>;

  auto const& grid   = mesh.grid;
  auto const  nz     = grid.extent(2);
  auto const  is_top = mesh.comm().is_last(2);

  if (!is_top) return; // only applies to the top boundary

  auto const us = subview(vec_u.x, ALL, ALL, nz - 1).view();
  auto const vs = subview(vec_u.y, ALL, ALL, nz - 1).view();
  auto const ws = subview(vec_u.z, ALL, ALL, index_range(nz - 4, nz - 1))
                    .view(); // shape (nx, ny, 3)

  MDView<Real** [3], default_memory_pool> const u_x(
    Kokkos::view_alloc("ux", Kokkos::WithoutInitializing),
    us.extent(0),
    us.extent(1));
  MDView<Real** [3], default_memory_pool> const u_y(
    Kokkos::view_alloc("uy", Kokkos::WithoutInitializing),
    us.extent(0),
    us.extent(1));
  MDView<Real** [2], default_memory_pool> const tmp2d(
    Kokkos::view_alloc("tmp", Kokkos::WithoutInitializing),
    us.extent(0),
    us.extent(1));

  auto                stream1 = get_next_stream();
  auto                stream2 = get_next_stream();
  LoopPolicy<2> const policy1(stream1, {0, 0}, {us.extent(0), us.extent(1)});
  LoopPolicy<2> const policy2(stream2, {0, 0}, {us.extent(0), us.extent(1)});

  // compute d(u,v,w)/dx and d(u,v,w)/dy
  for (auto const& [stream, dst] :
       {std::pair{stream1, u_x}, std::pair{stream2, u_y}}) {
    Kokkos::deep_copy(stream, subview(dst, ALL, ALL, 0), us);
    Kokkos::deep_copy(stream, subview(dst, ALL, ALL, 1), vs);
    Kokkos::deep_copy(
      stream, subview(dst, ALL, ALL, 2), subview(ws, ALL, ALL, end(ws, 2) - 1));
  }

  spectral::ddx(u_x, u_x, grid, stream1);
  spectral::ddy(u_y, u_y, grid, stream2);

  auto const& ex = mesh.ex;
  auto const& ey = mesh.ey;
  parallel_for(
    "metric coeff", policy1, KOKKOS_LAMBDA(int i, int j) {
      auto reh       = alps::rnorm3d((Real)1, ex(i, j), ey(i, j));
      tmp2d(i, j, 0) = Kokkos::hypot((Real)1, ex(i, j)) * reh;
      tmp2d(i, j, 1) = Kokkos::hypot((Real)1, ey(i, j)) * reh;
    });
  stream2.fence(); // subsequent dealias cannot execute concurrently with ddy
  spectral::dealias(tmp2d, grid, stream1);

  auto const& sigma1 = this->tau_1;
  auto const& sigma2 = this->tau_2;
  parallel_for(
    "scaled gradient", policy1, KOKKOS_LAMBDA(int i, int j) {
      auto ex2 = alps::square(ex(i, j));
      auto ey2 = alps::square(ey(i, j));
      auto uz  = (ex2 - 1) * u_x(i, j, 2) + ex(i, j) * ey(i, j) * u_y(i, j, 2)
              + 2 * ex(i, j) * u_x(i, j, 0) + ey(i, j) * u_y(i, j, 0)
              + ey(i, j) * u_x(i, j, 1);
      auto vz = (ey2 - 1) * u_y(i, j, 2) + ex(i, j) * ey(i, j) * u_x(i, j, 2)
              + 2 * ey(i, j) * u_y(i, j, 1) + ex(i, j) * u_x(i, j, 1)
              + ex(i, j) * u_y(i, j, 0);
      auto reh2      = 1 + ex2 + ey2;
      tmp2d(i, j, 0) = tmp2d(i, j, 0) * sigma1(i, j) / nu + uz / reh2;
      tmp2d(i, j, 1) = tmp2d(i, j, 1) * sigma2(i, j) / nu + vz / reh2;
    });
  spectral::dealias(tmp2d, grid, stream1);

  // dw/dzeta at the surface
  boundary_functor_t const functor(
    subview(vec_invJuz, ALL, ALL, 2),
    ws,
    subview(mesh.dzw_h, index_range(nz - 4, nz - 2)).view(),
    1,
    RightBoundary());
  parallel_for("dw/dzeta top", policy2, functor);

  auto const& invJ = mesh.invJ;
  policy2.space().fence();
  parallel_for(
    "surface uz", policy1, KOKKOS_LAMBDA(int i, int j) {
      auto wz = vec_invJuz(i, j, 2);
      vec_invJuz(i, j, 0) =
        (tmp2d(i, j, 0) * invJ(i, j) - ex(i, j) * wz) * invJ(i, j);
      vec_invJuz(i, j, 1) =
        (tmp2d(i, j, 1) * invJ(i, j) - ey(i, j) * wz) * invJ(i, j);
      vec_invJuz(i, j, 2) = wz * invJ(i, j);
    });
  spectral::dealias(vec_invJuz, grid, policy1.space());
  policy1.space().fence();
}

namespace {
template<class T, class S>
KOKKOS_FORCEINLINE_FUNCTION auto zetax(T zeta, S exr)
{
  return TopWaveMesh::zeta_x(zeta, exr);
}

template<class T, class S>
KOKKOS_FORCEINLINE_FUNCTION auto zetay(T zeta, S eyr)
{
  return TopWaveMesh::zeta_y(zeta, eyr);
}
} // anonymous namespace

void FreeSurfaceBC::get_surface_uhat(Vector3Field<Real***> const&    vec_invJu,
                                     MDView<Real const** [3]> const& vec_invJuz,
                                     HaloView<Real const***> const&  pp,
                                     Real                            dt,
                                     TopWaveMesh const&              mesh) const
{
  using Kokkos::parallel_for;

  auto const is_top = mesh.comm().is_last(2);
  if (!is_top) return;

  auto const  nx  = mesh.extent(0);
  auto const  ny  = mesh.extent(1);
  auto const  nz  = mesh.extent(2);
  auto const& dz  = mesh.dz;
  auto const& dzw = mesh.dzw;
  auto const& exr = mesh.exr;
  auto const& eyr = mesh.eyr;
  auto const& zz  = mesh.zz;
  auto const& zw  = mesh.zw;

  HaloView<Real***, default_memory_pool> const invJp(
    Kokkos::view_alloc("J^{-1}p at top", Kokkos::WithoutInitializing),
    {0, nx - 1},
    {0, ny - 1},
    {nz - 4, nz - 1});
  HaloView<Real***, default_memory_pool> const px(
    Kokkos::view_alloc("dpx", Kokkos::WithoutInitializing),
    {0, nx - 1},
    {0, ny - 1},
    {nz - 3, nz - 1});
  HaloView<Real***, default_memory_pool> const py(
    Kokkos::view_alloc("dpy", Kokkos::WithoutInitializing),
    {0, nx - 1},
    {0, ny - 1},
    {nz - 3, nz - 1});
  auto const& ue = this->ue_;
  auto const& ve = this->ve_;
  auto const& we = this->we_;

  auto event   = get_device_event();
  auto stream1 = get_next_stream();
  auto stream2 = get_next_stream();

  // Calculate Φ=J^{-1}p
  const auto& invJ = mesh.invJ;
  parallel_for(
    "scale p_top",
    LoopPolicy<3>(stream1, {0, 0, nz - 4}, {nx, ny, nz}),
    KOKKOS_LAMBDA(int i, int j, int k) {
      invJp(i, j, k) = pp(i, j, k) * invJ(i, j);
    });
  enqueue(event, stream1);

  spectral::ddx(px.view(),
                subview(invJp, ALL, ALL, index_range(nz - 3, nz)).view(),
                mesh.grid,
                stream1);
  wait_for(event, stream2);
  spectral::ddy(py.view(),
                subview(invJp, ALL, ALL, index_range(nz - 3, nz)).view(),
                mesh.grid,
                stream2);
  stream1.fence();
  stream2.fence();

  // estimate J^{-1}u = J^{-1}\hat{u} - dt * grad(J^{-1}p)
  auto const& u = vec_invJu.x;
  auto const& v = vec_invJu.y;
  auto const& w = vec_invJu.z;
  auto const& J = mesh.J;
  parallel_for(
    "estimated u-dt*grad(p)",
    LoopPolicy<2>(stream1, {0, 0}, {nx, ny}),
    KOKKOS_LAMBDA(int i, int j) {
      auto pw_0 = // zw(nz - 4)
        itp2node(invJp(i, j, nz - 3),
                 invJp(i, j, nz - 4),
                 dzw(nz - 4) / 2 / dz(nz - 4));
      auto pw_1 = // zw(nz - 3)
        itp2node(invJp(i, j, nz - 3),
                 invJp(i, j, nz - 2),
                 dzw(nz - 4) / 2 / dz(nz - 3));
      auto pw_2 = invJp(i, j, nz - 2); // zw(nz - 2)

      auto exrij = exr(i, j);
      auto eyrij = eyr(i, j);
      auto dp =
        px(i, j, nz - 3)
        + (zetax(zw(nz - 3), exrij) * pw_1 - zetax(zw(nz - 4), exrij) * pw_0)
            / dzw(nz - 4);
      ue(i, j, 0) = u(i, j, nz - 3) - dt * dp; // J^{-1}u at zz(nz - 3)

      dp = py(i, j, nz - 3)
         + (zetay(zw(nz - 3), eyrij) * pw_1 - zetay(zw(nz - 4), eyrij) * pw_0)
             / dzw(nz - 4);
      ve(i, j, 0) = v(i, j, nz - 3) - dt * dp; // J^{-1}v at zz(nz - 3)

      dp = px(i, j, nz - 2)
         + (zetax(zw(nz - 2), exrij) * pw_2 - zetax(zw(nz - 3), exrij) * pw_1)
             / dzw(nz - 3);
      ue(i, j, 1) = u(i, j, nz - 2) - dt * dp; // J^{-1}u at zz(nz - 2)

      dp = py(i, j, nz - 2)
         + (zetay(zw(nz - 2), eyrij) * pw_2 - zetay(zw(nz - 3), eyrij) * pw_1)
             / dzw(nz - 3);
      ve(i, j, 1) = v(i, j, nz - 2) - dt * dp; // J^{-1}v at zz(nz - 2)

      auto dz0    = dz(nz - 2);
      auto dz1    = dz(nz - 3);
      auto alpha  = dz0 / (dz0 + dz1);
      auto coeff2 = (1 + alpha) / dz0;
      auto coeff1 = -1 / alpha / dz1;
      auto coeff0 = alpha / dz1;
      px(i, j, nz - 1) +=
        zetax(zz(nz - 3), exrij) * invJp(i, j, nz - 3) * coeff0
        + zetax(zz(nz - 2), exrij) * invJp(i, j, nz - 2) * coeff1
        + zetax(zz(nz - 1), exrij) * invJp(i, j, nz - 1) * coeff2;
      py(i, j, nz - 1) +=
        zetay(zz(nz - 3), eyrij) * invJp(i, j, nz - 3) * coeff0
        + zetay(zz(nz - 2), eyrij) * invJp(i, j, nz - 2) * coeff1
        + zetay(zz(nz - 1), eyrij) * invJp(i, j, nz - 1) * coeff2;

      auto pz     = (invJp(i, j, nz - 2) - invJp(i, j, nz - 3)) / dz(nz - 3);
      we(i, j, 0) = w(i, j, nz - 3) - dt * pz * J(i, j); // J^{-1}w at zw(nz
                                                         // - 3)
    });

  spectral::dealias(ue, mesh.grid, stream1);
  spectral::dealias(ve, mesh.grid, stream1);
  spectral::dealias(
    subview(we, ALL, ALL, index_range(0, 1)), mesh.grid, stream1);

  // use d(J^{-1}u)/dzeta and d(J^{-1}v)/dzeta to extrapolate boundary velocity
  parallel_for(
    "compute estimated u and v",
    LoopPolicy<2>(stream1, {0, 0}, {nx, ny}),
    KOKKOS_LAMBDA(int i, int j) {
      auto alpha  = dz(nz - 2) + dz(nz - 3);
      auto beta   = alpha + dz(nz - 2);
      auto coeff1 = alpha * alpha / dz(nz - 3) / beta;
      auto coeff2 = -dz(nz - 2) * dz(nz - 2) / dz(nz - 3) / beta;
      auto coeff0 = dz(nz - 2) * alpha / beta;
      ue(i, j, 2) = coeff1 * ue(i, j, 1) + coeff2 * ue(i, j, 0)
                  + coeff0 * vec_invJuz(i, j, 0);
      ve(i, j, 2) = coeff1 * ve(i, j, 1) + coeff2 * ve(i, j, 0)
                  + coeff0 * vec_invJuz(i, j, 1);
      u(i, j, nz - 1) = ue(i, j, 2) + px(i, j, nz - 1) * dt; // u_hat at top
      v(i, j, nz - 1) = ve(i, j, 2) + py(i, j, nz - 1) * dt; // v_hat at top
    });

  // // ================================================================
  // // An alternative BC implementation, same as the original FST code
  // get_surface_invJw(subview(w, ALL, ALL, index_range(nz - 3, nz -
  // 1)).view(),
  //                   subview(u, ALL, ALL, index_range(nz - 3, nz)).view(),
  //                   subview(v, ALL, ALL, index_range(nz - 3, nz)).view(),
  //                   mesh,
  //                   stream1);
  // parallel_for(
  //   "compute estimated w",
  //   LoopPolicy<2>(stream1, {0, 0}, {nx, ny}),
  //   KOKKOS_LAMBDA(int i, int j) {
  //     auto dz0    = dz(nz - 2);
  //     auto dz1    = dz(nz - 3);
  //     auto alpha  = dz0 / (dz0 + dz1);
  //     auto coeff2 = (1 + alpha) / dz0;
  //     auto coeff1 = -1 / alpha / dz1;
  //     auto coeff0 = alpha / dz1;
  //
  //     auto pz = invJp(i, j, nz - 3) * coeff0 + invJp(i, j, nz - 2) *
  //     coeff1
  //             + invJp(i, j, nz - 1) * coeff2;
  //     we(i, j, 1) = w(i, j, nz - 2) - dt * pz * J(i, j);
  //   });
  // // ================================================================

  get_surface_invJw(we, ue, ve, mesh, stream1);
  parallel_for(
    "compute estimated w",
    LoopPolicy<2>(stream1, {0, 0}, {nx, ny}),
    KOKKOS_LAMBDA(int i, int j) {
      // w = \hat{w} seems to be necessary for stability
      w(i, j, nz - 2) = we(i, j, 1); // + dt * pz * J(i, j);
    });

  stream1.fence();
}

void FreeSurfaceBC::get_surface_invJw(
  MDView<Real** [2]> const&            invJw,
  MDView<Real const** [3]> const&      invJu,
  MDView<Real const** [3]> const&      invJv,
  TopWaveMesh const&                   mesh,
  Kokkos::DefaultExecutionSpace const& stream)
{
  using Kokkos::parallel_for;

  MDView<Real** [1]> const ux_vy("ux+vy", invJu.extent(0), invJu.extent(1));
  spectral::ddx(
    ux_vy, subview(invJu, ALL, ALL, index_range(1, 2)), mesh.grid, stream);
  spectral::ddy_and_add(
    ux_vy, subview(invJv, ALL, ALL, index_range(1, 2)), mesh.grid, stream);

  auto const& dz  = mesh.dz;
  auto const& dzw = mesh.dzw;
  auto const& zw  = mesh.zw;
  auto const& J   = mesh.J;
  auto const& exr = mesh.exr;
  auto const& eyr = mesh.eyr;
  auto const  nz  = mesh.extent(2);
  parallel_for(
    "compute estimated w",
    LoopPolicy<2>(stream, {0, 0}, {invJu.extent(0), invJu.extent(1)}),
    KOKKOS_LAMBDA(int i, int j) {
      auto invJu1 = // zw(nz - 3)
        itp2node(invJu(i, j, 1), invJu(i, j, 0), dzw(nz - 3) / 2 / dz(nz - 3));
      auto invJv1 = // zw(nz - 3)
        itp2node(invJv(i, j, 1), invJv(i, j, 0), dzw(nz - 3) / 2 / dz(nz - 3));
      auto invJW1 = zetax(zw(nz - 3), exr(i, j)) * invJu1
                  + zetay(zw(nz - 3), eyr(i, j)) * invJv1
                  + J(i, j) * invJw(i, j, 0);
      auto invJW0    = (0 - ux_vy(i, j, 0)) * dzw(nz - 3) + invJW1;
      auto invJu0    = zetax(zw(nz - 2), exr(i, j)) * invJu(i, j, 2);
      auto invJv0    = zetay(zw(nz - 2), eyr(i, j)) * invJv(i, j, 2);
      invJw(i, j, 1) = (invJW0 - invJu0 - invJv0) / J(i, j);
    });
  spectral::dealias(
    subview(invJw, ALL, ALL, index_range(1, 2)), mesh.grid, stream);
  stream.fence();
}

void FreeSurfaceBC::get_surface_pressure(MDView<Real**> const&       p_top,
                                         MDView<Real const**> const& eta,
                                         Real                        Fr2,
                                         Real                        RWe,
                                         Real                        nu,
                                         TopWaveMesh const&          mesh,
                                         bool update_stored_p_top) const
{
  using Kokkos::parallel_for;

  auto const is_top = mesh.comm().is_last(2);
  if (!is_top) return;

  auto const  nx   = mesh.extent(0);
  auto const  ny   = mesh.extent(1);
  auto const& grid = mesh.grid;

  auto event   = get_device_event();
  auto stream1 = get_next_stream();
  auto stream2 = get_next_stream();

  // views for velocity at the surface
  MDView<Real** [3]> const u_xi("u_xi_surface", nx, ny);
  MDView<Real** [3]> const u_psi("u_psi_surface", nx, ny);
  // copy top velocity for computing derivatives
  // the velocity is scaled by J because the stored velocity is J^{-1}u
  auto const  policy = LoopPolicy<2>(stream1, {0, 0}, {nx, ny});
  auto const& J      = mesh.J;
  auto const& ue     = subview(this->ue_, ALL, ALL, ue_.extent(2) - 1);
  auto const& ve     = subview(this->ve_, ALL, ALL, ve_.extent(2) - 1);
  auto const& we     = subview(this->we_, ALL, ALL, we_.extent(2) - 1);
  parallel_for(
    "[pbc] copy surface velocity", policy, KOKKOS_LAMBDA(int i, int j) {
      u_xi(i, j, 0) = ue(i, j) * J(i, j);
      u_xi(i, j, 1) = ve(i, j) * J(i, j);
      u_xi(i, j, 2) = we(i, j) * J(i, j);
    });
  spectral::dealias(u_xi, grid, stream1);
  enqueue(event, stream1);
  wait_for(event, stream2); // wait for dealias
  Kokkos::deep_copy(stream2, u_psi, u_xi);

  spectral::ddx(u_xi, u_xi, grid, stream1);
  spectral::ddy(u_psi, u_psi, grid, stream2);
  enqueue(event, stream2);
  wait_for(event, stream1); // wait for ddy

  // p_top with viscous stress
  auto const& ex = mesh.ex;
  auto const& ey = mesh.ey;
  parallel_for(
    policy, KOKKOS_LAMBDA(int i, int j) {
      // eta_x*eta_y*(u_psi+v_xi)
      auto p1 = ex(i, j) * ey(i, j) * (u_psi(i, j, 0) + u_xi(i, j, 1));
      // -eta_x*w_xi - eta_y*w_psi
      auto p2 = -ex(i, j) * u_xi(i, j, 2) - ey(i, j) * u_psi(i, j, 2);
      // -(1+eta_x^2)*v_psi - (1+eta_y^2)*u_xi
      auto p3 = -(1 + ex(i, j) * ex(i, j)) * u_psi(i, j, 1)
              - (1 + ey(i, j) * ey(i, j)) * u_xi(i, j, 0);
      auto reh2   = 1 + (ex(i, j) * ex(i, j) + ey(i, j) * ey(i, j));
      p_top(i, j) = (p1 + p2 + p3) / reh2 * 2 * nu + eta(i, j) / Fr2;
    });
  stream1.fence();

  // surface tension
  if (RWe > 0) {
    MDView<Real** [3], default_memory_pool> const e2(
      Kokkos::view_alloc("eta 2nd derivatives", Kokkos::WithoutInitializing),
      nx,
      ny);
    MDView<Real** [1]> const ex1(mesh.ex.data(), nx, ny);
    MDView<Real** [1]> const ey1(mesh.ey.data(), nx, ny);
    // eta_xx
    spectral::ddx(subview(e2, ALL, ALL, index_range(0, 1)), ex1, grid, stream1);
    // eta_xy
    spectral::ddx(subview(e2, ALL, ALL, index_range(1, 2)), ey1, grid, stream1);
    // eta_yy
    spectral::ddy(subview(e2, ALL, ALL, index_range(2, 3)), ey1, grid, stream1);

    parallel_for(
      "add surface tension", policy, KOKKOS_LAMBDA(int i, int j) {
        auto kappa = (1 + ex(i, j) * ex(i, j)) * e2(i, j, 2)
                   + (1 + ey(i, j) * ey(i, j)) * e2(i, j, 0)
                   - 2 * ex(i, j) * ey(i, j) * e2(i, j, 1);
        auto reh  = alps::rnorm3d(1, ex(i, j), ey(i, j));
        auto reh2 = 1 + (ex(i, j) * ex(i, j) + ey(i, j) * ey(i, j));
        p_top(i, j) -= RWe * kappa / (reh2 * reh);
      });

    stream1.fence();
  }

  spectral::dealias(MDView<Real** [1]>(Kokkos::view_wrap(p_top.data()), nx, ny),
                    mesh.grid,
                    stream1);
  auto const& p_top_0_ = this->p_top_0;
  if (update_stored_p_top) {
    parallel_for(
      "apply p_top", policy, KOKKOS_LAMBDA(int i, int j) {
        auto p_n       = p_top(i, j);
        p_top(i, j)    = (p_n + p_top_0_(i, j)) / 2;
        p_top_0_(i, j) = p_n;
      });
  } else {
    parallel_for(
      "apply p_top", policy, KOKKOS_LAMBDA(int i, int j) {
        p_top(i, j) = (p_top(i, j) + p_top_0_(i, j)) / 2;
      });
  }

  stream1.fence();
}

FreeSurfaceBCOptions FreeSurfaceBCOptions::parse_from(ConfigTable const& config)
{
  FreeSurfaceBCOptions opts;
  opts.hyper_viscosity_n = config.get_value_or("hyperviscosity_n", 2);
  if (opts.hyper_viscosity_n < 1 || opts.hyper_viscosity_n > 5) {
    throw std::invalid_argument("Invalid hyper viscosity power");
  }

  opts.hyper_viscosity_c = config.get_value_or("hyperviscosity_c", 1.0);

  return opts;
}

void FreeSurfaceSolver::apply_bottom_bc(
  const Kokkos::DefaultExecutionSpace& space) const
{
  bool        bc_set = false;
  auto const& flow   = this->flow_field;
  if (auto const* bc = dynamic_cast<NoSlipWall const*>(flow.bottom_bc.get());
      bc != nullptr) {
    solver::apply_bottom_bc(flow, *bc, space);
    bc_set = true;
  }
  if (auto const* bc = dynamic_cast<GradientWall const*>(flow.bottom_bc.get());
      bc != nullptr) {
    solver::apply_bottom_bc(flow, *bc, space);
    bc_set = true;
  }

  if (!bc_set) {
    logger->warn("No bottom boundary condition set.");
  }

  space.fence();
}

} // namespace alps::solver
