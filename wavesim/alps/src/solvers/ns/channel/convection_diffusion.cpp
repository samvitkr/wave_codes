#include "convection_diffusion.h"

#include <common/container/matrix_field.h>
#include <common/container/view_utils.h>
#include <common/device/device_traits.h>
#include <common/kokkos_abstraction/exec_policy.h>
#include <common/kokkos_abstraction/pool_space.h>
#include <common/runtime/async_utils.h>
#include <linear_algebra/axpy.h>
#include <solvers/field/flow_field.h>
#include <solvers/mesh/mesh.h>
#include <solvers/ns/convection_invoke.h>
#include <solvers/operators/grad.h>
#include <spectral/spectral.h>

#include <Kokkos_Core.hpp>
#include <Kokkos_Profiling_ScopedRegion.hpp>

namespace alps::solver {

struct ChannelConvectionFunctor
{
  struct Bottom_0
  {}; // k == 0
  struct Top_NzMinus2
  {}; // k == nz - 2
  struct Top_NzMinus1
  {}; // k == nz - 1
  struct General
  {};

  ChannelConvectionFunctor(const Tensor33Field<Real***>& fluxes,
                           const Mesh&                   mesh,
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
    , dzw{mesh.dzw}
    , nz{mesh.extent(2)}
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

  HaloView<Real const*> dzw;
  int                   nz;

  KOKKOS_FUNCTION void operator()(const Bottom_0& /*tag*/, int i, int j) const
  {
    fxx(i, j, 0) = -u(i, j, 0) * u(i, j, 0);
    auto r       = -u(i, j, 0) * v(i, j, 0);
    fxy(i, j, 0) = r;
    fyx(i, j, 0) = r;
    fyy(i, j, 0) = -v(i, j, 0) * v(i, j, 0);
    fxz(i, j, 0) = -w(i, j, 0) * u(i, j, 0);
    fzx(i, j, 0) = fxz(i, j, 0);
    fyz(i, j, 0) = -w(i, j, 0) * v(i, j, 0);
    fzy(i, j, 0) = fyz(i, j, 0);
    fzz(i, j, 0) = -w(i, j, 0) * w(i, j, 0);
  }

  KOKKOS_FUNCTION void
  operator()(const Top_NzMinus2& /*tag*/, int i, int j) const
  {
    fxx(i, j, nz - 2) = -u(i, j, nz - 2) * u(i, j, nz - 2);
    auto r            = -u(i, j, nz - 2) * v(i, j, nz - 2);
    fxy(i, j, nz - 2) = r;
    fyx(i, j, nz - 2) = r;
    fyy(i, j, nz - 2) = -v(i, j, nz - 2) * v(i, j, nz - 2);
    r                 = -w(i, j, nz - 2);
    fxz(i, j, nz - 2) = r * u(i, j, nz - 2 + 1);
    fzx(i, j, nz - 2) = fxz(i, j, nz - 2);
    fyz(i, j, nz - 2) = r * v(i, j, nz - 2 + 1);
    fzy(i, j, nz - 2) = fyz(i, j, nz - 2);
    r                 = (w(i, j, nz - 2) + w(i, j, nz - 2 - 1)) / 2;
    fzz(i, j, nz - 2) = -r * r;
  }

  KOKKOS_FUNCTION void
  operator()(const Top_NzMinus1& /*tag*/, int i, int j) const
  {
    fxx(i, j, nz - 1) = -u(i, j, nz - 1) * u(i, j, nz - 1);
    auto r            = -u(i, j, nz - 1) * v(i, j, nz - 1);
    fxy(i, j, nz - 1) = r;
    fyx(i, j, nz - 1) = r;
    fyy(i, j, nz - 1) = -v(i, j, nz - 1) * v(i, j, nz - 1);
    fzz(i, j, nz - 1) = -w(i, j, nz - 1 - 1) * w(i, j, nz - 1 - 1);
  }

  KOKKOS_FUNCTION void
  operator()(const General& /*tag*/, int i, int j, int k) const
  {
    auto offset = &u(i, j, k) - u.data(); // (i, j, k)
    auto uijk   = u.data()[offset];
    auto vijk   = v.data()[offset];
    auto wijk   = w.data()[offset];

    fxx.data()[offset] = -uijk * uijk;
    auto r             = -uijk * vijk;
    fxy.data()[offset] = r;
    fyx.data()[offset] = r;
    fyy.data()[offset] = -vijk * vijk;
    r                  = wijk + w.data()[offset - w.stride(2)];
    fzz.data()[offset] = -r * r / 4;
    auto beta          = dzw(k - 1) / (dzw(k - 1) + dzw(k));
    fxz.data()[offset] =
      -wijk
      * Kokkos::fma(beta,
                    u.data()[offset + u.stride(2)], // (i, j, k + 1)
                    Kokkos::fma(-beta, uijk, uijk));
    fzx.data()[offset] = fxz.data()[offset];
    fyz.data()[offset] =
      -wijk
      * Kokkos::fma(beta,
                    v.data()[offset + v.stride(2)], // (i, j, k + 1)
                    Kokkos::fma(-beta, vijk, vijk));
    fzy.data()[offset] = fyz.data()[offset];
  }
};

} // namespace alps::solver

namespace alps::solver {
void calculate_convection_fluxes(const Tensor33Field<Real***>& fluxes,
                                 const FlowField&              flow)
{
  invoke_convection_fluxes_functor<ChannelConvectionFunctor>(fluxes, flow);
}

void add_viscous_fluxes_from_u(Tensor33Field<Real***> const& fluxes,
                               const FlowField&              flow,
                               Real                          nu,
                               Real                          gamma)
{
  using Kokkos::parallel_for;
  using std::make_tuple;
  using std::tie;

  auto const region = Kokkos::Profiling::ScopedRegion("viscous");

  const Mesh& mesh        = flow.mesh;
  const auto& grid        = mesh.grid;
  const auto& [u, v, w]   = tie(flow.u.x, flow.u.y, flow.u.z);
  const auto [nx, ny, nz] = local_extents(u);
  const auto coeff        = nu * gamma;
  const auto is_top       = grid.comm().is_last(2);
  const auto is_bottom    = grid.comm().is_first(2);

  auto        stream1 = get_next_stream();
  auto        stream2 = get_next_stream();
  std::vector streams{stream1, stream2};

  MDView<Real***, default_memory_pool> tmp3x(
    Kokkos::view_alloc("tmp_x", Kokkos::WithoutInitializing),
    create_local_layout(grid.pencil));
  MDView<Real***, default_memory_pool> tmp3y(
    Kokkos::view_alloc("tmp_y", Kokkos::WithoutInitializing),
    create_local_layout(grid.pencil));

  LoopPolicy<1> axpy_policy_x(stream1, 0, tmp3x.span());
  LoopPolicy<1> axpy_policy_y(stream2, 0, tmp3y.span());

  for (const auto& [f, flux_x, flux_y] : {tie(u, fluxes.xx, fluxes.xy),
                                          tie(v, fluxes.yx, fluxes.yy),
                                          tie(w, fluxes.zx, fluxes.zy)}) {
    auto inner        = create_inner_view(f).view();
    auto flux_x_inner = create_inner_view(flux_x).view();
    auto flux_y_inner = create_inner_view(flux_y).view();
    spectral::ddx(tmp3x, inner, grid, stream1);
    axpy_with_policy(coeff, tmp3x, flux_x_inner, axpy_policy_x);
    spectral::ddy(tmp3y, inner, grid, stream2);
    axpy_with_policy(coeff, tmp3y, flux_y_inner, axpy_policy_y);
  }

  if (is_bottom) {
    for (const auto& [f, flux, dd] : {tie(u, fluxes.xz, mesh.dz_h),
                                      tie(v, fluxes.yz, mesh.dz_h),
                                      tie(w, fluxes.zz, mesh.dzw_h)}) {
      using functor_t = FDBoundaryFunctor<Real, Kokkos::DefaultExecutionSpace>;
      auto policy =
        LoopPolicy<2, functor_t::Add>(get_next_stream(), {0, 0}, {nx, ny});
      functor_t functor(subview(flux, ALL, ALL, 0).view(),
                        subview(f, ALL, ALL, index_range(0, 3)).view(),
                        subview(dd, index_range(0, 2)).view(),
                        mesh.hbar / coeff,
                        LeftBoundary());
      parallel_for(policy, functor);
      streams.push_back(policy.space());
    }
  }
  if (is_top) {
    for (const auto& [f, flux, dd, flux_idx, f_idx] :
         {make_tuple(u, fluxes.xz, mesh.dz_h, nz - 2, nz - 3),
          make_tuple(v, fluxes.yz, mesh.dz_h, nz - 2, nz - 3),
          make_tuple(w, fluxes.zz, mesh.dzw_h, nz - 1, nz - 4)}) {
      using functor_t = FDBoundaryFunctor<Real, Kokkos::DefaultExecutionSpace>;
      auto policy =
        LoopPolicy<2, functor_t::Add>(get_next_stream(), {0, 0}, {nx, ny});
      functor_t functor(
        subview(flux, ALL, ALL, flux_idx).view(),
        subview(f, ALL, ALL, index_range(f_idx, f_idx + 3)).view(),
        subview(dd, index_range(f_idx, f_idx + 2)).view(),
        mesh.hbar / coeff,
        RightBoundary());
      parallel_for(policy, functor);
      streams.push_back(policy.space());
    }
  }

  auto constexpr tile_size = []() -> Kokkos::Array<range_index_t, 3> {
    using Device = Kokkos::DefaultExecutionSpace;
    if constexpr (is_cuda_execution_space_v<Device>) return {16, 1, 8};
    if constexpr (is_hip_execution_space_v<Device>) return {64, 1, 4};
    return {0, 0, 0};
  }();

  for (const auto& [name, flux_, f_, dd_, z_end] :
       {make_tuple("viscous_xz", fluxes.xz, u, mesh.dz, is_top ? nz - 2 : nz),
        make_tuple("viscous_yz", fluxes.yz, v, mesh.dz, is_top ? nz - 2 : nz),
        // For the node-location variable w, the index is shifted by 1
        // Also, the top boundary node for fluxes.zz is k = nz - 1
        make_tuple("viscous_zz",
                   fluxes.zz,
                   std::decay_t<decltype(w)>(
                     w.view(), {begin(w, 0), begin(w, 1), begin(w, 2) + 1}),
                   HaloView<Real*>(mesh.dzw.view(), {begin(mesh.dzw, 0) + 1}),
                   is_top ? nz - 1 : nz)}) {
    const int   z_begin = is_bottom ? 1 : 0;
    const auto& flux    = flux_;
    const auto& f       = f_;
    const auto& dd      = dd_;
    auto        policy  = LoopPolicy<3>(
      get_next_stream(), {0, 0, z_begin}, {nx, ny, z_end}, tile_size);
    const auto hbar = mesh.hbar;
    parallel_for(
      name, policy, KOKKOS_LAMBDA(int i, int j, int k) {
        auto alpha = coeff / (dd(k) * hbar);
        flux(i, j, k) += alpha * (f(i, j, k + 1) - f(i, j, k));
      });
    streams.push_back(policy.space());
  }

  fence(std::move(streams));
}
} // namespace alps::solver
