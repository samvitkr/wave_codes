#include "scalar_advection_diffusion.h"

#include <common/container/vector_field.h>
#include <common/container/view_utils.h>
#include <common/kokkos_abstraction/exec_policy.h>
#include <common/kokkos_abstraction/pool_space.h>
#include <common/runtime/async_utils.h>
#include <linear_algebra/axpy.h>
#include <solvers/field/flow_field.h>
#include <solvers/mesh/mesh.h>
#include <solvers/ns/advection_invoke.h>
#include <solvers/operators/grad.h>
#include <spectral/spectral.h>

#include <Kokkos_Core.hpp>
#include <Kokkos_Profiling_ScopedRegion.hpp>

namespace alps::solver {

struct ChannelAdvectionFunctor
{
  struct Bottom_0
  {}; // k == 0
  struct Top_NzMinus2
  {}; // k == nz - 2
  struct Top_NzMinus1
  {}; // k == nz - 1
  struct General
  {};

  ChannelAdvectionFunctor(const Vector3Field<Real***>&   fluxes,
                          const Mesh&                    mesh,
                          const Vector3Field<Real***>&   vec_u,
                          const HaloView<Real const***>& scalar)
    : fx{fluxes.x}
    , fy{fluxes.y}
    , fz{fluxes.z}
    , u{vec_u.x}
    , v{vec_u.y}
    , w{vec_u.z}
    , f{scalar}
    , dzw{mesh.dzw}
    , nz{mesh.extent(2)}
  {
    if (!check_same_layout_and_offset(fx, fy, fz, u, v, w, f)) {
      throw std::runtime_error(
        "Mismatched extents in vec_u, scalar and fluxes");
    }
  }

  using view_t       = HaloView<Real***>;
  using const_view_t = view_t::const_type;
  view_t       fx, fy, fz;
  const_view_t u, v, w;
  const_view_t f;

  HaloView<Real const*> dzw;
  int                   nz;

  KOKKOS_FUNCTION void operator()(const Bottom_0& /*tag*/, int i, int j) const
  {
    fx(i, j, 0) = -u(i, j, 0) * f(i, j, 0);
    fy(i, j, 0) = -v(i, j, 0) * f(i, j, 0);
    fz(i, j, 0) = -w(i, j, 0) * f(i, j, 0);
  }

  KOKKOS_FUNCTION void
  operator()(const Top_NzMinus2& /*tag*/, int i, int j) const
  {
    fx(i, j, nz - 2) = -u(i, j, nz - 2) * f(i, j, nz - 2);
    fy(i, j, nz - 2) = -v(i, j, nz - 2) * f(i, j, nz - 2);
    fz(i, j, nz - 2) = -w(i, j, nz - 2) * f(i, j, nz - 2 + 1);
  }

  KOKKOS_FUNCTION void
  operator()(const Top_NzMinus1& /*tag*/, int i, int j) const
  {
    fx(i, j, nz - 1) = -u(i, j, nz - 1) * f(i, j, nz - 1);
    fy(i, j, nz - 1) = -v(i, j, nz - 1) * f(i, j, nz - 1);
  }

  KOKKOS_FUNCTION void
  operator()(const General& /*tag*/, int i, int j, int k) const
  {
    auto offset       = &f(i, j, k) - f.data();
    fx.data()[offset] = -u.data()[offset] * f.data()[offset];
    fy.data()[offset] = -v.data()[offset] * f.data()[offset];
    auto beta         = dzw(k - 1) / (dzw(k - 1) + dzw(k));
    fz.data()[offset] =
      -w.data()[offset] * itp2node(f.data()[offset], f(i, j, k + 1), beta);
  }
};

} // namespace alps::solver

namespace alps::solver {
void calculate_advection_fluxes(const Vector3Field<Real***>&   fluxes,
                                const HaloView<Real const***>& f,
                                const FlowField&               flow)
{
  invoke_advection_fluxes_functor<ChannelAdvectionFunctor>(fluxes, f, flow);
}

void add_diffusion_fluxes(Vector3Field<Real***> const&   fluxes,
                          const HaloView<Real const***>& f,
                          const FlowField&               flow,
                          Real                           D,
                          Real                           gamma)
{
  using Kokkos::parallel_for;
  using std::make_tuple;
  using std::tie;

  auto const region = Kokkos::Profiling::ScopedRegion("diffusion " + f.label());

  const Mesh& mesh        = flow.mesh;
  const auto& grid        = mesh.grid;
  const auto& [u, v, w]   = tie(flow.u.x, flow.u.y, flow.u.z);
  const auto [nx, ny, nz] = local_extents(u);
  const auto coeff        = D * gamma;
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

  {
    LoopPolicy<1> axpy_policy_x(stream1, 0, tmp3x.span());
    LoopPolicy<1> axpy_policy_y(stream2, 0, tmp3y.span());
    auto          inner        = create_inner_view(f).view();
    auto          flux_x_inner = create_inner_view(fluxes.x).view();
    auto          flux_y_inner = create_inner_view(fluxes.y).view();
    spectral::ddx(tmp3x, inner, grid, stream1);
    axpy_with_policy(coeff, tmp3x, flux_x_inner, axpy_policy_x);
    spectral::ddy(tmp3y, inner, grid, stream2);
    axpy_with_policy(coeff, tmp3y, flux_y_inner, axpy_policy_y);
  }

  if (is_bottom) {
    using functor_t = FDBoundaryFunctor<Real, Kokkos::DefaultExecutionSpace>;
    auto policy =
      LoopPolicy<2, functor_t::Add>(get_next_stream(), {0, 0}, {nx, ny});
    functor_t functor(subview(fluxes.z, ALL, ALL, 0).view(),
                      subview(f, ALL, ALL, index_range(0, 3)).view(),
                      subview(mesh.dz_h, index_range(0, 2)).view(),
                      mesh.hbar / coeff,
                      LeftBoundary());
    parallel_for(policy, functor);
    streams.push_back(policy.space());
  }
  if (is_top) {
    using functor_t = FDBoundaryFunctor<Real, Kokkos::DefaultExecutionSpace>;
    auto policy =
      LoopPolicy<2, functor_t::Add>(get_next_stream(), {0, 0}, {nx, ny});
    functor_t functor(
      subview(fluxes.z, ALL, ALL, nz - 2).view(),
      subview(f, ALL, ALL, index_range(nz - 3, nz)).view(),
      subview(mesh.dz_h, index_range(nz - 3, nz - 3 + 2)).view(),
      mesh.hbar / coeff,
      RightBoundary());
    parallel_for(policy, functor);
    streams.push_back(policy.space());
  }

  constexpr auto tile_size = []() -> Kokkos::Array<std::int64_t, 3> {
    using Device = Kokkos::DefaultExecutionSpace;
    if constexpr (is_cuda_execution_space_v<Device>) return {16, 1, 8};
    if constexpr (is_hip_execution_space_v<Device>) return {64, 1, 4};
    return {0, 0, 0};
  }();

  const int   z_begin = is_bottom ? 1 : 0;
  const int   z_end   = is_top ? nz - 2 : nz;
  const auto& flux    = fluxes.z;
  const auto& dd      = mesh.dz;
  auto        policy  = LoopPolicy<3>(
    get_next_stream(), {0, 0, z_begin}, {nx, ny, z_end}, tile_size);
  const auto hbar = mesh.hbar;
  parallel_for(
    "diffusion " + f.label() + "_z",
    policy,
    KOKKOS_LAMBDA(int i, int j, int k) {
      auto alpha = coeff / (dd(k) * hbar);
      flux(i, j, k) += alpha * (f(i, j, k + 1) - f(i, j, k));
    });
  streams.push_back(policy.space());

  fence(std::move(streams));
}

} // namespace alps::solver
