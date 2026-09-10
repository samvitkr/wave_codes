#include <common/container/matrix_field.h>
#include <common/container/view_utils.h>
#include <common/kokkos_abstraction/exec_policy.h>
#include <common/kokkos_abstraction/pool_space.h>
#include <common/runtime/async_utils.h>
#include <solvers/mesh/mesh.h>
#include <solvers/operators/grad.h>
#include <spectral/spectral.h>

namespace alps::solver {

void calculate_strain_rate(const SymmTensor33Field<Real***>& Sij,
                           const Vector3Field<Real***>&      u_vec,
                           const Mesh&                       mesh)
{
  using Kokkos::parallel_for;

  const auto& grid = mesh.grid;
  const auto  hbar = mesh.hbar;
  const auto& dz   = mesh.dz;
  const auto& dzw  = mesh.dzw;
  //  const auto  is_top    = grid.comm().is_last(2);
  //  const auto  is_bottom = grid.comm().is_first(2);

  const auto u      = u_vec.x;
  const auto v      = u_vec.y;
  const auto w      = u_vec.z;
  const auto begins = local_begins(u);
  const auto ends   = local_ends(u);
  const auto nz     = local_extent(u, 2);

  MDView<Real***, default_memory_pool> tmp_x("Sij_tmpx",
                                             create_local_layout(grid.pencil));

  auto stream1        = get_next_stream();
  auto stream2        = get_next_stream();
  auto constexpr tile = []() -> Kokkos::Array<std::int64_t, 3> {
    using Device = Kokkos::DefaultExecutionSpace;
    if constexpr (is_cuda_execution_space_v<Device>) return {16, 1, 8};
    if constexpr (is_hip_execution_space_v<Device>) return {64, 1, 4};
    return {0, 0, 0};
  }();

  {
    auto u_inner = create_inner_view(u).view();
    spectral::ddx(create_inner_view(Sij.xx).view(), u_inner, grid, stream1);
    spectral::ddy(create_inner_view(Sij.xy).view(), u_inner, grid, stream2);

    const auto&   Sxz = Sij.xz;
    LoopPolicy<3> policy(stream1, begins, ends, tile);
    parallel_for(
      "Sij add uz", policy, KOKKOS_LAMBDA(int i, int j, int k) {
        auto alpha   = 1 / (dz(k) * hbar);
        Sxz(i, j, k) = (u(i, j, k + 1) - u(i, j, k)) * alpha;
      });
    stream1.fence();
    stream2.fence();
  }

  {
    auto u_inner = create_inner_view(v).view();
    spectral::ddx(tmp_x, u_inner, grid, stream1);
    spectral::ddy(create_inner_view(Sij.yy).view(), u_inner, grid, stream2);

    const auto&   Syx = Sij.xy;
    const auto&   Syz = Sij.yz;
    LoopPolicy<3> policy(stream1, begins, ends, tile);
    parallel_for(
      "Sij add gradv", policy, KOKKOS_LAMBDA(int i, int j, int k) {
        auto alpha   = 1 / (dz(k) * hbar);
        Syx(i, j, k) = (Syx(i, j, k) + tmp_x(i, j, k)) / 2;
        Syz(i, j, k) = (v(i, j, k + 1) - v(i, j, k)) * alpha;
      });
    stream1.fence();
    stream2.fence();
  }

  // Before finishing S13 and S23, first re-calculate the boundary values of
  // du/dz and dv/dz using one-side stencils
  if (grid.comm().is_first(2)) {
    for (const auto& [f, flux, dd] :
         {std::tie(u, Sij.xz, mesh.dz_h), std::tie(v, Sij.yz, mesh.dz_h)}) {
      using functor_t = FDBoundaryFunctor<Real, Kokkos::DefaultExecutionSpace>;
      auto policy =
        LoopPolicy<2, functor_t::Add>(stream1, {0, 0}, {ends[0], ends[1]});
      functor_t functor(subview(flux, ALL, ALL, 0).view(),
                        subview(f, ALL, ALL, index_range(0, 3)).view(),
                        subview(dd, index_range(0, 2)).view(),
                        mesh.hbar,
                        LeftBoundary());
      parallel_for(policy, functor);
    }
  }
  if (grid.comm().is_last(2)) {
    for (const auto& [f, flux, dd, flux_idx, f_idx] :
         {std::make_tuple(u, Sij.xz, mesh.dz_h, nz - 2, nz - 3),
          std::make_tuple(v, Sij.yz, mesh.dz_h, nz - 2, nz - 3)}) {
      using functor_t = FDBoundaryFunctor<Real, Kokkos::DefaultExecutionSpace>;
      auto policy =
        LoopPolicy<2, functor_t::Add>(stream2, {0, 0}, {ends[0], ends[1]});
      functor_t functor(
        subview(flux, ALL, ALL, flux_idx).view(),
        subview(f, ALL, ALL, index_range(f_idx, f_idx + 3)).view(),
        subview(dd, index_range(f_idx, f_idx + 2)).view(),
        mesh.hbar,
        RightBoundary());
      parallel_for(policy, functor);
    }
  }
  stream1.fence();
  stream2.fence();

  {
    auto u_inner = create_inner_view(w).view();
    spectral::ddx(tmp_x, u_inner, grid, stream1);
    spectral::ddy_and_add(
      create_inner_view(Sij.zy).view(), u_inner, grid, stream2);
    auto event = get_device_event();
    enqueue(event, stream2);
    wait_for(
      event,
      stream1); // stream1 must wait for the computation of dw/dy to finish

    const auto&   Szx = Sij.xz;
    const auto&   Szy = Sij.yz;
    const auto&   Szz = Sij.zz;
    LoopPolicy<3> policy(stream1, begins, ends, tile);
    parallel_for(
      "Sij add gradw", policy, KOKKOS_LAMBDA(int i, int j, int k) {
        Szx(i, j, k) = (Szx(i, j, k) + tmp_x(i, j, k)) / 2;
        Szy(i, j, k) /= 2;
        auto beta    = 1 / (dzw(k - 1) * hbar);
        Szz(i, j, k) = (w(i, j, k) - w(i, j, k - 1)) * beta;
      });
    stream1.fence();
  }
}
} // namespace alps::solver
