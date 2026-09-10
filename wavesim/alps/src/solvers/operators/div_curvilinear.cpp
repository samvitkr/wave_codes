#include "div_curvilinear.h"

#include <common/base/macros.h>
#include <common/container/view_types.h>
#include <common/container/view_utils.h>
#include <common/kokkos_abstraction/exec_policy.h>
#include <common/kokkos_abstraction/reducer.h>
#include <common/runtime/async_utils.h>
#include <decomp/ghost_cell_exchange.h>
#include <spectral/spectral.h>

#include <Kokkos_Profiling_ScopedRegion.hpp>

using Kokkos::parallel_for;

namespace alps::solver {

namespace detail {
template<typename MT>
[[nodiscard]] MDView<Real***> div_impl(const Vector3Field<Real***>& vec_invJu,
                                       const MT&                    mesh,
                                       CenterPt /*unused*/,
                                       std::string_view label)
{
  const auto& u         = vec_invJu.x;
  const auto& v         = vec_invJu.y;
  const auto& w         = vec_invJu.z;
  const auto  is_top    = mesh.grid.comm().is_last(2);
  const auto  is_bottom = mesh.grid.comm().is_first(2);
  const auto  nx        = local_extent(u, 0);
  const auto  ny        = local_extent(u, 1);
  const auto  nz        = local_extent(u, 2);
  const auto& zw        = mesh.zw;
  const auto& dz        = mesh.dz;
  const auto& dzw       = mesh.dzw;
  const auto& exr       = mesh.exr;
  const auto& eyr       = mesh.eyr;
  const auto& J         = mesh.J;

  MDView<Real***, default_memory_pool> div_u(
    Kokkos::view_alloc(std::string(label), Kokkos::WithoutInitializing),
    create_local_layout(mesh.grid.pencil));

  /* update the upper halos of velocities */
  auto requests = async_update_halo_lower_z(mesh.grid.x_pencil(), w, 1);
  requests.push(async_update_halo_z(mesh.grid.x_pencil(), u, 2));
  requests.push(async_update_halo_z(mesh.grid.x_pencil(), v, 3));

  auto stream = get_next_stream();
  // calculate ∂(J^{-1}u)/∂ξ
  spectral::ddx(div_u, create_inner_view(u).view(), mesh.grid, stream);

  // calculate ∂(J^{-1}v)/∂𝜓
  spectral::ddy_and_add(div_u, create_inner_view(v).view(), mesh.grid, stream);

  requests.waitall(); // wait for the velocity halo update

  // Calculate ∂(J^{-1}W)/∂ζ
  constexpr auto tile = []() -> Kokkos::Array<std::int64_t, 3> {
    using Device = Kokkos::DefaultExecutionSpace;
    if constexpr (is_cuda_execution_space_v<Device>) return {16, 1, 8};
    if constexpr (is_hip_execution_space_v<Device>) return {64, 1, 2};
    return {0, 0, 0};
  }();
  LoopPolicy<3> policy(
    stream, {0, 0, is_bottom ? 1 : 0}, {nx, ny, is_top ? nz - 1 : nz}, tile);
  parallel_for(
    "dW/dzeta", policy, KOKKOS_LAMBDA(int i, int j, int k) {
      auto inv_dzk  = 1 / dz(k);
      auto inv_dzk1 = 1 / dz(k - 1);
      auto beta     = dzw(k - 1) / 2 * inv_dzk;
      auto u_w      = itp2node(u(i, j, k), u(i, j, k + 1), beta);
      auto v_w      = itp2node(v(i, j, k), v(i, j, k + 1), beta);
      auto beta1    = dzw(k - 1) / 2 * inv_dzk1;
      auto u_w1     = itp2node(u(i, j, k), u(i, j, k - 1), beta1);
      auto v_w1     = itp2node(v(i, j, k), v(i, j, k - 1), beta1);
      auto tmp_k    = MT::zeta_x(zw(k), exr(i, j)) * u_w
                 + MT::zeta_y(zw(k), eyr(i, j)) * v_w + J(i, j) * w(i, j, k);
      auto tmp_k1 = MT::zeta_x(zw(k - 1), exr(i, j)) * u_w1
                  + MT::zeta_y(zw(k - 1), eyr(i, j)) * v_w1
                  + J(i, j) * w(i, j, k - 1);
      auto alpha = 1 / dzw(k - 1);
      div_u(i, j, k) += (tmp_k - tmp_k1) * alpha;
    });
  stream.fence();

  return div_u;
}

template<typename MT>
std::pair<Real, std::array<int, 3>>
max_div_impl(const Vector3Field<Real***>& vec_u,
             const MT&                    mesh,
             CenterPt /*unused*/)
{
  Kokkos::Profiling::ScopedRegion region("Max div");

  const auto& u         = vec_u.x;
  const auto& v         = vec_u.y;
  const auto& w         = vec_u.z;
  const auto  is_top    = mesh.grid.comm().is_last(2);
  const auto  is_bottom = mesh.grid.comm().is_first(2);
  const auto  nx        = local_extent(u, 0);
  const auto  ny        = local_extent(u, 1);
  const auto  nz        = local_extent(u, 2);
  const auto& dz        = mesh.dz;
  const auto& zw        = mesh.zw;
  const auto& dzw       = mesh.dzw;
  const auto& eta_x     = mesh.ex;
  const auto& eta_y     = mesh.ey;
  const auto& invJ      = mesh.invJ;

  MDView<Real***, default_memory_pool> div_u(
    Kokkos::view_alloc("div_u", Kokkos::WithoutInitializing),
    create_local_layout(mesh.grid.pencil));
  MDView<Real***, default_memory_pool> tmp(
    Kokkos::view_alloc("div tmp", Kokkos::WithoutInitializing), div_u.layout());

  using member_t = GridPolicy<>::member_type;
  using Kokkos::TeamVectorRange;

  auto      stream  = get_next_stream();
  const int z_begin = is_bottom ? 1 : 0;
  const int z_end   = is_top ? nz - 1 : nz;
  const int z_size  = z_end - z_begin;
  // calculate d(J^{-1}u)/dξ and d(J^{-1}v)/d𝜓
  auto policy1 = GridPolicy<>(stream, ny * z_size, Kokkos::AUTO());
  parallel_for(
    "(u,v)/J", policy1, KOKKOS_LAMBDA(member_t team) {
      int k = team.league_rank() / ny + z_begin;
      int j = team.league_rank() % ny;
      parallel_for(
        TeamVectorRange(team, nx), KOKKOS_TR_LAMBDA(int& i) {
          div_u(i, j, k) = u(i, j, k) * invJ(i, j);
          tmp(i, j, k)   = v(i, j, k) * invJ(i, j);
        });
    });
  spectral::ddx(div_u, div_u, mesh.grid, stream);
  spectral::ddy_and_add(div_u, tmp, mesh.grid, stream);

  // calculate d(J^{-1}W)/dζ
  if (!check_same_layout_and_offset(u, v, w)) {
    throw std::runtime_error("Mismatched extents in vec_u");
  }
  Kokkos::parallel_for(
    "dW/dz",
    GridPolicy<>(stream, ny * z_size, Kokkos::AUTO()),
    KOKKOS_LAMBDA(member_t team) {
      int  j      = team.league_rank() / z_size;
      int  k      = team.league_rank() % z_size + z_begin;
      auto alpha  = 1 / dzw(k - 1);
      auto beta   = dzw(k - 1) / 2 / dz(k);
      auto beta1  = dzw(k - 1) / 2 / dz(k - 1);
      auto offset = &u(0, j, k) - u.data();
      auto stride = static_cast<int>(u.stride(2));
      Kokkos::parallel_for(
        Kokkos::TeamThreadRange(team, nx), KOKKOS_TR_LAMBDA(int& i) {
          auto u_w =
            itp2node(u.data()[i + offset], u.data()[i + offset + stride], beta);
          auto v_w =
            itp2node(v.data()[i + offset], v.data()[i + offset + stride], beta);
          auto tmp_k = u_w * MT::invJ_zeta_x(zw(k), eta_x(i, j))
                     + v_w * MT::invJ_zeta_y(zw(k), eta_y(i, j))
                     + w.data()[i + offset] * MT::invJ_zeta_z();
          auto u_w1 = itp2node(
            u.data()[i + offset], u.data()[i + offset - stride], beta1);
          auto v_w1 = itp2node(
            v.data()[i + offset], v.data()[i + offset - stride], beta1);
          auto tmp_k1 = u_w1 * MT::invJ_zeta_x(zw(k - 1), eta_x(i, j))
                      + v_w1 * MT::invJ_zeta_y(zw(k - 1), eta_y(i, j))
                      + w.data()[i + offset - stride] * MT::invJ_zeta_z();
          div_u(i, j, k) += (tmp_k - tmp_k1) * alpha;
        });
    });

  spectral::dealias(div_u, mesh.grid, stream);

  auto valid_div_u =
    Kokkos::subview(div_u, Kokkos::ALL, Kokkos::ALL, std::pair(z_begin, z_end));
  auto div_flat           = Kokkos::View<Real*>(valid_div_u.data(),
                                      (z_end - z_begin) * div_u.stride(2));
  using max_div_reducer_t = Kokkos::MaxLoc<Real, int>;
  max_div_reducer_t::value_type maxDivAndLoc{};
  Kokkos::parallel_reduce(
    "max div reduction",
    LoopPolicy<1>(stream, 0, div_flat.extent(0)),
    KOKKOS_LAMBDA(int i, max_div_reducer_t::value_type& l_max) {
      const auto t = Kokkos::abs(div_flat(i));
      if (t > l_max.val) {
        l_max.val = t;
        l_max.loc = i;
      }
    },
    max_div_reducer_t(maxDivAndLoc));

  stream.fence();
  const std::array<int, 3> loc_ijk = [&](int loc) {
    auto k = loc / (int)div_u.stride(2);
    loc -= k * (int)div_u.stride(2);
    auto j = loc / (int)div_u.stride(1);
    auto i = loc % (int)div_u.stride(1);
    return std::array<int, 3>{i, j, k + z_begin};
  }(maxDivAndLoc.loc);
  return {maxDivAndLoc.val, loc_ijk};
}
} // namespace detail

MDView<Real***> div(const Vector3Field<Real***>& vec_invJu,
                    const BottomWaveMesh&        mesh,
                    CenterPt,
                    std::string_view label)
{
  auto const region = Kokkos::Profiling::ScopedRegion(std::string(label));
  return detail::div_impl(vec_invJu, mesh, CenterPt{}, label);
}

MDView<Real***> div(const Vector3Field<Real***>& vec_invJu,
                    const TopWaveMesh&           mesh,
                    CenterPt,
                    std::string_view label)
{
  auto const region = Kokkos::Profiling::ScopedRegion(std::string(label));
  return detail::div_impl(vec_invJu, mesh, CenterPt{}, label);
}

std::pair<Real, std::array<int, 3>> max_div(const Vector3Field<Real***>& vec_u,
                                            const BottomWaveMesh&        mesh,
                                            CenterPt)
{
  return detail::max_div_impl(vec_u, mesh, CenterPt{});
}

std::pair<Real, std::array<int, 3>>
max_div(const Vector3Field<Real***>& vec_u, const TopWaveMesh& mesh, CenterPt)
{
  return detail::max_div_impl(vec_u, mesh, CenterPt{});
}
} // namespace alps::solver
