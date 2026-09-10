#include "tdma.h"

#include <common/base/macros.h>
#include <common/container/view_utils.h>
#include <common/device/device_traits.h>
#include <common/kokkos_abstraction/exec_policy.h>
#include <common/kokkos_abstraction/pool_space.h>
#include <solvers/mesh/curvilinear_mesh.h>
#include <solvers/ns/curvilinear_common/diffusion_cn_zeta_coeff.h>
#include <spectral/spectral.h>

#include <Kokkos_Core.hpp>
#include <mpipp/point2point.h>

namespace alps::solver {
namespace {
using Kokkos::ALL;
using Kokkos::subview;

template<class View>
void recv_view(View const&                          view,
               int                                  source,
               int                                  tag,
               mpipp::communicator const&           comm,
               Kokkos::DefaultExecutionSpace const& stream)
{
  if constexpr (!alps::mpi_can_access_v<default_memory_pool>) {
    auto host = Kokkos::create_mirror_view(
      Kokkos::WithoutInitializing,
      alps::PoolSpace<Kokkos::SharedHostPinnedSpace>(),
      view);
    mpipp::recv(nonstd::span(host.data(), host.span()),
                source,
                tag,
                comm,
                mpipp::status_ignore);
    Kokkos::deep_copy(stream, view, host);
    stream.fence(); // ensure host view is available before return
  } else {
    mpipp::recv(nonstd::span(view.data(), view.span()),
                source,
                tag,
                comm,
                mpipp::status_ignore);
  }
}

template<class View>
void send_view(View const&                          view,
               int                                  dest,
               int                                  tag,
               mpipp::communicator const&           comm,
               Kokkos::DefaultExecutionSpace const& stream)
{
  if constexpr (!alps::mpi_can_access_v<default_memory_pool>) {
    auto host = Kokkos::create_mirror_view(
      Kokkos::WithoutInitializing,
      alps::PoolSpace<Kokkos::SharedHostPinnedSpace>(),
      view);
    Kokkos::deep_copy(stream, host, view);
    stream.fence();
    mpipp::send(nonstd::span(host.data(), host.span()), dest, tag, comm);
  } else {
    stream.fence();
    mpipp::send(nonstd::span(view.data(), view.span()), dest, tag, comm);
  }
}

template<typename CoeffFunctor>
struct ForwardFunctor
{
  struct BottomBoundary
  {};

  MDView<Real***>       dst;
  MDView<Real const***> rhs;
  MDView<Real***>       du;
  MDView<Real***>       recv_buf;
  CoeffFunctor          coeff;

  ForwardFunctor(MDView<Real***> const&       dst_,
                 MDView<Real const***> const& rhs_,
                 MDView<Real***> const&       du_,
                 MDView<Real***> const&       recv_buf_,
                 CoeffFunctor                 coeff_)
    : dst{dst_}
    , rhs{rhs_}
    , du{du_}
    , recv_buf{recv_buf_}
    , coeff{std::move(coeff_)}
  {}

  KOKKOS_FUNCTION void operator()(int i, int j) const
  {
    auto c = coeff.coeff(i, j, 0);
    Real y = rhs(i, j, 0);
    c.diag -= c.lower * recv_buf(i, j, 0);
    y -= c.lower * recv_buf(i, j, 1);
    du(i, j, 0)  = c.upper / c.diag;
    dst(i, j, 0) = y / c.diag;
  }

  KOKKOS_FUNCTION void operator()(BottomBoundary /*tag*/, int i, int j) const
  {
    auto c       = coeff.coeff(i, j, 0);
    du(i, j, 0)  = c.upper / c.diag;
    dst(i, j, 0) = rhs(i, j, 0) / c.diag;
  }

  KOKKOS_FUNCTION void operator()(GridPolicy<>::member_type const& team) const
  {
    int const j = team.league_rank();
    for (int k = 1; k < du.extent_int(2); ++k) {
      Kokkos::parallel_for(
        Kokkos::TeamThreadRange(team, du.extent_int(0)), [&](int& i) {
          auto       c = coeff.coeff(i, j, k);
          Real const d = c.diag - c.lower * du(i, j, k - 1);
          du(i, j, k)  = c.upper / d;
          dst(i, j, k) = (rhs(i, j, k) - c.lower * dst(i, j, k - 1)) / d;
        });
    }
  }
};
} // namespace

template<typename CoeffFunctor>
void solve_diffusion_cn_zeta_tdma(MDView<Real***> const&               dst,
                                  MDView<Real const***> const&         rhs,
                                  CoeffFunctor const&                  coeff,
                                  int                                  nx,
                                  int                                  ny,
                                  int                                  n_eqns,
                                  mpipp::communicator const&           comm,
                                  Kokkos::DefaultExecutionSpace const& stream)
{
  auto const rank      = comm.rank();
  auto const is_bottom = rank == 0;
  auto const is_top    = rank == comm.size() - 1;

  MDView<Real***, default_memory_pool> du(
    Kokkos::view_alloc("cn_du", Kokkos::WithoutInitializing), nx, ny, n_eqns);
  MDView<Real***, default_memory_pool> recv_buf("cn recv", nx, ny, 2);
  MDView<Real***, default_memory_pool> send_buf("cn send", nx, ny, 2);

  ForwardFunctor<CoeffFunctor> forward_functor(dst, rhs, du, recv_buf, coeff);

  auto policy2d = LoopPolicy<2>(stream, {0, 0}, {nx, ny});
  if (!is_bottom) {
    recv_view(recv_buf, rank - 1, 0, comm, stream);
    Kokkos::parallel_for("cn_forward_first", policy2d, forward_functor);
  } else {
    Kokkos::parallel_for(
      "cn_forward_first_bottom",
      LoopPolicy<2, typename ForwardFunctor<CoeffFunctor>::BottomBoundary>(
        stream, {0, 0}, {nx, ny}),
      forward_functor);
  }
  Kokkos::parallel_for(
    "cn_forward", GridPolicy<>(stream, ny, Kokkos::AUTO()), forward_functor);

  if (!is_top) {
    Kokkos::parallel_for(
      "cn_pack_send", policy2d, KOKKOS_LAMBDA(int i, int j) {
        send_buf(i, j, 0) = du(i, j, n_eqns - 1);
        send_buf(i, j, 1) = dst(i, j, n_eqns - 1);
      });
    send_view(send_buf, rank + 1, 0, comm, stream);
  }

  if (!is_top) {
    auto value_recv = subview(recv_buf, ALL, ALL, 1);
    recv_view(value_recv, rank + 1, 1, comm, stream);
    Kokkos::parallel_for(
      "cn_back_recv_top", policy2d, KOKKOS_LAMBDA(int i, int j) {
        dst(i, j, n_eqns - 1) =
          dst(i, j, n_eqns - 1) - du(i, j, n_eqns - 1) * value_recv(i, j);
      });
  }

  Kokkos::parallel_for(
    "cn_backward",
    GridPolicy<>(stream, ny, Kokkos::AUTO()),
    KOKKOS_LAMBDA(GridPolicy<>::member_type team) {
      int const j = team.league_rank();
      for (int k = n_eqns - 2; k >= 0; --k) {
        Kokkos::parallel_for(Kokkos::TeamThreadRange(team, nx), [&](int& i) {
          dst(i, j, k) = dst(i, j, k) - du(i, j, k) * dst(i, j, k + 1);
        });
      }
    });

  if (!is_bottom) {
    auto x_bot = subview(dst, ALL, ALL, 0);
    send_view(x_bot, rank - 1, 1, comm, stream);
  }

  stream.fence();
}

template void solve_diffusion_cn_zeta_tdma<
  DiffusionCNZetaCoeffProvider<CenterPt, BottomWaveMesh>>(
  MDView<Real***> const&,
  MDView<Real const***> const&,
  DiffusionCNZetaCoeffProvider<CenterPt, BottomWaveMesh> const&,
  int,
  int,
  int,
  mpipp::communicator const&,
  Kokkos::DefaultExecutionSpace const&);
template void solve_diffusion_cn_zeta_tdma<
  DiffusionCNZetaCoeffProvider<NodePt, BottomWaveMesh>>(
  MDView<Real***> const&,
  MDView<Real const***> const&,
  DiffusionCNZetaCoeffProvider<NodePt, BottomWaveMesh> const&,
  int,
  int,
  int,
  mpipp::communicator const&,
  Kokkos::DefaultExecutionSpace const&);

} // namespace alps::solver
