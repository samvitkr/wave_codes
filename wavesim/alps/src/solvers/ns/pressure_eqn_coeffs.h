//
// Created by xuanx004 on 6/29/24.
//

#pragma once

#include <common/container/view_types.h>
#include <common/device/device_traits.h>
#include <common/kokkos_abstraction/exec_policy.h>
#include <common/real_type.h>
#include <common/runtime/async_utils.h>
#include <solvers/ns/pressure_bctype.h>
#include <spectral/spectral.h>

namespace alps::solver {

template<typename ValueT, typename MeshT>
void set_pressure_eqn_coefficients(MDView<ValueT***> const&             d,
                                   MDView<ValueT***> const&             dl,
                                   MDView<ValueT***> const&             du,
                                   PressureBCType                       top_bc,
                                   MeshT const&                         mesh,
                                   Kokkos::DefaultExecutionSpace const& stream)
{
  using Kokkos::parallel_for;
  using Kokkos::TeamThreadRange;

  Kokkos::Array<int, 3> ends{d.extent_int(0), d.extent_int(1), d.extent_int(2)};

  auto const& dz   = mesh.dz;
  auto const& dzw  = mesh.dzw;
  auto const  hbar = mesh.hbar;

  auto const is_top    = mesh.grid.comm().is_last(2);
  auto const is_bottom = mesh.grid.comm().is_first(2);

  auto const pex      = mesh.pex;
  auto const pey      = mesh.pey;
  auto const offset_x = mesh.grid.offset(1, Pencil::Y);

  auto const policy = [&] {
    if constexpr (is_cuda_execution_space_v<Kokkos::DefaultExecutionSpace>
                  || is_hip_execution_space_v<Kokkos::DefaultExecutionSpace>) {
      return GridPolicy<>(stream, d.extent_int(1) * d.extent_int(2), 128);
    }
    return GridPolicy<>(
      stream, d.extent_int(1) * d.extent_int(2), Kokkos::AUTO());
  }();
  auto k2 = KOKKOS_LAMBDA(Real kx, Real ky, Real alpha)
  {
    return (kx * alpha) * (kx * alpha) + (ky * alpha) * (ky * alpha);
  };
  parallel_for(
    "init_peqn_coeff",
    policy,
    KOKKOS_LAMBDA(GridPolicy<>::member_type const& team) {
      auto const k = team.league_rank() / d.extent_int(1);
      auto const m = team.league_rank() % d.extent_int(1);

      auto const ax = pex * int((m + offset_x) / 2);
      if (k == 0 && is_bottom) {
        parallel_for(TeamThreadRange(team, 0, d.extent_int(0)), [&](int l) {
          auto ay     = pey * int(l / 2);
          auto c      = -k2(ax, ay, dzw(0) * hbar);
          dl(l, m, k) = 0;
          d(l, m, k)  = ValueT(-dzw(0) / dz(0));
          du(l, m, k) = ValueT(c / (2 + dz(1) / dz(0)) + dzw(0) / dz(0));
        });
      } else if (k == 1 && is_bottom) {
        parallel_for(TeamThreadRange(team, 0, d.extent_int(0)), [&](int l) {
          auto ay = pey * int(l / 2);
          auto c  = -k2(ax, ay, dzw(0) * hbar);
          dl(l, m, k) =
            ValueT((2 + dz(1) / dz(0)) * (dzw(0) / (dz(0) + dz(1))));
          d(l, m, k) = ValueT(c - (2 + dz(1) / dz(0)) * (dzw(0) / dz(1)));
          du(l, m, k) =
            ValueT((2 * dz(0) / dz(1) + 1) * (dzw(0) / (dz(0) + dz(1))));
        });
      } else if (k == ends[2] - 2 && is_top) {
        auto dz0 = dz(ends[2] - 2);
        auto dz1 = dz(ends[2] - 3);
        auto dzk = dzw(ends[2] - 3);
        parallel_for(TeamThreadRange(team, 0, d.extent_int(0)), [&](int l) {
          auto ay     = pey * int(l / 2);
          auto c      = -k2(ax, ay, dzk * hbar);
          dl(l, m, k) = ValueT((2 * dz0 / dz1 + 1) * (dzk / (dz0 + dz1)));
          d(l, m, k)  = ValueT(c - (2 + dz1 / dz0) * (dzk / dz1));
          du(l, m, k) = ValueT((2 + dz1 / dz0) * (dzk / (dz0 + dz1)));
        });
      } else if (k == ends[2] - 1 && is_top) {
        parallel_for(TeamThreadRange(team, 0, d.extent_int(0)), [&](int l) {
          auto ay  = pey * int(l / 2);
          auto dzk = dzw(ends[2] - 3);
          auto c   = -k2(ax, ay, dzk * hbar);
          if (c < 0) {
            const auto nz = ends[2];
            dl(l, m, k) =
              ValueT(-c / (2 + dz(nz - 3) / dz(nz - 2)) - dzk / dz(nz - 2));
            d(l, m, k)  = ValueT(dzk / dz(nz - 2));
            du(l, m, k) = 0;
          } else {
            // wavenumber 0 is Dirichlet to fix the gauge
            dl(l, m, k) = 0;
            d(l, m, k)  = -1;
            du(l, m, k) = 0;
          }
        });
      } else {
        parallel_for(TeamThreadRange(team, 0, d.extent_int(0)), [&](int l) {
          auto ay     = pey * int(l / 2);
          auto c      = -k2(ax, ay, dzw(k - 1) * hbar);
          dl(l, m, k) = ValueT(dzw(k - 1) / dz(k - 1));
          d(l, m, k) =
            ValueT(c - (dzw(k - 1) / dz(k - 1) + dzw(k - 1) / dz(k)));
          du(l, m, k) = ValueT(dzw(k - 1) / dz(k));
        });
      }
    });
  if (top_bc == PressureBCType::DIRICHLET && is_top) {
    // the last equation only has diagonal term
    Kokkos::deep_copy(
      stream, subview(dl, Kokkos::ALL, Kokkos::ALL, ends[2] - 1), 0);
    Kokkos::deep_copy(
      stream, subview(d, Kokkos::ALL, Kokkos::ALL, ends[2] - 1), 1);
  }
}

} // namespace alps::solver
