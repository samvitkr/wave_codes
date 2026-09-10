#include "diffusion_cn.h"

#include <common/kokkos_abstraction/exec_policy.h>
#include <common/kokkos_abstraction/pool_space.h>
#include <solvers/mesh/mesh.h>
#include <solvers/poisson/tridiagonal_solver.h>
#include <spectral/spectral.h>

#include <Kokkos_Core.hpp>

namespace alps::solver {

DiffusionCNEqn<CenterPt>::DiffusionCNEqn(Mesh const&       mesh,
                                         double            alph,
                                         DiffusionCNBCType top_bc_type,
                                         DiffusionCNBCType bottom_bc_type)
  : mesh_{mesh}
  , coeff_d(MDView<Real***, default_memory_pool>(
      Kokkos::view_alloc("diff_cn_coeff_d", Kokkos::WithoutInitializing),
      mesh.grid
        .get_r2c_xy_output_layout<Real, Kokkos::DefaultExecutionSpace>()))
  , coeff_l(MDView<Real***, default_memory_pool>(
      Kokkos::view_alloc("diff_cn_coeff_l", Kokkos::WithoutInitializing),
      coeff_d.layout()))
  , coeff_u(MDView<Real***, default_memory_pool>(
      Kokkos::view_alloc("diff_cn_coeff_u", Kokkos::WithoutInitializing),
      coeff_d.layout()))
  , solver{create_tridiagonal_solver<Real>(coeff_d.extent_int(0),
                                           coeff_d.extent_int(1),
                                           coeff_d.extent_int(2),
                                           mesh.grid.comm().axis_comm[2])}
  , alpha{alph}
  , top_bc{top_bc_type}
  , bottom_bc{bottom_bc_type}
{}

void DiffusionCNEqn<CenterPt>::initialize()
{
  set_cn_operator_coefficients(coeff_d,
                               coeff_l,
                               coeff_u,
                               alpha,
                               top_bc,
                               bottom_bc,
                               mesh_,
                               CenterPt{},
                               Kokkos::DefaultExecutionSpace());
  Kokkos::DefaultExecutionSpace().fence();
  solver->setup(coeff_d, coeff_l, coeff_u);
}

void DiffusionCNEqn<CenterPt>::solve(
  MDView<Real***> const&               rhs,
  Kokkos::DefaultExecutionSpace const& stream) const
{
  if (solver->is_coeff_overwritten) {
    set_cn_operator_coefficients(coeff_d,
                                 coeff_l,
                                 coeff_u,
                                 alpha,
                                 top_bc,
                                 bottom_bc,
                                 mesh_,
                                 CenterPt{},
                                 stream);
  }
  solver->comm_.barrier();
  stream.fence();
  solver->solve(rhs, coeff_d, coeff_l, coeff_u);
}

DiffusionCNEqn<CenterPt>::~DiffusionCNEqn() = default;

DiffusionCNEqn<NodePt>::DiffusionCNEqn(Mesh const& mesh, double alph)
  : mesh_{mesh}
  , coeff_d([&] {
    auto layout =
      mesh.grid.get_r2c_xy_output_layout<Real, Kokkos::DefaultExecutionSpace>();
    // node variable has one less DOF at the top
    auto is_top = mesh.grid.comm().is_last(2);
    auto nz     = is_top ? layout.dimension[2] - 1 : layout.dimension[2];
    return MDView<Real***, default_memory_pool>(
      Kokkos::view_alloc("diff_cn_node_coeff_d", Kokkos::WithoutInitializing),
      layout.dimension[0],
      layout.dimension[1],
      nz);
  }())
  , coeff_l(MDView<Real***, default_memory_pool>(
      Kokkos::view_alloc("diff_cn_node_coeff_l", Kokkos::WithoutInitializing),
      coeff_d.layout()))
  , coeff_u(MDView<Real***, default_memory_pool>(
      Kokkos::view_alloc("diff_cn_node_coeff_u", Kokkos::WithoutInitializing),
      coeff_d.layout()))
  , solver{create_tridiagonal_solver<Real>(coeff_d.extent_int(0),
                                           coeff_d.extent_int(1),
                                           coeff_d.extent_int(2),
                                           mesh.grid.comm().axis_comm[2])}
  , alpha{alph}
{}

void DiffusionCNEqn<NodePt>::initialize()
{
  set_cn_operator_coefficients(coeff_d,
                               coeff_l,
                               coeff_u,
                               alpha,
                               top_bc,
                               bottom_bc,
                               mesh_,
                               NodePt{},
                               Kokkos::DefaultExecutionSpace());
  Kokkos::DefaultExecutionSpace().fence();
  solver->setup(coeff_d, coeff_l, coeff_u);
}

void DiffusionCNEqn<NodePt>::solve(
  MDView<Real***> const&               rhs,
  Kokkos::DefaultExecutionSpace const& stream) const
{
  if (solver->is_coeff_overwritten) {
    set_cn_operator_coefficients(coeff_d,
                                 coeff_l,
                                 coeff_u,
                                 alpha,
                                 top_bc,
                                 bottom_bc,
                                 mesh_,
                                 NodePt{},
                                 stream);
  }
  solver->comm_.barrier();
  stream.fence();
  solver->solve(rhs, coeff_d, coeff_l, coeff_u);
}

DiffusionCNEqn<NodePt>::~DiffusionCNEqn() = default;

namespace {
GridPolicy<> cn_coeffs_grid_policy(Kokkos::DefaultExecutionSpace const& stream,
                                   Kokkos::Array<int, 3> const&         ends)
{
  if constexpr (is_cuda_execution_space_v<Kokkos::DefaultExecutionSpace>
                || is_hip_execution_space_v<Kokkos::DefaultExecutionSpace>) {
    return {stream, ends[1] * ends[2], 128};
  }
  return {stream, ends[1] * ends[2], Kokkos::AUTO()};
}
} // anonymous namespace

void set_cn_operator_coefficients(MDView<Real***> const& coeff_d,
                                  MDView<Real***> const& coeff_l,
                                  MDView<Real***> const& coeff_u,
                                  double const           alpha,
                                  DiffusionCNBCType      top_bc,
                                  DiffusionCNBCType      bottom_bc,
                                  Mesh const&            mesh,
                                  CenterPt               tag,
                                  Kokkos::DefaultExecutionSpace const& stream)
{
  (void)tag;
  auto constexpr Dirichlet = DiffusionCNBCType::Dirichlet;
  auto constexpr Neumann   = DiffusionCNBCType::Neumann;

  Kokkos::Array<int, 3> const ends{
    coeff_d.extent_int(0), coeff_d.extent_int(1), coeff_d.extent_int(2)};

  auto const& dz   = mesh.dz;
  auto const& dzw  = mesh.dzw;
  auto const  hbar = mesh.hbar;

  auto const is_top    = mesh.grid.comm().is_last(2);
  auto const is_bottom = mesh.grid.comm().is_first(2);

  auto const pex      = mesh.pex;
  auto const pey      = mesh.pey;
  auto const offset_x = mesh.grid.offset(1, Pencil::Y);

  auto const policy = cn_coeffs_grid_policy(stream, ends);
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdouble-promotion"
#pragma GCC diagnostic ignored "-Wfloat-conversion"
  auto const functor = KOKKOS_LAMBDA(GridPolicy<>::member_type const& team)
  {
    auto const k = team.league_rank() / ends[1];
    auto const m = team.league_rank() % ends[1];

    auto const ax = pex * int((m + offset_x) / 2);
    if (k == 0 && is_bottom && bottom_bc == Dirichlet) {
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, ends[0]), [&](int l) {
        coeff_l(l, m, k) = 0;
        coeff_d(l, m, k) = 1;
        coeff_u(l, m, k) = 0;
      });
    } else if (k == 0 && is_bottom && bottom_bc == Neumann) {
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, ends[0]), [&](int l) {
        coeff_l(l, m, k) = 0;
        coeff_d(l, m, k) = -1 / (dz(0) * hbar);
        coeff_u(l, m, k) = 1 / (dz(0) * hbar);
      });
    } else if (k == 1 && is_bottom && bottom_bc == Dirichlet) {
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, ends[0]), [&](int l) {
        auto ay          = pey * int(l / 2);
        auto c           = -ax * ax - ay * ay;
        coeff_l(l, m, k) = -(2 + dz(1) / dz(0)) / ((dz(0) + dz(1)) * hbar)
                         / (dzw(0) * hbar) * alpha;
        coeff_d(l, m, k) =
          1 - alpha * c
          + (2 + dz(1) / dz(0)) / (dz(1) * hbar) / (dzw(0) * hbar) * alpha;
        coeff_u(l, m, k) = -(2 * dz(0) / dz(1) + 1) / ((dz(0) + dz(1)) * hbar)
                         / (dzw(0) * hbar) * alpha;
      });
    } else if (k == 1 && is_bottom && bottom_bc == Neumann) {
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, ends[0]), [&](int l) {
        auto ay          = pey * int(l / 2);
        auto c           = -ax * ax - ay * ay;
        coeff_l(l, m, k) = 0;
        coeff_d(l, m, k) =
          1 - alpha * c + alpha / (dz(1) * hbar) / (dzw(0) * hbar);
        coeff_u(l, m, k) = -alpha / (dz(1) * hbar) / (dzw(0) * hbar);
      });
    } else if (k == ends[2] - 2 && is_top && top_bc == Dirichlet) {
      auto dz0 = dz(ends[2] - 2);
      auto dz1 = dz(ends[2] - 3);
      auto dw  = dzw(ends[2] - 3);
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, ends[0]), [&](int l) {
        auto ay = pey * int(l / 2);
        auto c  = -ax * ax - ay * ay;
        coeff_l(l, m, k) =
          -(2 * dz0 / dz1 + 1) / ((dz0 + dz1) * hbar) / (dw * hbar) * alpha;
        coeff_d(l, m, k) =
          1 - alpha * c + (2 + dz1 / dz0) / (dz1 * hbar) / (dw * hbar) * alpha;
        coeff_u(l, m, k) =
          -(2 + dz1 / dz0) / ((dz0 + dz1) * hbar) / (dw * hbar) * alpha;
      });
    } else if (k == ends[2] - 2 && is_top && top_bc == Neumann) {
      auto dz1 = dz(ends[2] - 3);
      auto dw  = dzw(ends[2] - 3);
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, ends[0]), [&](int l) {
        auto ay          = pey * int(l / 2);
        auto c           = -ax * ax - ay * ay;
        coeff_l(l, m, k) = -alpha / (dz1 * hbar) / (dw * hbar);
        coeff_d(l, m, k) = 1 - alpha * c + alpha / (dz1 * hbar) / (dw * hbar);
        coeff_u(l, m, k) = 0;
      });
    } else if (k == ends[2] - 1 && is_top && top_bc == Dirichlet) {
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, ends[0]), [&](int l) {
        coeff_l(l, m, k) = 0;
        coeff_d(l, m, k) = 1;
        coeff_u(l, m, k) = 0;
      });
    } else if (k == ends[2] - 1 && is_top && top_bc == Neumann) {
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, ends[0]), [&](int l) {
        coeff_l(l, m, k) = -1 / (dz(ends[2] - 2) * hbar);
        coeff_d(l, m, k) = 1 / (dz(ends[2] - 2) * hbar);
        coeff_u(l, m, k) = 0;
      });
    } else {
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, ends[0]), [&](int l) {
        auto ay          = pey * int(l / 2);
        auto c           = -ax * ax - ay * ay;
        coeff_l(l, m, k) = -alpha / (dz(k - 1) * hbar) / (dzw(k - 1) * hbar);
        coeff_d(l, m, k) =
          1 - alpha * c
          + (1 / dz(k - 1) + 1 / dz(k)) / hbar / (dzw(k - 1) * hbar) * alpha;
        coeff_u(l, m, k) = -alpha / (dz(k) * hbar) / (dzw(k - 1) * hbar);
      });
    }
  };
  Kokkos::parallel_for("set_CN_coeff", policy, functor);
#pragma GCC diagnostic pop
}

void set_cn_operator_coefficients(MDView<Real***> const& coeff_d,
                                  MDView<Real***> const& coeff_l,
                                  MDView<Real***> const& coeff_u,
                                  double const           alpha,
                                  DiffusionCNBCType      top_bc,
                                  DiffusionCNBCType      bottom_bc,
                                  Mesh const&            mesh,
                                  NodePt                 tag,
                                  Kokkos::DefaultExecutionSpace const& stream)
{
  (void)tag;
  if (top_bc != DiffusionCNBCType::Dirichlet
      || bottom_bc != DiffusionCNBCType::Dirichlet) {
    throw std::runtime_error(
      "node-centered diffusion equation is only implemented for Dirichlet BCs");
  }

  Kokkos::Array<int, 3> ends{
    coeff_d.extent_int(0), coeff_d.extent_int(1), coeff_d.extent_int(2)};

  auto const& dz   = mesh.dz;
  auto const& dzw  = mesh.dzw;
  auto const  hbar = mesh.hbar;

  auto const is_top    = mesh.grid.comm().is_last(2);
  auto const is_bottom = mesh.grid.comm().is_first(2);

  auto const pex      = mesh.pex;
  auto const pey      = mesh.pey;
  auto const offset_x = mesh.grid.offset(1, Pencil::Y);

  auto const policy = cn_coeffs_grid_policy(stream, ends);
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdouble-promotion"
#pragma GCC diagnostic ignored "-Wfloat-conversion"
  auto const functor = KOKKOS_LAMBDA(GridPolicy<>::member_type const& team)
  {
    auto const k = team.league_rank() / ends[1];
    auto const m = team.league_rank() % ends[1];

    auto const ax = pex * int((m + offset_x) / 2);
    if (k == 0 && is_bottom) {
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, ends[0]), [&](int l) {
        coeff_l(l, m, k) = 0;
        coeff_d(l, m, k) = 1;
        coeff_u(l, m, k) = 0;
      });
    } else if (k >= ends[2] - 2 && is_top) {
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, ends[0]), [&](int l) {
        coeff_l(l, m, k) = 0;
        coeff_d(l, m, k) = 1;
        coeff_u(l, m, k) = 0;
      });
    } else {
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, ends[0]), [&](int l) {
        auto ay          = pey * int(l / 2);
        auto c           = -ax * ax - ay * ay;
        coeff_l(l, m, k) = -alpha / (dzw(k - 1) * hbar) / (dz(k) * hbar);
        coeff_d(l, m, k) =
          1 - alpha * c
          + alpha * (1 / dzw(k - 1) + 1 / dzw(k)) / hbar / (dz(k) * hbar);
        coeff_u(l, m, k) = -alpha / (dzw(k) * hbar) / (dz(k) * hbar);
      });
    }
  };
  Kokkos::parallel_for("set_CN_coeff", policy, functor);
#pragma GCC diagnostic pop
}

} // namespace alps::solver
