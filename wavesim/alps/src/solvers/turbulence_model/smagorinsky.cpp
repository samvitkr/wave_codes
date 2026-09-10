//
// Created by xuananqing on 3/7/23.
//

#include "smagorinsky.h"

#include "detail/smagorinsky_report.h"

#include <common/base/macros.h>
#include <common/kokkos_abstraction/exec_policy.h>
#include <common/math.h>
#include <common/runtime/async_utils.h>
#include <decomp/ghost_cell_exchange.h>
#include <spectral/spectral.h>

#include <Kokkos_Profiling_ScopedRegion.hpp>
#include <fmt/format.h>

namespace alps::solver {

ConstantSmagorinskyOptions
ConstantSmagorinskyOptions::parse_from(ConfigTable const& config)
{
  ConstantSmagorinskyOptions options{};

  options.C0      = config.get_value_or("Cs0", (double)options.C0);
  options.kappa   = config.get_value_or("kappa", (double)options.kappa);
  options.wall_z0 = config.get_value_or("wall_z0", (double)options.wall_z0);
  options.wall_damp_exp =
    config.get_value_or("wall_damp_exp", (int)options.wall_damp_exp);
  options.horizontal_filter_width_scale =
    config.get_value_or("horizontal_filter_width_scale",
                        (double)options.horizontal_filter_width_scale);
  options.enable_debug_C0 = config.get_value_or("debug_C0", false);

  if (options.horizontal_filter_width_scale < 1) {
    throw std::runtime_error(
      fmt::format("horizontal_filter_width_scale must be >= 1, got {}",
                  options.horizontal_filter_width_scale));
  }

  return options;
}

ConstantSmagorinsky::ConstantSmagorinsky(ConstantSmagorinskyOptions options)
  : options_{std::move(options)}
{}

std::string ConstantSmagorinsky::show() const
{
  return fmt::format("{} (C0 = {}, kappa = {}, z0 = {}, wall_damp_exp = {})",
                     name,
                     options_.C0,
                     options_.kappa,
                     options_.wall_z0,
                     options_.wall_damp_exp);
}

/// Calculate a SGS mixing length that decreases close to wall
/* See Bou-Zeid, Meneveau & Parlange 2005 (PoF) */
void calculate_scaled_delta_Parlange(
  const MDView<Real*, Kokkos::HostSpace>& scaled_delta,
  const ConstantSmagorinsky&              sgs_model,
  const FlowField&                        flow)
{
  constexpr auto SHEAR_FREE_THRESHOLD{std::numeric_limits<Real>::epsilon()};
  constexpr auto PI = Kokkos::numbers::pi_v<double>;

  const Mesh& mesh = flow.mesh;

  // The filter width considers the effective resolution from dealias
  const auto dxf =
    Real(PI / (double)mesh.pex / int(mesh.grid.global_extent(0) / 2)
         * sgs_model.options_.horizontal_filter_width_scale);
  const auto dyf =
    Real(PI / (double)mesh.pey / int(mesh.grid.global_extent(1) / 2)
         * sgs_model.options_.horizontal_filter_width_scale);
  const auto C0 = Real(sgs_model.options_.C0);

  // Calculate the Smagorinsky coefficient on the host
  const auto determine_stress_free = [](const auto& bc) {
    if (auto const* ptr_bc = dynamic_cast<GradientWall const*>(bc.get());
        ptr_bc != nullptr) {
      return static_cast<bool>(Kokkos::hypot(ptr_bc->grad_1, +ptr_bc->grad_2)
                               < SHEAR_FREE_THRESHOLD);
    }
    if (auto const* ptr_bc =
          dynamic_cast<TangentialStressWall const*>(bc.get());
        ptr_bc != nullptr) {
      return static_cast<bool>(Kokkos::hypot(ptr_bc->tau_1, ptr_bc->tau_2)
                               < SHEAR_FREE_THRESHOLD);
    }
    return false;
  };
  bool const top_bc_is_stress_free{determine_stress_free(flow.top_bc)};
  bool const bottom_bc_is_stress_free{determine_stress_free(flow.bottom_bc)};

  // calculate a SGS mixing length that decreases close to wall
  // See Bou-Zeid, Meneveau & Parlange 2005 (PoF)
  if (top_bc_is_stress_free && bottom_bc_is_stress_free) {
    // both top and bottom stress free
    for (auto k = 0; k < scaled_delta.extent_int(0); ++k) {
      const auto dzf  = mesh.dzw_h(k - 1) * mesh.hbar;
      scaled_delta(k) = Kokkos::cbrt(dxf * dyf * dzf);
    }
  } else {
    if (top_bc_is_stress_free) {
      // top stress free, bottom wall
      for (auto k = 0; k < scaled_delta.extent_int(0); ++k) {
        scaled_delta(k) = mesh.zz_h(k) * mesh.hbar; // distance to bottom
      }
    } else if (bottom_bc_is_stress_free) {
      // bottom stress free, top wall
      for (auto k = 0; k < scaled_delta.extent_int(0); ++k) {
        scaled_delta(k) = (1 - mesh.zz_h(k)) * mesh.hbar; // distance to top
      }
    } else {
      // top and bottom both walls
      for (auto k = 0; k < scaled_delta.extent_int(0); ++k) {
        const auto z = Kokkos::min(mesh.zz_h(k), 1 - mesh.zz_h(k))
                     * mesh.hbar; // distance to nearest wall
        scaled_delta(k) = z;
      }
    }

    for (auto k = 0; k < scaled_delta.extent_int(0); ++k) {
      const auto dzf   = mesh.dzw_h(k - 1) * mesh.hbar;
      const auto delta = Kokkos::cbrt(dxf * dyf * dzf);
      const auto w_exp = sgs_model.options_.wall_damp_exp;
      scaled_delta(k)  = static_cast<Real>(Kokkos::pow(
        Kokkos::pow(C0, w_exp)
            * Kokkos::pow(sgs_model.options_.kappa * (double)scaled_delta(k)
                            + sgs_model.options_.wall_z0,
                          -w_exp)
          + Kokkos::pow(delta, -w_exp),
        -1.0 / w_exp));
    }
  }
}

void calculate_eddy_viscosity(HaloView<Real***> const&          nu,
                              ConstantSmagorinsky const&        sgs_model,
                              SymmTensor33Field<Real***> const& Sij,
                              FlowField const&                  flow)
{
  using Kokkos::parallel_for;
  using member_t = GridPolicy<>::member_type;

  auto const region = Kokkos::Profiling::ScopedRegion("SGS viscosity");

  const Mesh& mesh = flow.mesh;
  const Grid& grid = mesh.grid;
  const auto  nx   = grid.extent(0);
  const auto  ny   = grid.extent(1);
  const auto  nz   = grid.extent(2);

  const auto C0 = Real(sgs_model.options_.C0);

  // calculate a SGS mixing length that decreases close to wall
  // See Bou-Zeid, Meneveau & Parlange 2005 (PoF)
  MDView<Real*, Kokkos::HostSpace> const scaled_delta("delta", nz);
  calculate_scaled_delta_Parlange(scaled_delta, sgs_model, flow);

  const auto stream1 = get_next_stream();
  const auto stream2 = get_next_stream();

  const auto s_delta =
    Kokkos::create_mirror_view(default_memory_pool(), scaled_delta);
  Kokkos::deep_copy(stream1, s_delta, scaled_delta);

  // Update the ghost cells of S13 and S23 before using them
  auto reqs = async_update_halo_lower_z(grid.partition(), Sij.xz, 1);
  reqs.push(async_update_halo_lower_z(grid.partition(), Sij.yz, 2));
  reqs.waitall();

  using alps::square;
  const auto z_begin = grid.comm().is_first(2) ? 1 : 0;
  const auto z_end   = grid.comm().is_last(2) ? nz - 1 : nz;
  parallel_for(
    "nu_sgs const",
    GridPolicy<>(stream1, ny * (z_end - z_begin), Kokkos::AUTO()),
    KOKKOS_LAMBDA(member_t team) {
      const int k = team.league_rank() / ny + z_begin;
      const int j = team.league_rank() % ny;
      parallel_for(
        Kokkos::TeamThreadRange(team, nx), KOKKOS_TR_LAMBDA(int& i) {
          auto S2 = square(Sij.xx(i, j, k)) + square(Sij.yy(i, j, k))
                  + square(Sij.zz(i, j, k)) + 2 * square(Sij.xy(i, j, k));
          const auto S13 = itp2center(Sij.xz(i, j, k - 1), Sij.xz(i, j, k));
          const auto S23 = itp2center(Sij.yz(i, j, k - 1), Sij.yz(i, j, k));
          S2 += 2 * (square(S13) + square(S23));

          nu(i, j, k) = Kokkos::sqrt(2 * S2) * square(C0 * s_delta(k));
        });
    });
  if (grid.comm().is_first(2)) {
    // bottom processor
    Kokkos::deep_copy(stream2, subview(nu, ALL, ALL, 0).view(), 0);
  }
  if (grid.comm().is_last(2)) {
    Kokkos::deep_copy(stream2, subview(nu, ALL, ALL, z_end).view(), 0);
  }

  stream1.fence();
  stream2.fence();

  if (sgs_model.options_.enable_debug_C0) {
    std::string            logger_name = "smagorinsky_Cs_file_logger";
    RotatingFileSinkConfig file_config{"logs/smagorinsky_Cs.log"};
    detail::report_smagorinsky_C0(
      sgs_model, scaled_delta, flow.time, mesh, logger_name, file_config);
  }
}

} // namespace alps::solver
