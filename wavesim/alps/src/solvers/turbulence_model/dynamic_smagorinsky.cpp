//
// Created by xuananqing on 3/7/23.
//

#include "dynamic_smagorinsky.h"

#include "detail/dynamic_smagorinsky_functors.h"
#include "detail/nut_report.h"
#include "detail/smagorinsky_report.h"

#include <common/base/macros.h>
#include <common/kokkos_abstraction/exec_policy.h>
#include <common/kokkos_abstraction/pool_space.h>
#include <common/runtime/async_utils.h>
#include <decomp/ghost_cell_exchange.h>
#include <spectral/spectral.h>

#include <Kokkos_Profiling_ScopedRegion.hpp>
#include <fmt/format.h>
#include <mpipp/collectives.h>

namespace alps::solver {

DynamicSmagorinskyOptions
DynamicSmagorinskyOptions::parse_from(ConfigTable const& config)
{
  DynamicSmagorinskyOptions options{};

  options.horizontal_filter_width_scale =
    config.get_value_or("horizontal_filter_width_scale",
                        (double)options.horizontal_filter_width_scale);
  options.enable_debug_C0 =
    config.get_value_or("debug_C0", (bool)options.enable_debug_C0);

  // validate options
  if (options.horizontal_filter_width_scale < 1) {
    throw std::runtime_error(
      fmt::format("horizontal_filter_width_scale must be >= 1, got {}",
                  options.horizontal_filter_width_scale));
  }

  return options;
}

DynamicSmagorinsky::DynamicSmagorinsky(
  alps::solver::DynamicSmagorinskyOptions options)
  : options_{std::move(options)}
{}

void calculate_Cs_Germano_Lilly(HaloView<Real***> const&          S_mag,
                                DynamicSmagorinsky const&         sgs_model,
                                SymmTensor33Field<Real***> const& Sij,
                                FlowField const&                  flow,
                                bool                              keep_Smag)
{
  using LM_floating_t = DynamicSmagorinsky::LM_floating_t;
  using Kokkos::parallel_for;
  constexpr auto PI = Kokkos::numbers::pi_v<double>;

  const Mesh& mesh      = flow.mesh;
  const Grid& grid      = mesh.grid;
  const auto  nx_global = grid.global_extent(0);
  const auto  ny_global = grid.global_extent(1);
  const auto  nx        = grid.extent(0);
  const auto  ny        = grid.extent(1);
  const auto  nz        = grid.extent(2);
  const auto  is_top    = grid.comm().is_last(2);
  const auto  is_bottom = grid.comm().is_first(2);
  const auto  z_begin   = is_bottom ? 1 : 0;
  const auto  z_end     = is_top ? nz - 1 : nz;

  const auto dxf = Real(PI / (double)mesh.pex / int(nx_global / 2)
                        * sgs_model.options_.horizontal_filter_width_scale);
  const auto dyf = Real(PI / (double)mesh.pey / int(ny_global / 2)
                        * sgs_model.options_.horizontal_filter_width_scale);
  // cutoff wavenumber for test filter
  const auto cutoff_kxt = int(Kokkos::nearbyint(
    (int)(nx_global / 4) / sgs_model.options_.horizontal_filter_width_scale));
  const auto cutoff_kyt = int(Kokkos::nearbyint(
    (int)(ny_global / 4) / sgs_model.options_.horizontal_filter_width_scale));
  // test filter widths
  const auto dxt         = Real(PI / (double)mesh.pex / cutoff_kxt);
  const auto dyt         = Real(PI / (double)mesh.pey / cutoff_kyt);
  const auto test_ratio  = Kokkos::cbrt((dxt / dxf) * (dyt / dyf));
  const auto test_filter = [cutoff_kxt, cutoff_kyt, &grid](auto const& view,
                                                           auto const& stream) {
    spectral::cutoff_xy(view, cutoff_kxt, cutoff_kyt, grid, stream);
  };

  const auto     stream1   = get_next_stream();
  constexpr auto tile_size = []() -> Kokkos::Array<range_index_t, 3> {
    using Device = Kokkos::DefaultExecutionSpace;
    if constexpr (is_cuda_execution_space_v<Device>) return {16, 1, 8};
    if constexpr (is_hip_execution_space_v<Device>) return {64, 1, 2};
    return {0, 0, 0};
  }();

  if (sgs_model.C0_Delta2.span() < 1u) {
    sgs_model.C0_Delta2 = MDView<Real*, default_memory_pool>(
      Kokkos::view_alloc("C0*Delta^2", Kokkos::WithoutInitializing), nz);
  }

  // Obtain uᵢuⱼ and test-filtered velocity
  // Lij will be computed at cell centers
  // no need to allocate halo for Lij
  const SymmTensor33Field<LM_floating_t***, default_memory_pool> Lij(
    Kokkos::view_alloc("Lij", Kokkos::WithoutInitializing), local_extents(Sij));

  const auto z_pair = std::pair{z_begin, z_end};

  auto reqs = async_update_halo_lower_z(grid.partition(), Sij.xz, 1);
  reqs.push(async_update_halo_lower_z(grid.partition(), Sij.yz, 2));

  {
    const Vector3Field<LM_floating_t***, default_memory_pool> u_test(
      Kokkos::view_alloc("u test filter", Kokkos::WithoutInitializing),
      local_extents(flow.u)); // no halo for u_test
    const auto calc_Lij =
      detail::Lij_functor<LM_floating_t>(flow.u, u_test, Lij);

    // calculate Lij = uᵢuⱼ and copy u for test filter
    Kokkos::deep_copy(stream1, u_test, 0);
    Kokkos::deep_copy(stream1, Lij, 0);
    calc_Lij.execute(LoopPolicy<3, detail::Lij_Step1_tag>(
      stream1, {0, 0, z_begin}, {nx, ny, z_end}, tile_size));

    test_filter(subview(u_test.x, ALL, ALL, z_pair).view(), stream1);
    test_filter(subview(u_test.y, ALL, ALL, z_pair).view(), stream1);
    test_filter(subview(u_test.z, ALL, ALL, z_pair).view(), stream1);
    test_filter(subview(Lij.xx, ALL, ALL, z_pair).view(), stream1);
    test_filter(subview(Lij.xy, ALL, ALL, z_pair).view(), stream1);
    test_filter(subview(Lij.xz, ALL, ALL, z_pair).view(), stream1);
    test_filter(subview(Lij.yy, ALL, ALL, z_pair).view(), stream1);
    test_filter(subview(Lij.yz, ALL, ALL, z_pair).view(), stream1);
    test_filter(subview(Lij.zz, ALL, ALL, z_pair).view(), stream1);

    // calculate Lij = filter(uᵢuⱼ) - filter(uᵢ)filter(uⱼ)
    calc_Lij.execute(LoopPolicy<3, detail::Lij_Step2_tag>(
      stream1, {0, 0, z_begin}, {nx, ny, z_end}, tile_size));

    reqs.waitall(); // wait for Sij halo exchange
    stream1.fence();
  } // lifetime of u_test ends here

  // local plane sum reduction for LᵢⱼMᵢⱼ and MᵢⱼMᵢⱼ
  MDView<Real* [2], default_memory_pool> const local_sum("local sum", nz);
  {
    // Sij at test-filter scale and Mij will be computed at the cell centers
    const SymmTensor33Field<LM_floating_t***, default_memory_pool> Sij_test(
      Kokkos::view_alloc("Sij test", Kokkos::WithoutInitializing),
      local_extents(Sij));
    const SymmTensor33Field<LM_floating_t***, default_memory_pool> Mij(
      Kokkos::view_alloc("Mij", Kokkos::WithoutInitializing),
      local_extents(Sij)); // no halo for Mij

    const auto calc_Mij = detail::Mij_functor<LM_floating_t>(
      Sij, Sij_test, Mij, S_mag, (LM_floating_t)test_ratio);
    const auto calc_LijMij = detail::LijMij_functor<LM_floating_t>(
      Lij, Mij, local_sum, z_begin, z_end);

    // Calculate |S|Sᵢⱼ and store |S|
    Kokkos::deep_copy(stream1, Sij_test, 0);
    Kokkos::deep_copy(stream1, Mij, 0);
    calc_Mij.execute(LoopPolicy<3, detail::Mij_Step1_tag>(
      stream1, {0, 0, z_begin}, {nx, ny, z_end}, tile_size));

    test_filter(subview(Sij_test.xx, ALL, ALL, z_pair).view(), stream1);
    test_filter(subview(Sij_test.xy, ALL, ALL, z_pair).view(), stream1);
    test_filter(subview(Sij_test.xz, ALL, ALL, z_pair).view(), stream1);
    test_filter(subview(Sij_test.yy, ALL, ALL, z_pair).view(), stream1);
    test_filter(subview(Sij_test.yz, ALL, ALL, z_pair).view(), stream1);
    test_filter(subview(Sij_test.zz, ALL, ALL, z_pair).view(), stream1);
    test_filter(subview(Mij.xx, ALL, ALL, z_pair).view(), stream1);
    test_filter(subview(Mij.xy, ALL, ALL, z_pair).view(), stream1);
    test_filter(subview(Mij.xz, ALL, ALL, z_pair).view(), stream1);
    test_filter(subview(Mij.yy, ALL, ALL, z_pair).view(), stream1);
    test_filter(subview(Mij.yz, ALL, ALL, z_pair).view(), stream1);
    test_filter(subview(Mij.zz, ALL, ALL, z_pair).view(), stream1);
    // test_filter(subview(S_mag, ALL, ALL, z_pair).view(), stream1);

    calc_Mij.execute(LoopPolicy<3, detail::Mij_Step2_tag>(
      stream1, {0, 0, z_begin}, {nx, ny, z_end}, tile_size));
    calc_LijMij.execute(stream1);
    stream1.fence();

    if (keep_Smag) {
      sgs_model.Smag_test = std::move(Sij_test.xx);
    }
  } // lifetime of Mij, Sij_test ends here

  // perform global reduction on the host
  MDView<Real* [2], default_host_memory_pool> const global_sum("global sum",
                                                               nz);
  {
    const auto local_sum_h =
      Kokkos::create_mirror_view(default_host_memory_pool(), local_sum);
    Kokkos::deep_copy(stream1, local_sum_h, local_sum);
    stream1.fence();
    mpipp::allreduce(nonstd::span(local_sum_h.data(), nz * 2),
                     nonstd::span(global_sum.data(), nz * 2),
                     mpipp::plus<Real>(),
                     grid.comm().axis_comm[1]);
  }
  for (int k = 0; k < nz; ++k) { // CₛΔ²=LᵢⱼMᵢⱼ/(2MᵢⱼMᵢⱼ)
    global_sum(k, 0) = global_sum(k, 0) / global_sum(k, 1) / 2;
    if (Kokkos::isnan(global_sum(k, 0)) || global_sum(k, 0) < 0) {
      global_sum(k, 0) = 0;
    }
  }

  Kokkos::deep_copy(stream1, sgs_model.C0_Delta2, subview(global_sum, ALL, 0));

  stream1.fence();

  if (sgs_model.C0_Delta2.is_allocated()
      && sgs_model.options_.enable_debug_C0) {
    std::string            logger_name = "smagorinsky_Cs_file_logger";
    RotatingFileSinkConfig file_config{"logs/smagorinsky_Cs.log"};
    detail::report_dynamic_smagorinsky_C0(
      sgs_model, flow.time, mesh, logger_name, file_config);
  }
}

void calculate_eddy_viscosity(HaloView<Real***> const&          nu,
                              DynamicSmagorinsky const&         sgs_model,
                              SymmTensor33Field<Real***> const& Sij,
                              FlowField const&                  flow,
                              bool                              skip_update_C0)
{
  using LM_floating_t = DynamicSmagorinsky::LM_floating_t;
  using Kokkos::parallel_for;
  using member_t = GridPolicy<>::member_type;

  auto const region = Kokkos::Profiling::ScopedRegion("SGS viscosity");

  const Mesh& mesh      = flow.mesh;
  const Grid& grid      = mesh.grid;
  const auto  nx        = grid.extent(0);
  const auto  ny        = grid.extent(1);
  const auto  nz        = grid.extent(2);
  const auto  is_top    = grid.comm().is_last(2);
  const auto  is_bottom = grid.comm().is_first(2);
  const auto  z_begin   = is_bottom ? 1 : 0;
  const auto  z_end     = is_top ? nz - 1 : nz;

  const auto stream1 = get_next_stream();
  const auto stream2 = get_next_stream();

  if (!skip_update_C0 || !sgs_model.C0_Delta2.is_allocated()) {
    calculate_Cs_Germano_Lilly(nu, sgs_model, Sij, flow, sgs_model.keep_Smag);
  } else {
    // skip update_C0, compute Smag only
    detail::SijMag_functor<LM_floating_t> SijMag(Sij, nu);
    LoopPolicy<3> policy(stream1, {0, 0, z_begin}, {nx, ny, z_end});
    SijMag.execute(policy);
  }

  const auto& C0_Delta2 = sgs_model.C0_Delta2;
  if (!sgs_model.keep_Smag) {
    parallel_for(
      "nu_sgs const",
      GridPolicy<>(stream1, ny * (z_end - z_begin), Kokkos::AUTO()),
      KOKKOS_LAMBDA(member_t team) {
        const int k = team.league_rank() / ny + z_begin;
        const int j = team.league_rank() % ny;
        parallel_for(
          Kokkos::TeamThreadRange(team, nx), KOKKOS_TR_LAMBDA(int& i) {
            nu(i, j, k) = nu(i, j, k) * C0_Delta2(k);
          });
      });
  } else {
    if (!sgs_model.Smag.is_allocated()) {
      sgs_model.Smag = HaloView<Real***, default_memory_pool>(
        Kokkos::view_alloc("Smag", Kokkos::WithoutInitializing),
        nu.layout(),
        {0, 0, 0});
    }
    auto const& Smag = sgs_model.Smag;
    parallel_for(
      "nu_sgs const",
      GridPolicy<>(stream1, ny * (z_end - z_begin), Kokkos::AUTO()),
      KOKKOS_LAMBDA(member_t team) {
        const int k = team.league_rank() / ny + z_begin;
        const int j = team.league_rank() % ny;
        parallel_for(
          Kokkos::TeamThreadRange(team, nx), KOKKOS_TR_LAMBDA(int& i) {
            Smag(i, j, k) = nu(i, j, k);
            nu(i, j, k)   = nu(i, j, k) * C0_Delta2(k);
          });
      });
  }
  if (is_top) {
    Kokkos::deep_copy(stream2, subview(nu, ALL, ALL, z_end).view(), 0);
  }
  if (is_bottom) {
    Kokkos::deep_copy(stream2, subview(nu, ALL, ALL, 0).view(), 0);
  }

  stream1.fence();
  stream2.fence();

  if (sgs_model.options_.enable_debug_C0) {
    std::string const      logger_name = "eddy_viscosity_file_logger";
    RotatingFileSinkConfig file_config{"logs/eddy_viscosity.log"};
    detail::report_nu_t_stats(
      create_inner_view(nu).view(), flow.time, mesh, logger_name, file_config);
  }
}

std::string DynamicSmagorinsky::show() const
{
  return std::string(name);
}

void DynamicSmagorinsky::cleanup_Smag() const
{
  if (Smag.is_allocated()) {
    Smag = {};
  }
  if (Smag_test.is_allocated()) {
    Smag_test = {};
  }
}
} // namespace alps::solver
