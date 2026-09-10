//
// Created by xuananqing on 3/7/23.
//

#include "dynamic_smagorinsky.h"

#include "detail/dynamic_smagorinsky_scalar_functors.h"
#include "detail/nut_report.h"
#include "detail/smagorinsky_report.h"

#include <common/base/macros.h>
#include <common/kokkos_abstraction/exec_policy.h>
#include <common/kokkos_abstraction/pool_space.h>
#include <common/runtime/async_utils.h>
#include <decomp/ghost_cell_exchange.h>
#include <spectral/spectral.h>

#include <Kokkos_Profiling_ScopedRegion.hpp>
#include <mpipp/collectives.h>

namespace alps::solver {

DynamicSmagorinskyScalar::DynamicSmagorinskyScalar(
  DynamicSmagorinsky& sgs_model)
  : sgs_model_(&sgs_model)
{
  sgs_model.keep_Smag = true;
}

DynamicSmagorinskyScalarOptions
DynamicSmagorinskyScalarOptions::parse_from(ConfigTable const& /*config*/)
{
  return DynamicSmagorinskyScalarOptions{};
}

// calculate Sc_{sgs}^{-1} C_s^2
void calculate_invScCs_DS(DynamicSmagorinskyScalar const& scalar_sgs_model,
                          HaloView<Real const***> const&  f,
                          Vector3Field<Real***> const&    grad_f,
                          FlowField const&                flow)
{
  using LM_floating_t = DynamicSmagorinsky::LM_floating_t;
  using Kokkos::parallel_for;
  constexpr auto PI = Kokkos::numbers::pi_v<double>;

  const auto& sgs_model = *scalar_sgs_model.sgs_model_;

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
    if constexpr (is_hip_execution_space_v<Device>) return {64, 1, 4};
    return {0, 0, 0};
  }();

  if (scalar_sgs_model.C0_Delta2.span() < 1u) {
    scalar_sgs_model.C0_Delta2 = MDView<Real*, default_memory_pool>(
      Kokkos::view_alloc("invSc*C0*Delta^2", Kokkos::WithoutInitializing), nz);
  }

  // Obtain uᵢuⱼ and test-filtered velocity
  const Vector3Field<LM_floating_t***, default_memory_pool> Ki(
    Kokkos::view_alloc("Ki", Kokkos::WithoutInitializing),
    mesh.extents(),
    {0, 0, 0}); // Ki will be computed at cell centers

  const auto z_pair = std::pair{z_begin, z_end};

  auto reqs = async_update_halo_lower_z(grid.partition(), grad_f.z, 1);

  {
    MDView<LM_floating_t***, default_memory_pool> f_test(
      "f test filter", nx, ny, nz);
    const Vector3Field<LM_floating_t***, default_memory_pool> u_test(
      Kokkos::view_alloc("u test filter", Kokkos::WithoutInitializing),
      local_extents(flow.u));
    const auto calc_Ki =
      detail::Ki_functor<LM_floating_t>(f, flow.u, f_test, u_test, Ki);

    calc_Ki.execute(LoopPolicy<3, detail::Ki_Step1_tag>(
      stream1, {0, 0, z_begin}, {nx, ny, z_end}, tile_size));

    test_filter(subview(f_test, ALL, ALL, z_pair), stream1);
    test_filter(subview(u_test.x, ALL, ALL, z_pair).view(), stream1);
    test_filter(subview(u_test.y, ALL, ALL, z_pair).view(), stream1);
    test_filter(subview(u_test.z, ALL, ALL, z_pair).view(), stream1);
    test_filter(subview(Ki.x, ALL, ALL, z_pair).view(), stream1);
    test_filter(subview(Ki.y, ALL, ALL, z_pair).view(), stream1);
    test_filter(subview(Ki.z, ALL, ALL, z_pair).view(), stream1);

    // calculate Lij = filter(uᵢuⱼ) - filter(uᵢ)filter(uⱼ)
    calc_Ki.execute(LoopPolicy<3, detail::Ki_Step2_tag>(
      stream1, {0, 0, z_begin}, {nx, ny, z_end}, tile_size));

    reqs.waitall(); // wait for grad_f halo exchange
    stream1.fence();
  } // lifetime of u_test ends here

  // Sij at test-filter scale and Mij will be computed at the cell centers
  const Vector3Field<LM_floating_t***, default_memory_pool> Xi(
    Kokkos::view_alloc("Xi", Kokkos::WithoutInitializing),
    mesh.extents(),
    {0, 0, 0});
  {
    const Vector3Field<LM_floating_t***, default_memory_pool> grad_f_test(
      Kokkos::view_alloc("grad(f) test", Kokkos::WithoutInitializing),
      local_extents(grad_f));
    const auto calc_Xi =
      detail::Xi_functor<LM_floating_t>(grad_f,
                                        grad_f_test,
                                        Xi,
                                        sgs_model.Smag,
                                        sgs_model.Smag_test,
                                        (LM_floating_t)test_ratio);

    // Calculate |S|Sᵢⱼ and store |S|
    calc_Xi.execute(LoopPolicy<3, detail::Xi_Step1_tag>(
      stream1, {0, 0, z_begin}, {nx, ny, z_end}, tile_size));

    test_filter(subview(grad_f_test.x, ALL, ALL, z_pair).view(), stream1);
    test_filter(subview(grad_f_test.y, ALL, ALL, z_pair).view(), stream1);
    test_filter(subview(grad_f_test.z, ALL, ALL, z_pair).view(), stream1);
    test_filter(subview(Xi.x, ALL, ALL, z_pair).view(), stream1);
    test_filter(subview(Xi.y, ALL, ALL, z_pair).view(), stream1);
    test_filter(subview(Xi.z, ALL, ALL, z_pair).view(), stream1);
    // test_filter(subview(S_mag, ALL, ALL, z_pair).view(), stream1);

    calc_Xi.execute(LoopPolicy<3, detail::Xi_Step2_tag>(
      stream1, {0, 0, z_begin}, {nx, ny, z_end}, tile_size));

    stream1.fence();
  } // lifetime of Sij_test ends here

  // Calculate KᵢXᵢ and XᵢXᵢ and their local plane sum reduction
  MDView<Real* [2], default_memory_pool> const local_sum("local sum", nz);

  const auto calc_KiXi =
    detail::KiXi_functor<LM_floating_t>(Ki, Xi, local_sum, z_begin, z_end);
  calc_KiXi.execute(stream1);

  // perform global reduction on the host
  MDView<Real* [2], Kokkos::HostSpace> const global_sum("global sum", nz);
  {
    const auto local_sum_h =
      Kokkos::create_mirror_view(Kokkos::HostSpace(), local_sum);
    Kokkos::deep_copy(stream1, local_sum_h, local_sum);
    stream1.fence();
    mpipp::allreduce(nonstd::span(local_sum_h.data(), long(nz) * 2),
                     nonstd::span(global_sum.data(), long(nz) * 2),
                     mpipp::plus<Real>(),
                     grid.comm().axis_comm[1]);
  }
  for (int k = 0; k < nz; ++k) { // (1/Scₜ)CₛΔ²=KᵢXᵢ/(XᵢXᵢ)
    global_sum(k, 0) = global_sum(k, 0) / global_sum(k, 1);
    if (Kokkos::isnan(global_sum(k, 0)) || global_sum(k, 0) < 0) {
      global_sum(k, 0) = 0;
    }
  }

  Kokkos::deep_copy(
    stream1, scalar_sgs_model.C0_Delta2, subview(global_sum, ALL, 0));

  stream1.fence();

  if (scalar_sgs_model.C0_Delta2.is_allocated()
      && sgs_model.options_.enable_debug_C0) {
    auto const logger_name = "scalar_smagorinsky_Cs_file_logger_" + f.label();
    RotatingFileSinkConfig file_config{"logs/scalar_smagorinsky_Cs_" + f.label()
                                       + ".log"};
    detail::report_dynamic_smagorinsky_C0(
      scalar_sgs_model, flow.time, mesh, logger_name, file_config);
  }
}

void calculate_eddy_diffusivity(
  HaloView<Real***> const&        nuD,
  DynamicSmagorinskyScalar const& scalar_sgs_model,
  ScalarField const&              scalar,
  Vector3Field<Real***> const&    grad_f,
  FlowField const&                flow,
  bool                            skip_update_C0)
{
  using Kokkos::parallel_for;
  using member_t = GridPolicy<>::member_type;

  auto const& f = scalar.array;
  auto const  region =
    Kokkos::Profiling::ScopedRegion("SGS diffusivity " + scalar.label());

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

  const auto& sgs_model = *scalar_sgs_model.sgs_model_;
  if (!sgs_model.Smag.is_allocated()) {
    throw std::runtime_error("Velocity Smagorinsky tensor not computed");
  }

  if (!skip_update_C0 || !scalar_sgs_model.C0_Delta2.is_allocated()) {
    if (!sgs_model.Smag_test.is_allocated()) {
      throw std::runtime_error("Velocity Smagorinsky tensor not computed");
    }
    calculate_invScCs_DS(scalar_sgs_model, f, grad_f, flow);
  }

  const auto& C0_Delta2 = scalar_sgs_model.C0_Delta2;
  const auto& Smag      = scalar_sgs_model.sgs_model_->Smag;
  parallel_for(
    "nu_sgs const",
    GridPolicy<>(stream1, ny * (z_end - z_begin), Kokkos::AUTO()),
    KOKKOS_LAMBDA(member_t team) {
      const int k = team.league_rank() / ny + z_begin;
      const int j = team.league_rank() % ny;
      parallel_for(
        Kokkos::TeamThreadRange(team, nx), KOKKOS_TR_LAMBDA(int& i) {
          nuD(i, j, k) = Smag(i, j, k) * C0_Delta2(k);
        });
    });
  if (is_top) {
    Kokkos::deep_copy(stream2, subview(nuD, ALL, ALL, z_end).view(), 0);
  }
  if (is_bottom) {
    Kokkos::deep_copy(stream2, subview(nuD, ALL, ALL, 0).view(), 0);
  }

  stream1.fence();
  stream2.fence();

  if (scalar_sgs_model.sgs_model_->options_.enable_debug_C0) {
    std::string const logger_name =
      "eddy_diffusivity_file_logger_c" + std::to_string(scalar.unique_id);
    RotatingFileSinkConfig file_config{
      "logs/eddy_diffusivity_c" + std::to_string(scalar.unique_id) + ".log"};
    detail::report_nu_t_stats(
      create_inner_view(nuD).view(), flow.time, mesh, logger_name, file_config);
  }
}

std::string DynamicSmagorinskyScalar::show() const
{
  return std::string{name};
}

} // namespace alps::solver
