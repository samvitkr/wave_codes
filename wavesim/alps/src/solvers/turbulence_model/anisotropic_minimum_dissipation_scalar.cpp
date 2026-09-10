#include "anisotropic_minimum_dissipation.h"

#include "detail/nut_report.h"
#include <common/container/view_utils.h>
#include <common/kokkos_abstraction/exec_policy.h>
#include <common/math.h>
#include <common/runtime/async_utils.h>
#include <decomp/ghost_cell_exchange.h>
#include <spectral/spectral.h>

#include <Kokkos_MathematicalConstants.hpp>
#include <Kokkos_Profiling_ScopedRegion.hpp>
#include <fmt/format.h>

namespace alps::solver {

namespace {
template<typename T>
KOKKOS_FORCEINLINE_FUNCTION constexpr auto term_KK(T theta_x,
                                                   T theta_y,
                                                   T theta_z,
                                                   T ux,
                                                   T uy,
                                                   T uz,
                                                   T delta_x,
                                                   T delta_y,
                                                   T delta_z,
                                                   T Cxy,
                                                   T Cz)
{
  return theta_x * ux * square(delta_x) * Cxy
       + theta_y * uy * square(delta_y) * Cxy
       + theta_z * uz * square(delta_z) * Cz;
}
} // namespace

AnisotropicMinimumDissipationScalarOptions
AnisotropicMinimumDissipationScalarOptions::parse_from(
  ConfigTable const& /*config*/)
{
  return {};
}

AnisotropicMinimumDissipationScalar::AnisotropicMinimumDissipationScalar(
  AnisotropicMinimumDissipation& sgs_model)
  : sgs_model_(&sgs_model)
{
  sgs_model.keep_grad = true;
}

/*
 * The anisotropic minimum dissipation model by Abkar, Bae and Moin (2016)
 */
void calculate_eddy_diffusivity(
  HaloView<Real***> const&                   De,
  AnisotropicMinimumDissipationScalar const& scalar_sgs_model,
  ScalarField const&                         f,
  Vector3Field<Real***> const&               grad_f,
  FlowField const&                           flow)
{
  using alps::square;
  using Kokkos::parallel_for;

  auto const region =
    Kokkos::Profiling::ScopedRegion("SGS diffusivity " + f.label());

  constexpr auto PI = Kokkos::numbers::pi_v<double>;

  auto const&                         sgs_model = *scalar_sgs_model.sgs_model_;
  Tensor33Field<Real const***> const& grad_u    = sgs_model.grad_u;
  if (!grad_u.xx.is_allocated()) {
    throw std::runtime_error("Expecting grad_u stored in sgs_model");
  }

  auto const& mesh      = flow.mesh;
  auto const& grid      = mesh.grid;
  auto const  nx_global = mesh.global_extent(0);
  auto const  ny_global = mesh.global_extent(1);
  auto const  nx        = mesh.extent(0);
  auto const  ny        = mesh.extent(1);
  auto const  nz        = mesh.extent(2);

  auto const  dx  = 2 * PI / (double)mesh.pex / nx_global;
  auto const  dy  = 2 * PI / (double)mesh.pey / ny_global;
  auto const& dzw = mesh.dzw;

  // The filter width is the effective resolution from dealiasing
  auto const dxf = Real(dx * sgs_model.options_.horizontal_filter_width_scale);
  auto const dyf = Real(dy * sgs_model.options_.horizontal_filter_width_scale);
  auto const C0  = Real(sgs_model.options_.C0);
  auto const Cz  = Real(sgs_model.options_.Cz);

  auto const stream1 = get_next_stream();
  auto const stream2 = get_next_stream();

  auto reqs = async_update_halo_lower_z(grid.partition(), grad_f.x, 1);
  reqs.push(async_update_halo_lower_z(grid.partition(), grad_f.y, 2));
  reqs.push(async_update_halo_lower_z(grid.partition(), grad_f.z, 3));

  auto constexpr tile_size = []() -> Kokkos::Array<range_index_t, 3> {
    using Device = Kokkos::DefaultExecutionSpace;
    if constexpr (is_cuda_execution_space_v<Device>) return {32, 1, 4};
    if constexpr (is_hip_execution_space_v<Device>) return {64, 1, 2};
    return {0, 0, 0};
  }();
  auto const is_top      = grid.comm().is_last(2);
  auto const is_bottom   = grid.comm().is_first(2);
  auto constexpr z_begin = 1;
  auto const z_end       = is_top ? nz - 1 : nz;

  if (!check_same_layout_and_offset(grad_f.x,
                                    grad_f.y,
                                    grad_f.z,
                                    grad_u.xx,
                                    grad_u.xy,
                                    grad_u.xz,
                                    grad_u.yx,
                                    grad_u.yy,
                                    grad_u.yz,
                                    grad_u.zx,
                                    grad_u.zy,
                                    grad_u.zz,
                                    De)) {
    // kernel below calculates address offsets for all fields instead of
    // indexing directly to reduce register pressure, thus requiring all fields
    // to have the same extents and begins
    throw std::runtime_error("Mismatched extents in grad_f and grad_u");
  }

  auto const amd_lambda = KOKKOS_LAMBDA(int i, int j, int k)
  {
    auto const dzf = dzw(k - 1);

    auto offset  = &grad_f.x(i, j, k) - grad_f.x.data();     // (i, j, k)
    auto offset1 = &grad_f.x(i, j, k - 1) - grad_f.x.data(); // (i, j, k - 1)

    // cell-centered variables
    auto const theta_x = grad_f.x.data()[offset];
    auto const theta_y = grad_f.y.data()[offset];
    // node-centered variables
    auto const theta_z =
      itp2center(grad_f.z.data()[offset1], grad_f.z.data()[offset]);

    auto ux = grad_u.xx.data()[offset];
    auto uy = grad_u.xy.data()[offset];
    auto uz = itp2center(grad_u.xz.data()[offset1], grad_u.xz.data()[offset]);
    auto uik_thetaik_thetai =
      term_KK(theta_x, theta_y, theta_z, ux, uy, uz, dxf, dyf, dzf, C0, Cz)
      * theta_x;

    ux = grad_u.yx.data()[offset];
    uy = grad_u.yy.data()[offset];
    uz = itp2center(grad_u.yz.data()[offset1], grad_u.yz.data()[offset]);
    uik_thetaik_thetai +=
      term_KK(theta_x, theta_y, theta_z, ux, uy, uz, dxf, dyf, dzf, C0, Cz)
      * theta_y;

    ux = itp2center(grad_u.zx.data()[offset1], grad_u.zx.data()[offset]);
    uy = itp2center(grad_u.zy.data()[offset1], grad_u.zy.data()[offset]);
    uz = grad_u.zz.data()[offset];
    uik_thetaik_thetai +=
      term_KK(theta_x, theta_y, theta_z, ux, uy, uz, dxf, dyf, dzf, C0, Cz)
      * theta_z;
    auto theta2_l = square(theta_x) + square(theta_y) + square(theta_z);

    auto D_sgs =
      -uik_thetaik_thetai
      / (theta2_l + Kokkos::Experimental::epsilon_v<decltype(theta2_l)>);
    De(i, j, k) = Kokkos::max(D_sgs, static_cast<decltype(D_sgs)>(0));
  };
  parallel_for(
    "De amd " + f.label(),
    LoopPolicy<3>(stream1, {0, 0, z_begin}, {nx, ny, z_end}, tile_size),
    amd_lambda);

  reqs.waitall();
  if (is_top) {
    Kokkos::deep_copy(stream2, subview(De, ALL, ALL, nz - 1).view(), 0);
  }
  if (is_bottom) {
    Kokkos::deep_copy(stream2, subview(De, ALL, ALL, 0).view(), 0);
  } else {
    parallel_for("De amd 0",
                 LoopPolicy<3>(stream2, {0, 0, 0}, {nx, ny, 1}, tile_size),
                 amd_lambda);
  }
  stream1.fence();
  stream2.fence();

  if (sgs_model.options_.enable_debug_C0) {
    std::string const logger_name =
      "eddy_diffusivity_file_logger_c" + std::to_string(f.unique_id);
    RotatingFileSinkConfig file_config{"logs/eddy_diffusivity_c"
                                       + std::to_string(f.unique_id) + ".log"};
    detail::report_nu_t_stats(
      create_inner_view(De).view(), flow.time, mesh, logger_name, file_config);
  }
}

std::string AnisotropicMinimumDissipationScalar::show() const
{
  return fmt::format("{} (Poincare constant = {}, {})",
                     name,
                     this->sgs_model_->options_.C0,
                     this->sgs_model_->options_.Cz);
}

} // namespace alps::solver
