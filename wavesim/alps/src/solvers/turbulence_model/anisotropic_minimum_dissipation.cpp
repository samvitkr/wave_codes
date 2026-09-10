//
// Created by xuanx004 on 3/8/23.
//

#include "anisotropic_minimum_dissipation.h"

#include "detail/nut_report.h"
#include <common/container/view_utils.h>
#include <common/kokkos_abstraction/exec_policy.h>
#include <common/math.h>
#include <common/runtime/async_utils.h>
#include <decomp/ghost_cell_exchange.h>
#include <solvers/operators/grad.h>
#include <spectral/spectral.h>

#include <Kokkos_MathematicalConstants.hpp>
#include <Kokkos_Profiling_ScopedRegion.hpp>
#include <fmt/format.h>

namespace alps::solver {

namespace {
template<class T>
KOKKOS_FORCEINLINE_FUNCTION constexpr auto term_KK(T du,
                                                   T dv,
                                                   T dw,
                                                   T S11,
                                                   T S22,
                                                   T S33,
                                                   T twiceS12,
                                                   T twiceS13,
                                                   T twiceS23,
                                                   T delta_k,
                                                   T C_k)
{
  return (S11 * square(du) + S22 * square(dv) + S33 * square(dw)
          + du * dv * twiceS12 + du * dw * twiceS13 + dv * dw * twiceS23)
       * delta_k * delta_k * C_k;
}
} // namespace

AnisotropicMinimumDissipationOptions
AnisotropicMinimumDissipationOptions::parse_from(ConfigTable const& config)
{
  AnisotropicMinimumDissipationOptions options{};

  options.C0 = config.get_value_or("PoincareCxy", (double)options.C0);
  options.Cz = config.get_value_or("PoincareCz", (double)options.Cz);
  options.horizontal_filter_width_scale =
    config.get_value_or("horizontal_filter_width_scale",
                        (double)options.horizontal_filter_width_scale);
  options.enable_debug_C0 =
    config.get_value_or("debug_C0", (bool)options.enable_debug_C0);

  // validate options
  if (options.C0 <= 0) {
    throw std::runtime_error(
      fmt::format("PoincareCxy must be positive, got {}", options.C0));
  }
  if (options.Cz <= 0) {
    throw std::runtime_error(
      fmt::format("PoincareCz must be positive, got {}", options.Cz));
  }
  if (options.horizontal_filter_width_scale < 1) {
    throw std::runtime_error(
      fmt::format("horizontal_filter_width_scale must be >= 1, got {}",
                  options.horizontal_filter_width_scale));
  }

  return options;
}

AnisotropicMinimumDissipation::AnisotropicMinimumDissipation(
  AnisotropicMinimumDissipationOptions options)
  : options_{std::move(options)}
{}

template<bool has_buoyancy>
struct amd_functor
{
  HaloView<Real***>       nu;
  HaloView<Real const***> uxx;
  HaloView<Real const***> uxy;
  HaloView<Real const***> uxz;
  HaloView<Real const***> uyx;
  HaloView<Real const***> uyy;
  HaloView<Real const***> uyz;
  HaloView<Real const***> uzx;
  HaloView<Real const***> uzy;
  HaloView<Real const***> uzz;

  HaloView<Real const***> bx;
  HaloView<Real const***> by;
  HaloView<Real const***> bz;

  HaloView<Real const*> dzw;
  Real                  dxf;
  Real                  dyf;
  Real                  C0;
  Real                  Cz;

  amd_functor(HaloView<Real***> const&      nu_,
              Tensor33Field<Real***> const& grad_u,
              Vector3Field<Real***> const*  grad_b,
              Mesh const&                   mesh,
              Real                          dxf_,
              Real                          dyf_,
              Real                          C0_,
              Real                          Cz_)
    : nu{nu_}
    , uxx{grad_u.xx}
    , uxy{grad_u.xy}
    , uxz{grad_u.xz}
    , uyx{grad_u.yx}
    , uyy{grad_u.yy}
    , uyz{grad_u.yz}
    , uzx{grad_u.zx}
    , uzy{grad_u.zy}
    , uzz{grad_u.zz}
    , dzw{mesh.dzw}
    , dxf{dxf_}
    , dyf{dyf_}
    , C0{C0_}
    , Cz{Cz_}
  {
    // the functor uses pointer offset to access the data, in order to reduce
    // register pressure, so we need to ensure that all fields have the same
    // layout and offset
    if (!check_same_layout_and_offset(
          uxx, uxy, uxz, uyx, uyy, uyz, uzx, uzy, uzz, nu)) {
      throw std::runtime_error("Mismatched extents in grad_u and nu");
    }

    if constexpr (has_buoyancy) {
      bx = grad_b->x;
      by = grad_b->y;
      bz = grad_b->z;

      if (!check_same_layout_and_offset(uxx, bx, by, bz)) {
        throw std::runtime_error("Mismatched extents in grad_b and grad_u");
      }
    }
  }

  KOKKOS_FUNCTION void operator()(int i, int j, int k) const
  {
    const auto dzf = dzw(k - 1);

    auto offset  = &uxx(i, j, k) - uxx.data();     // (i, j, k)
    auto offset1 = &uxx(i, j, k - 1) - uxx.data(); // (i, j, k - 1)

    // for cell-centered variables
    const auto ux = uxx.data()[offset];
    const auto vx = uyx.data()[offset];
    const auto uy = uxy.data()[offset];
    const auto vy = uyy.data()[offset];
    const auto wz = uzz.data()[offset];

    // for node-centered variables
    const auto uz = itp2center(uxz.data()[offset1], uxz.data()[offset]);
    const auto vz = itp2center(uyz.data()[offset1], uyz.data()[offset]);
    const auto wx = itp2center(uzx.data()[offset1], uzx.data()[offset]);
    const auto wy = itp2center(uzy.data()[offset1], uzy.data()[offset]);

    // trace of grad(u)
    auto tr_grad_u = square(ux) + square(vx) + square(uy) + square(vy)
                   + square(wz) + square(uz) + square(vz) + square(wx)
                   + square(wy);

    const auto twiceS12 = uy + vx;
    const auto twiceS13 = uz + wx;
    const auto twiceS23 = vz + wy;
    auto       uik_ujk_Sij =
      term_KK(ux, vx, wx, ux, vy, wz, twiceS12, twiceS13, twiceS23, dxf, C0);
    uik_ujk_Sij +=
      term_KK(uy, vy, wy, ux, vy, wz, twiceS12, twiceS13, twiceS23, dyf, C0);
    uik_ujk_Sij +=
      term_KK(uz, vz, wz, ux, vy, wz, twiceS12, twiceS13, twiceS23, dzf, Cz);

    if constexpr (has_buoyancy) {
      auto bxc = itp2center(bx.data()[offset1], bx.data()[offset]);
      uik_ujk_Sij -= bxc * wx * dxf * dxf * C0;
      auto byc = itp2center(by.data()[offset1], by.data()[offset]);
      uik_ujk_Sij -= byc * wy * dyf * dyf * C0;
      uik_ujk_Sij -= bz.data()[offset] * wz * dzf * dzf * Cz;
    }

    auto nu_sgs = -uik_ujk_Sij
                / (tr_grad_u
                   + Kokkos::Experimental::epsilon_v<
                     decltype(tr_grad_u)>); // prevent divide by zero
    nu.data()[offset] = Kokkos::max(nu_sgs, static_cast<decltype(nu_sgs)>(0));
  }
};

namespace detail {
/*
 * The anisotropic minimum dissipation model by Abkar, Bae and Moin (2016).
 * Boussinesq effect not implemented yet.
 */
template<bool has_buoyancy>
void calculate_eddy_viscosity_impl(
  HaloView<Real***> const&             nu,
  AnisotropicMinimumDissipation const& sgs_model,
  Tensor33Field<Real***> const&        grad_u,
  Mesh const&                          mesh,
  double                               time,
  Vector3Field<Real***> const*         grad_b)
{
  using alps::square;
  using Kokkos::parallel_for;

  constexpr auto PI = Kokkos::numbers::pi_v<double>;

  auto const region = Kokkos::Profiling::ScopedRegion("SGS viscosity");

  const Grid& grid      = mesh.grid;
  const auto  nx_global = mesh.global_extent(0);
  const auto  ny_global = mesh.global_extent(1);
  const auto  nx        = nu.extent_int(0);
  const auto  ny        = nu.extent_int(1);
  const auto  nz        = mesh.extent(2);

  const auto dx = 2 * PI / (double)mesh.pex / nx_global;
  const auto dy = 2 * PI / (double)mesh.pey / ny_global;

  // The filter width considers the effective resolution from dealias
  const auto dxf = Real(dx * sgs_model.options_.horizontal_filter_width_scale);
  const auto dyf = Real(dy * sgs_model.options_.horizontal_filter_width_scale);
  const auto C0  = Real(sgs_model.options_.C0);
  const auto Cz  = Real(sgs_model.options_.Cz);

  const auto stream1 = get_next_stream();
  const auto stream2 = get_next_stream();

  auto reqs = async_update_halo_lower_z(grid.partition(), grad_u.xz, 1);
  reqs.push(async_update_halo_lower_z(grid.partition(), grad_u.yz, 2));
  reqs.push(async_update_halo_lower_z(grid.partition(), grad_u.zx, 3));
  reqs.push(async_update_halo_lower_z(grid.partition(), grad_u.zy, 4));
  if constexpr (has_buoyancy) {
    reqs.push(async_update_halo_lower_z(grid.partition(), grad_b->x, 5));
    reqs.push(async_update_halo_lower_z(grid.partition(), grad_b->y, 6));
  }

  auto constexpr tile_size = []() -> Kokkos::Array<range_index_t, 3> {
    using Device = Kokkos::DefaultExecutionSpace;
    if constexpr (is_cuda_execution_space_v<Device>) return {32, 1, 4};
    if constexpr (is_hip_execution_space_v<Device>) return {64, 1, 2};
    return {0, 0, 0};
  }();
  const auto is_top    = grid.comm().is_last(2);
  const auto is_bottom = grid.comm().is_first(2);

  auto const functor =
    amd_functor<has_buoyancy>{nu, grad_u, grad_b, mesh, dxf, dyf, C0, Cz};
  parallel_for("nu_sgs amd",
               LoopPolicy<3>(
                 stream1, {0, 0, 1}, {nx, ny, is_top ? nz - 1 : nz}, tile_size),
               functor);

  reqs.waitall();
  if (is_top) {
    Kokkos::deep_copy(stream2, subview(nu, ALL, ALL, nz - 1).view(), 0);
  }
  if (is_bottom) {
    Kokkos::deep_copy(stream2, subview(nu, ALL, ALL, 0).view(), 0);
  } else {
    parallel_for("nu_sgs amd 0",
                 LoopPolicy<3>(stream2, {0, 0, 0}, {nx, ny, 1}, tile_size),
                 functor);
  }
  stream1.fence();
  stream2.fence();

  if (sgs_model.options_.enable_debug_C0) {
    std::string const      logger_name = "eddy_viscosity_file_logger";
    RotatingFileSinkConfig file_config{"logs/eddy_viscosity.log"};
    detail::report_nu_t_stats(
      create_inner_view(nu).view(), time, mesh, logger_name, file_config);
  }

  if (sgs_model.keep_grad) {
    sgs_model.grad_u = grad_u;
  }
}
} // namespace detail

void calculate_eddy_viscosity(
  HaloView<Real***> const&                     nu,
  AnisotropicMinimumDissipation const&         sgs_model,
  Tensor33Field<Real***> const&                grad_u,
  FlowField const&                             flow,
  BoussinesqForce<ChannelFlowSolverAB2> const* buoyancy)
{
  if (buoyancy != nullptr) {
    // calculate buoyancy force and its gradient
    auto const stream        = get_next_stream();
    auto const [nx, ny, nz]  = flow.mesh.extents();
    auto const buoyancy_flux = HaloView<Real***, default_memory_pool>(
      "buoyancy_flux", {0, nx - 1}, {0, ny - 1}, {-1, nz});
    // the above should initialize the array to zero, the following is just to
    // make sure it is zero-initialized
    Kokkos::deep_copy(stream, buoyancy_flux.view(), 0);
    // calculate the buoyance term then take gradient
    add_buoyancy_force_without_ref_scalar(buoyancy_flux,
                                          buoyancy->get_scalar()->array,
                                          buoyancy->get_Ri(),
                                          flow.mesh,
                                          stream);
    stream.fence();
    update_halo_lower_z(flow.mesh.grid.partition(), buoyancy_flux, 1);
    auto const grad_b = grad(buoyancy_flux, flow.mesh, NodePt());
    detail::calculate_eddy_viscosity_impl<true>(
      nu, sgs_model, grad_u, flow.mesh, flow.time, &grad_b);
  } else {
    detail::calculate_eddy_viscosity_impl<false>(
      nu, sgs_model, grad_u, flow.mesh, flow.time, nullptr);
  }
}

void AnisotropicMinimumDissipation::cleanup_grad() const
{
  grad_u = {};
}

std::string AnisotropicMinimumDissipation::show() const
{
  return fmt::format(
    "{} (Poincare constant = {}, {})", name, options_.C0, options_.Cz);
}
} // namespace alps::solver
