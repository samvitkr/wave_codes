#include "solver.h"

#include "rhs.h"
#include "series.h"
#include "smooth.h"
#include "stats.h"

#include <common/base/logging.h>
#include <common/device/devices.h>
#include <common/kokkos_abstraction/pool_space.h>
#include <common/runtime/async_utils.h>
#include <io/hdf5.h>
#include <spectral/spectral.h>

#include <Kokkos_Profiling_ScopedRegion.hpp>

namespace alps::solver::hos {

HOSSolver::HOSSolver(const HOSField& wave_field,
                     Real            Fr2_,
                     Real            We_,
                     int             expansion_order)
  : solution{wave_field}
  , grid{wave_field.grid}
  , Fr2{Fr2_}
  , We{We_}
  , order{expansion_order}
  , integrator_{ORK256()}
  , smoother_{default_smoother()}
  , logger{get_logger("hos")}
{
  if (Fr2 <= 0) {
    throw std::runtime_error("Fr2 must be positive");
  }
  if (We < 0) {
    throw std::runtime_error("We must be non-negative");
  }
  if (We <= std::numeric_limits<Real>::epsilon()) {
    // if We is too small, the surface tension term is disabled
    if (grid.comm().rank() == 0) {
      logger->warn("Surface tension term is disabled (We = {})", We);
    }
    We = std::numeric_limits<Real>::infinity();
  }
  if (grid.extent(2) < expansion_order) {
    throw std::runtime_error("Expansion order is too large for grid");
  }

  auto r2c_layout =
    grid.get_r2c_xy_output_layout<Real, Kokkos::DefaultExecutionSpace>();
  wvn_hos = MDView<Real***, default_memory_pool>(
    Kokkos::view_alloc("wavenumber", Kokkos::WithoutInitializing),
    r2c_layout.dimension[0],
    r2c_layout.dimension[1],
    expansion_order);
  calc_wavenumbers(wvn_hos, (Real)grid.pex, (Real)grid.pey, grid);
}

std::unique_ptr<Smoother> HOSSolver::default_smoother()
{
  return std::make_unique<LowPassFilter>(0.5);
}

void HOSSolver::set_smoother(std::unique_ptr<Smoother> smoother)
{
  smoother_ = std::move(smoother);
}

HOSSolver::HOSSolver(const HOSField& wave_field, ConfigTable const& config)
  : HOSSolver(wave_field,
              config.get_value<Real>("Fr2"),
              config.get_value<Real>("We"),
              config.get_value<int>("order"))
{
  if (config.contains("integrator")) {
    auto const name = config.get_value<std::string>("integrator");
    integrator_     = make_integrator(name);
  }
  if (config.contains("smoother")) {
    set_smoother(make_smoother(config.extract_table("smoother")));
  }
}

void HOSSolver::do_step(Real dt, PaCallBackFcn const& pa_callback) const
{
  auto const region = Kokkos::Profiling::ScopedRegion("HOS step");
  std::visit([this, dt, pa_callback](
               auto const& method) { do_step(method, dt, pa_callback); },
             integrator_);
}

void HOSSolver::do_step(RK4 const& /*integrator*/,
                        Real const           dt,
                        PaCallBackFcn const& pa_callback) const
{
  auto const&                             sol = solution.value;
  auto const&                             pa  = solution.pa;
  MDView<Real** [2], default_memory_pool> tmp_sol(
    Kokkos::view_alloc("tmp sol", Kokkos::WithoutInitializing),
    sol.extent(0),
    sol.extent(1));
  MDView<Real** [2][4], default_memory_pool> Ft(
    Kokkos::view_alloc("Ft", Kokkos::WithoutInitializing),
    sol.extent(0),
    sol.extent(1));

  auto stream = get_next_stream();

  // k1 = rhs(t, F(t))
  if (pa_callback) {
    pa_callback(pa, get_time(), stream);
  }
  update_rhs(subview(Ft, ALL, ALL, ALL, 0), sol, pa, stream);

  std::array<std::pair<int, Real>, 3> rk_coeffs{
    std::pair{1, dt / 2}, std::pair{2, dt / 2}, std::pair{3, dt}};
  for (auto const& [ii, aih] : rk_coeffs) {
    // F_i = F(t) + ai * k_i * h
    update_intermediate_solution_rk4_async(
      tmp_sol, sol, subview(Ft, ALL, ALL, ALL, ii - 1), aih, stream);
    if (smoother_) {
      smoother_->apply(tmp_sol, grid, stream);
    }

    // k_{i+1} = rhs(t + DT, F_i)
    if (pa_callback) {
      pa_callback(pa, get_time() + (double)aih, stream);
    }
    update_rhs(subview(Ft, ALL, ALL, ALL, ii), tmp_sol, pa, stream);
  }

  // F(t+dt) = F(t) + (k1 + 2k2 + 2k3 + k4)/6*dt
  update_final_solution_rk4_async(sol, Ft, dt, stream);
  if (smoother_) {
    smoother_->apply(sol, grid, stream);
  }

  stream.fence();

  ALPS_CHECK_LAST_DEVICE_ERROR();

  if (!check_validity(solution.eta(), grid)) {
    throw std::runtime_error("NaN or Inf detected in solution");
  }
}

void HOSSolver::do_step(LowStorageRK2N const& integrator,
                        Real const            dt,
                        PaCallBackFcn const&  pa_callback) const
{
  auto const&                             sol = solution.value;
  auto const&                             pa  = solution.pa;
  MDView<Real** [2], default_memory_pool> Ft(
    Kokkos::view_alloc("Ft", Kokkos::WithoutInitializing),
    sol.extent(0),
    sol.extent(1));
  MDView<Real** [2], default_memory_pool> sol_W(
    Kokkos::view_alloc("tmp W", Kokkos::WithoutInitializing),
    sol.extent(0),
    sol.extent(1));

  auto stream = get_next_stream();

  for (int s = 0; s < integrator.S; ++s) {
    if (pa_callback) {
      pa_callback(pa, get_time() + integrator.c.at(s) * (double)dt, stream);
    }
    update_rhs(Ft, sol, pa, stream);
    update_solution_williamson_async(sol, sol_W, Ft, integrator, s, dt, stream);
    if (smoother_) {
      smoother_->apply(sol, grid, stream);
    }
  }

  stream.fence();

  ALPS_CHECK_LAST_DEVICE_ERROR();

  if (!check_validity(solution.eta(), grid)) {
    throw std::runtime_error("NaN or Inf detected in solution");
  }
}

void HOSSolver::do_step(LowStorageRK2C const& integrator,
                        Real const            dt,
                        PaCallBackFcn const&  pa_callback) const
{
  auto const&                             sol = solution.value;
  auto const&                             pa  = solution.pa;
  MDView<Real** [2], default_memory_pool> Ft(
    Kokkos::view_alloc("Ft", Kokkos::WithoutInitializing),
    sol.extent(0),
    sol.extent(1));
  MDView<Real** [2], default_memory_pool> sol_W(
    Kokkos::view_alloc("tmp W", Kokkos::WithoutInitializing),
    sol.extent(0),
    sol.extent(1));

  auto stream = get_next_stream();

  if (pa_callback) {
    pa_callback(pa, get_time() + integrator.c[0] * (double)dt, stream);
  }
  update_rhs(Ft, sol, pa, stream);
  for (int s = 0; s < integrator.S - 1; ++s) {
    update_solution_RK2C_async(sol, sol_W, Ft, integrator, s, dt, stream);
    if (smoother_) {
      smoother_->apply(sol, grid, stream);
    }
    if (pa_callback) {
      pa_callback(pa, get_time() + integrator.c.at(s + 1) * (double)dt, stream);
    }
    update_rhs(Ft, sol, pa, stream);
  }
  update_solution_RK2C_async(
    sol, sol_W, Ft, integrator, integrator.S - 1, dt, stream);
  if (smoother_) {
    smoother_->apply(sol, grid, stream);
  }

  stream.fence();

  ALPS_CHECK_LAST_DEVICE_ERROR();

  if (!check_validity(solution.eta(), grid)) {
    throw std::runtime_error("NaN or Inf detected in solution");
  }
}

void HOSSolver::update_rhs(MDView<Real** [2]> const&            dFdt,
                           MDView<Real const** [2]> const&      sol,
                           MDView<Real const**> const&          pa,
                           Kokkos::DefaultExecutionSpace const& space) const
{
  auto const region = Kokkos::Profiling::ScopedRegion("HOS RHS");

  MDView<Real**, default_memory_pool> ws(
    Kokkos::view_alloc("w_s", Kokkos::WithoutInitializing),
    sol.extent(0),
    sol.extent(1));

  get_ws(ws, sol, space);

  calc_evolution_rhs(dFdt, sol, ws, pa, Fr2, We, grid, space);

  ALPS_CHECK_LAST_DEVICE_ERROR();
}

void HOSSolver::get_ws(MDView<Real**> const&                ws,
                       MDView<const Real** [2]> const&      sol,
                       Kokkos::DefaultExecutionSpace const& space) const
{
  MDView<Real***, default_memory_pool> zp_hos(
    Kokkos::view_alloc("zp_hos", Kokkos::WithoutInitializing),
    sol.extent(0),
    sol.extent(1),
    order - 1);
  MDView<Real***, default_memory_pool> r_hat(
    Kokkos::view_alloc("r_hat", Kokkos::WithoutInitializing), wvn_hos.layout());

  taylor_series_coeff_async(
    zp_hos, subview(sol, ALL, ALL, 0), order, grid, space);

  surface_vp_expansion(
    r_hat, subview(sol, ALL, ALL, 1), zp_hos, wvn_hos, grid, space);

  surface_w(ws, r_hat, zp_hos, wvn_hos, grid, space);

  ALPS_CHECK_LAST_DEVICE_ERROR();
}

Real HOSSolver::estimate_max_omega() const
{
  auto kx_max = grid.pex * int(grid.global_extent(0) / 3); // dealias considered
  auto ky_max = grid.pey * int(grid.global_extent(1) / 3); // dealias considered
  auto k_max  = (Real)std::hypot(kx_max, ky_max);
  auto omega2 = k_max * (1 / Fr2 + k_max * k_max / We);

  return std::sqrt(omega2);
}

void HOSSolver::save(std::filesystem::path filename) const
{
  if (grid.comm().rank() == 0) {
    logger->info("Saving solution to file {}", filename.string());
  }
  solution.save(filename);
}

void HOSSolver::load(std::filesystem::path filename)
{
  Kokkos::fence();

  if (grid.comm().rank() == 0) {
    logger->info("Loading solution from file {}", filename.string());
  }

  if (!std::filesystem::exists(filename)) {
    throw std::runtime_error("Data file " + filename.string()
                             + " does not exist");
  }

  auto h5file = io::hdf5::open_file_with_mpi(
    filename.string(), HighFive::File::ReadOnly, grid.comm().raw_handle());

  if (h5file.exist("time")) {
    auto t_dset = h5file.getDataSet("time");
    t_dset.read(solution.time);
  }

  const BlockPartition& partition = grid.partition();

  const std::vector total_shape_xy{partition.global_extents[0],
                                   partition.global_extents[1]};
  const std::vector block_shape_xy{partition.extents[0], partition.extents[1]};
  const std::vector offset_xy{partition.offsets[0], partition.offsets[1]};
  io::hdf5::read_blocks(
    h5file, "eta", solution.eta(), total_shape_xy, block_shape_xy, offset_xy);

  io::hdf5::read_blocks(
    h5file, "vps", solution.vps(), total_shape_xy, block_shape_xy, offset_xy);

  Kokkos::fence();
}

std::ostream& operator<<(std::ostream& os, HOSSolver const& solver)
{
  auto const& grid = solver.grid;
  os << fmt::format("HOS solver (order = {})\n", solver.order);
  os << fmt::format(
    "Nx x Ny = {} x {}\n", grid.global_extent(0), grid.global_extent(1));
  os << fmt::format("Lx x Ly = {} x {}\n",
                    2 * Kokkos::numbers::pi_v<double> / grid.pex,
                    2 * Kokkos::numbers::pi_v<double> / grid.pey);
  os << fmt::format("Fr2 = {}, We = {}\n", solver.Fr2, solver.We);
  os << fmt::format("Smoother: {}\n",
                    solver.smoother_ ? solver.smoother_->info() : "none");
  os << fmt::format(
    "Integrator: {}",
    std::visit([](auto const& m) { return m.name; }, solver.integrator_));

  return os;
}

} // namespace alps::solver::hos
