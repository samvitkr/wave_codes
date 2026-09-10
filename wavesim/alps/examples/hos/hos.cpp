#include <apps/utils/git_version.h>
#include <apps/utils/progress_tracker.h>
#include <apps/utils/std_cli_opts.h>
#include <apps/utils/symlinks.h>
#include <common/base/logging.h>
#include <common/device/devices.h>
#include <common/math.h>
#include <common/program_options/config_table.h>
#include <common/runtime/manager.h>
#include <mpipp/collectives.h>
#include <solvers/hos/solver.h>
#include <solvers/hos/stats.h>
#include <spectral/spectral.h>

#include <filesystem>
#include <sstream>

namespace fs = std::filesystem;

constexpr auto app_name = "alps_hos";

namespace alps {
struct SaveOptions
{
  std::string restart_prefix;
  std::string solution_prefix;
  int         save_restart_frequency;
  int         save_solution_frequency;
  int         in_situ_stats_frequency;
};

static std::string time_to_number(double t)
{
  return fmt::format("{:014d}", std::llround(t * 1e8));
}

static void save_restart(solver::hos::HOSSolver const& solver,
                         SaveOptions const&            save_options,
                         ConfigTable const&            config)
{
  const auto save_filename =
    save_options.restart_prefix + time_to_number(solver.get_time()) + ".h5";
  solver.save(save_filename);
  if (solver.grid.comm().rank() == 0) {
    HighFive::File file(save_filename, HighFive::File::ReadWrite);
    io::hdf5::write_string(file, "config", config.to_string());
    io::hdf5::write_string(
      file, "app_version", alps::apps::git_version(app_name));

    apps::force_create_symlink(save_filename,
                               save_options.restart_prefix + ".h5");
  }
}

void integrate_solver(solver::hos::HOSSolver& solver,
                      int                     n_steps,
                      double                  max_time,
                      double                  dt,
                      SaveOptions const&      save_options,
                      ConfigTable const&      config)
{
  auto world = solver.grid.comm().comm;

  Kokkos::fence();
  mpipp::barrier(world);
  apps::ProgressTracker progress_tracker(
    n_steps, solver.get_time(), max_time); // Set up progress tracker

  start_profiling();

  double error_time = 0;
  while (!progress_tracker.is_complete()) {
    const auto step = progress_tracker.current();
    Kokkos::Profiling::pushRegion(fmt::format("Step {}", step));

    const auto new_time =
      alps::accumulate_compensated(solver.get_time(), dt, error_time);

    if (world.rank() == 0) {
      spdlog::info(
        "({}/{}) t = {:.8f}", step, progress_tracker.total(), new_time);
    }

    solver.do_step((Real)dt);
    solver.set_time(new_time);

    if (step % save_options.save_solution_frequency == 0) {
      auto filename = save_options.solution_prefix
                    + time_to_number(solver.get_time()) + ".h5";
      solver.solution.save(filename);
      if (world.rank() == 0) {
        HighFive::File file(filename, HighFive::File::ReadWrite);
        io::hdf5::write_string(file, "config", config.to_string());
        io::hdf5::write_string(
          file, "app_version", alps::apps::git_version(app_name));
      }
    }

    if (step % save_options.save_restart_frequency == 0) {
      save_restart(solver, save_options, config);
    }
    if (step % save_options.in_situ_stats_frequency == 0) {
      solver::hos::report_wave_stats(solver.solution, new_time);
    }

    progress_tracker.update(new_time);
    if (world.rank() == 0) {
      if (auto output = progress_tracker.output(); output) {
        spdlog::info("{}", output.value());
      }
    }

    Kokkos::Profiling::popRegion();
  }

  stop_profiling();

  if (world.rank() == 0) {
    spdlog::info("{}", progress_tracker.complete_output());
  }
}

void program(ConfigTable const&      root_config,
             alps::CLIOptions const& cli_options)
{
  alps::RuntimeManager::instance().init_runtimes();

  // Load values from TOML config file
  constexpr auto config_section = std::string_view("hos_solver");
  const auto     config         = root_config.extract_table(config_section);

  const auto n_proc      = config.get_value<int>("parallel_decompose");
  const auto grid_size   = config.get_value<std::vector<int>>("grid_size");
  const auto domain_size = config.get_value<std::vector<double>>("domain_size");
  const auto order       = config.get_value<int>("order");
  if (grid_size.size() < 2 || domain_size.size() < 2) {
    throw std::invalid_argument("Grid size and domain size must be 2D.");
  }

  const auto n_steps = config.get_value<int>("NStep");
  const auto max_time =
    config.get_value_or<double>("MaxTime", std::numeric_limits<double>::max());
  const auto        dt = config.get_value<double>("dt");
  SaveOptions const save_options{
    "restart_wave",
    "wave",
    config.get_value_or("NRestart", std::numeric_limits<int>::max()),
    config.get_value_or("NOutD", std::numeric_limits<int>::max()),
    config.get_value_or("NOutC", std::numeric_limits<int>::max())};

  const std::vector periodic{1, 1, 0};
  const std::vector proc_dims3{1, n_proc, 1};
  const MPIComm3D   comm(mpipp::COMM_WORLD(), proc_dims3, periodic);

  const Grid grid(PencilPlan(comm, {grid_size[0], grid_size[1], order}),
                  domain_size[0],
                  domain_size[1]);

  solver::hos::HOSField hosField(grid);

  solver::hos::HOSSolver solver(hosField, config);

  if (comm.rank() == 0) {
    std::stringstream ss;
    ss << solver;
    spdlog::info(ss.str());
    spdlog::info("max omega*dt ~ {:.2g}",
                 (double)solver.estimate_max_omega() * dt);
  }

  auto const restart_id =
    cli_options.has("restart-from")
      ? fmt::format("{:014d}", cli_options.get<std::int64_t>("restart-from"))
      : std::string{};
  fs::path const restart_filepath =
    save_options.restart_prefix + restart_id + ".h5";

  solver.load(restart_filepath);

  integrate_solver(solver, n_steps, max_time, dt, save_options, root_config);

  if (n_steps % save_options.save_restart_frequency != 0) {
    save_restart(solver, save_options, root_config);
  }
}
} // namespace alps

int alps_main(int argc, char* argv[])
{
  alps::CLIOptions app_options(app_name);

  try {
    alps::apps::add_std_cli_options(app_options);
    app_options.parse(argc, argv);

    if (app_options.has("h")) {
      // Print help message and exit
      std::cout << app_options.help() << std::endl;
      return 0;
    }

    // Load TOML configuration file
    const auto config = alps::ConfigTable::parse_from_file(
      app_options.get<std::string>("config"));

    alps::program(config, app_options);
  } catch (std::exception const& e) {
    spdlog::error("{}", e.what());
    return 1;
  } catch (...) {
    spdlog::error("Caught unknown errors.");
    return 1;
  }
  return 0;
}
