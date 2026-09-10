#include <apps/utils/git_version.h>
#include <apps/utils/progress_tracker.h>
#include <apps/utils/std_cli_opts.h>
#include <apps/utils/symlinks.h>
#include <common/base/logging.h>
#include <common/device/devices.h>
#include <common/program_options/config_table.h>
#include <common/runtime/manager.h>
#include <common/utils/to_lower_case.h>
#include <diagnostics/diagnostics.h>
#include <diagnostics/field_statistics.h>
#include <io/hdf5.h>
#include <mpipp/collectives.h>
#include <solvers/mesh/mesh.h>
#include <solvers/ns/channel/solver.h>
#include <spectral/spectral.h>

namespace fs = std::filesystem;

constexpr auto app_name = "alps_channel";

struct SaveOptions
{
  std::string restart_prefix;
  std::string solution_prefix;
  std::string in_situ_stats_prefix;
  int         save_restart_frequency;
  int         save_solution_frequency;
  int         in_situ_stats_frequency;
};

static std::string time_to_number(double t)
{
  return fmt::format("{:014d}", std::llround(t * 1e8));
}

namespace alps {
template<typename SolverType>
void integrate_solver(SolverType&        solver,
                      SaveOptions const& save_options,
                      ConfigTable const& config)
{
  // Setup diagnostics
  alps::diagnostics::StatisticalDiagnosticsManager statistical_output(
    solver.flow_field.mesh);
  statistical_output.attach_statistics(diagnostics::create_velocity_statistics(
    solver.flow_field,
    std::unordered_set{diagnostics::FieldStatisticsType::XY}));
  statistical_output.attach_statistics(diagnostics::create_scalar_statistics(
    solver.flow_field,
    std::unordered_set{diagnostics::FieldStatisticsType::XY}));

  auto world = mpipp::COMM_WORLD();

  Kokkos::fence();
  mpipp::barrier(world);
  apps::ProgressTracker progress_tracker(
    solver.options.n_steps,
    solver.get_time(),
    solver.options.max_time); // Set up progress tracker

  start_profiling();

  double error_time = 0;
  while (!progress_tracker.is_complete()) {
    const auto step = progress_tracker.current();
    Kokkos::Profiling::pushRegion(fmt::format("Step {}", step));

    const auto new_time =
      alps::accumulate_compensated(solver.get_time(), solver.dt, error_time);

    if (world.rank() == 0) {
      spdlog::info(
        "({}/{}) t = {:.8f}", step, progress_tracker.total(), new_time);
    }

    solver.calc_uhat(solver.calc_explicit_rhs());

    solver.set_time(new_time);
    solver.project();
    solver.correct();
    solver.validate();

    if (step % save_options.save_solution_frequency == 0) {
      auto filename = save_options.solution_prefix
                    + time_to_number(solver.get_time()) + ".h5";
      solver.flow_field.save(filename);
      if (world.rank() == 0) {
        HighFive::File file(filename, HighFive::File::ReadWrite);
        io::hdf5::write_string(file, "config", config.to_string());
        io::hdf5::write_string(
          file, "app_version", alps::apps::git_version(app_name));
      }
    }
    if (step % save_options.save_restart_frequency == 0) {
      const auto save_filename =
        save_options.restart_prefix + time_to_number(solver.get_time()) + ".h5";
      const auto aux_filename = save_options.restart_prefix
                              + time_to_number(solver.get_time()) + "_aux.h5";
      solver.save("grid.h5", save_filename, aux_filename);
      if (world.rank() == 0) {
        HighFive::File file(save_filename, HighFive::File::ReadWrite);
        io::hdf5::write_string(file, "config", config.to_string());
        io::hdf5::write_string(
          file, "app_version", alps::apps::git_version(app_name));

        apps::force_create_symlink(save_filename,
                                   save_options.restart_prefix + ".h5");
        apps::force_create_symlink(aux_filename,
                                   save_options.restart_prefix + "_aux.h5");
      }
    }
    if (step % save_options.in_situ_stats_frequency == 0) {
      statistical_output.calculate();
      statistical_output.write_to_tecplot(save_options.in_situ_stats_prefix
                                            + time_to_number(solver.get_time()),
                                          "monitor");
      statistical_output.cleanup();
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

template<typename SolverType>
void run_instance(ConfigTable const& root_config, CLIOptions const& cli_options)
{
  alps::RuntimeManager::instance().init_runtimes();

  // Load values from TOML config file
  constexpr auto config_section = std::string_view("default_solver");
  const auto     config         = root_config.extract_table(config_section);

  const auto proc_dims =
    config.get_value<std::vector<int>>("parallel_decompose");
  const auto grid_size   = config.get_value<std::vector<int>>("grid_size");
  const auto domain_size = config.get_value<std::vector<double>>("domain_size");
  if (proc_dims.size() < 2 || grid_size.size() < 3 || domain_size.size() < 3) {
    throw std::invalid_argument(
      "Not enough inputs in parallel_decompose, grid_size or domain_size");
  }
  SaveOptions const save_options{
    "restart",
    "Sol",
    "Stat",
    config.get_value_or("NRestart", std::numeric_limits<int>::max()),
    config.get_value_or("NOutD", std::numeric_limits<int>::max()),
    config.get_value_or("NOutC", std::numeric_limits<int>::max())};

  // Set up parallelization
  const std::vector periodic{1, 1, 0};
  const std::vector proc_dims3{1, proc_dims[0], proc_dims[1]};
  const MPIComm3D   comm(mpipp::COMM_WORLD(), proc_dims3, periodic);

  const Grid grid(PencilPlan(comm, grid_size), domain_size[0], domain_size[1]);

  solver::Mesh mesh(grid, (Real)domain_size[2]);

  solver::FlowField flow(mesh);
  flow.parse_bcs_from(config);
  flow.initialize_scalars(config);

  const auto solver_options = solver::ChannelSolverOptions::parse_from(config);
  SolverType ns_solver(flow, solver_options);

  auto const restart_id =
    cli_options.has("restart-from")
      ? fmt::format("{:014d}", cli_options.get<std::int64_t>("restart-from"))
      : std::string{};
  fs::path const restart_filepath =
    save_options.restart_prefix + restart_id + ".h5";
  fs::path const restart_aux_filepath =
    save_options.restart_prefix + restart_id + "_aux.h5";
  fs::path const grid_filepath = cli_options.get<std::string>("grid");

  ns_solver.load(grid_filepath, restart_filepath, restart_aux_filepath);

  if (comm.rank() == 0) spdlog::info("{}", ns_solver.info());

  ns_solver.initialize();
  Kokkos::fence();

  integrate_solver(ns_solver, save_options, root_config);

  // Force the solver to save the state at the end of the simulation.
  if (solver_options.n_steps % save_options.save_restart_frequency != 0) {
    const auto save_filename = save_options.restart_prefix
                             + time_to_number(ns_solver.get_time()) + ".h5";
    const auto aux_filename = save_options.restart_prefix
                            + time_to_number(ns_solver.get_time()) + "_aux.h5";
    ns_solver.save("grid.h5", save_filename, aux_filename);

    if (mpipp::COMM_WORLD().rank() == 0) {
      HighFive::File file(save_filename, HighFive::File::ReadWrite);
      io::hdf5::write_string(file, "config", root_config.to_string());
      io::hdf5::write_string(
        file, "app_version", alps::apps::git_version(app_name));

      apps::force_create_symlink(save_filename,
                                 save_options.restart_prefix + ".h5");
      apps::force_create_symlink(aux_filename,
                                 save_options.restart_prefix + "_aux.h5");
    }
  }
}

} // namespace alps

void run(alps::ConfigTable const& root_config,
         alps::CLIOptions const&  cli_options)
{
  using namespace alps;
  auto integrator = to_lower_case(
    root_config.get_value_or<std::string>("default_solver.integrator", "ab2"));
  if (integrator == "ab2cn") {
    run_instance<solver::ChannelFlowSolverAB2CN>(root_config, cli_options);
  } else if (integrator == "ab2") {
    run_instance<solver::ChannelFlowSolverAB2>(root_config, cli_options);
  } else {
    throw std::invalid_argument("Unknown integrator type: " + integrator);
  }
}

int alps_main(int argc, char* argv[])
{
  // Command line options
  alps::CLIOptions app_options(app_name);
  try {
    alps::apps::add_std_cli_options(app_options);
    app_options.add_option_with_default<std::string>(
      "grid", "Specify a grid file to be read in", "grid.h5", "filename");

    app_options.parse(argc, argv);

    if (app_options.has("h")) {
      // Print help message and exit
      std::cout << app_options.help() << std::endl;
      return 0;
    }

    // Load TOML configuration file
    const auto config = alps::ConfigTable::parse_from_file(
      app_options.get<std::string>("config"));

    run(config, app_options);
  } catch (std::exception const& e) {
    spdlog::error("{}", e.what());
    return 1;
  } catch (...) {
    spdlog::error("Caught unknown errors.");
    return 1;
  }

  return 0;
}
