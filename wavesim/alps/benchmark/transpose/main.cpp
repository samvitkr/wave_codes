#include <common/container/view_types.h>
#include <common/program_options/cli_options.h>
#include <common/runtime/manager.h>
#include <common/utils/bench.h>
#include <transpose/transposer_base.h>

#include <Kokkos_Core.hpp>

#include <algorithm>
#include <cctype>
#include <chrono>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

template<class T, class ES>
class transpose_case
{
 public:
  using transposer_t = std::unique_ptr<alps::transpose::TransposerBase<T, ES>>;

  alps::MDView<T***, ES> data_x;
  alps::MDView<T***, ES> data_y;

  transposer_t transposer_x2y;
  transposer_t transposer_y2x;

  transpose_case(mpipp::communicator const&         comm,
                 int                                nx,
                 int                                ny,
                 int                                nz,
                 alps::transpose::TransposerOptions options)
    : data_x(Kokkos::view_alloc("x_layout", Kokkos::WithoutInitializing),
             nx,
             ny / comm.size(),
             nz)
    , data_y(Kokkos::view_alloc("y_layout", Kokkos::WithoutInitializing),
             ny,
             nx / comm.size(),
             nz)
    , transposer_x2y(
        alps::transpose::create_transposer<T, ES>(comm, nx, ny, nz, options))
    , transposer_y2x(
        alps::transpose::create_transposer<T, ES>(comm, ny, nx, nz, options))
  {
    Kokkos::deep_copy(data_x, static_cast<T>(0));
    Kokkos::deep_copy(data_y, static_cast<T>(0));
  }

  void execute(ES const& space) const
  {
    transposer_x2y->execute(
      data_y, data_x, alps::transpose::TransposeOpAssign(), space);
    transposer_y2x->execute(
      data_x, data_y, alps::transpose::TransposeOpAssign(), space);
    space.fence();
  }
};

template<typename T, typename ExecSpace>
void run_benchmarks(int nx, int ny, int nz, ExecSpace const& space)
{
  static_assert(
    std::is_same_v<ExecSpace, Kokkos::DefaultExecutionSpace>
      || std::is_same_v<ExecSpace, Kokkos::DefaultHostExecutionSpace>,
    "run_benchmarks instantiated with an unexpected ExecSpace");

  auto comm       = mpipp::COMM_WORLD();
  auto space_name = std::string(ExecSpace::name());
  auto type_name  = std::string("double");
  if constexpr (std::is_same_v<T, float>) {
    type_name = "float";
  }

  using alps::transpose::TransposeMethod;
  using alps::transpose::TransposerOptions;
  std::vector<std::pair<std::string, TransposerOptions>> test_options;
  if (comm.size() == 1) {
    TransposerOptions options;
    options.method          = TransposeMethod::Single;
    options.tune_tile_sizes = true;
    test_options.emplace_back(type_name + " " + space_name + " single",
                              options);
  } else {
    TransposerOptions options;
    options.tune_tile_sizes = false;

    options.method = TransposeMethod::All2All;
    test_options.emplace_back(type_name + " " + space_name + " all2all",
                              options);

    options.method = TransposeMethod::Point2Point;
    test_options.emplace_back(type_name + " " + space_name + " p2p", options);

    if constexpr (std::is_same_v<ExecSpace,
                                 Kokkos::DefaultHostExecutionSpace>) {
      options.method = TransposeMethod::Point2PointSHM;
      test_options.emplace_back(type_name + " " + space_name + " p2pSHM",
                                options);
    }
  }
  for (const auto& options : test_options) {
    transpose_case<T, ExecSpace> benchCase(comm, nx, ny, nz, options.second);

    comm.barrier();
    using std::chrono_literals::operator""ms;
    alps::bench::Bench()
      .timeUnit(1ms, "ms")
      .warmup(1)
      .minEpochIterations(comm.size() > 1 ? 5 : 3)
      .run(options.first, [&]() { benchCase.execute(space); });
  }
  space.fence();
}

template<typename T>
void run_benchmarks_all_exec_spaces(int nx, int ny, int nz)
{
  run_benchmarks<T>(nx, ny, nz, Kokkos::DefaultExecutionSpace());
  if constexpr (!std::is_same_v<Kokkos::DefaultExecutionSpace,
                                Kokkos::DefaultHostExecutionSpace>) {
    run_benchmarks<T>(nx, ny, nz, Kokkos::DefaultHostExecutionSpace());
  }
}

inline std::string to_lower(std::string s)
{
  std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) {
    return static_cast<char>(std::tolower(c));
  });
  return s;
}

/// Parse a comma-separated list into a set of lowercase tokens.
/// Throws std::runtime_error for any unrecognised value.
inline std::set<std::string>
parse_csv_option(std::string const&           raw,
                 std::set<std::string> const& valid_names,
                 std::string const&           option_name,
                 std::string const&           valid_values)
{
  std::set<std::string> result;
  std::istringstream    ss(raw);
  std::string           token;
  while (std::getline(ss, token, ',')) {
    // trim whitespace
    token.erase(0, token.find_first_not_of(" \t"));
    token.erase(token.find_last_not_of(" \t") + 1);
    if (token.empty()) {
      continue;
    }
    std::string lower = to_lower(token);
    if (valid_names.find(lower) == valid_names.end()) {
      throw std::runtime_error("Unknown " + option_name + " value '" + token
                               + "'. Valid names: " + valid_values);
    }
    result.insert(lower);
  }
  return result;
}

int alps_main(int argc, char* argv[])
{
  alps::CLIOptions app_options("Benchmark");
  app_options.add_option<int>("nx", "Number of grid points in x-direction");
  app_options.add_option<int>("ny", "Number of grid points in y-direction");
  app_options.add_option<int>("nz", "Number of grid points in z-direction");
  app_options.add_option_with_default<std::string>(
    "real-type",
    "Comma-separated list of real types to benchmark: float,double "
    "(aliases: single,fp32,f32,fp64,f64)",
    "double");
  app_options.parse(argc, argv);

  const int nx{app_options.get<int>("nx")};
  const int ny{app_options.get<int>("ny")};
  const int nz{app_options.get<int>("nz")};

  static std::set<std::string> const valid_real_types{
    "float", "double", "single", "fp32", "f32", "fp64", "f64"};

  const std::set<std::string> real_types =
    parse_csv_option(app_options.get<std::string>("real-type"),
                     valid_real_types,
                     "--real-type",
                     "float, double (aliases: single, fp32, f32, fp64, f64)");

  alps::RuntimeManager::instance().init_runtimes();

  bool const run_float =
    (real_types.count("float") != 0u) || (real_types.count("single") != 0u)
    || (real_types.count("fp32") != 0u) || (real_types.count("f32") != 0u);
  bool const run_double = (real_types.count("double") != 0u)
                       || (real_types.count("fp64") != 0u)
                       || (real_types.count("f64") != 0u);

  if (run_float) {
    run_benchmarks_all_exec_spaces<float>(nx, ny, nz);
  }

  if (run_double) {
    run_benchmarks_all_exec_spaces<double>(nx, ny, nz);
  }

  return 0;
}
