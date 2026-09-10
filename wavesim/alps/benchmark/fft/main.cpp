#include <common/base/logging.h>
#include <common/container/view_types.h>
#include <common/program_options/cli_options.h>
#include <common/runtime/manager.h>
#include <common/utils/bench.h>
#include <decomp/mdcomm.h>
#include <spectral/spectral_base.h>

#include <Kokkos_Core.hpp>
#include <mpipp/environment.h>

#include <algorithm>
#include <cctype>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>

template<class T, class ES, class Backend>
class BenchmarkCase
{
 public:
  int nx, ny, nz;

  alps::MPIComm3D md_comm;

  alps::PencilPlan                                         grid;
  std::unique_ptr<alps::spectral::SpectralPlanBase<T, ES>> plan;

  alps::MDView<T***, ES> data1;
  alps::MDView<T***, ES> data2;

  BenchmarkCase(mpipp::communicator const& world,
                int                        npx,
                int                        npz,
                int                        nx_,
                int                        ny_,
                int                        nz_)
    : nx{nx_}
    , ny{ny_}
    , nz{nz_}
    , md_comm(world, {1, npx, npz}, {1, 1, 0})
    , grid(md_comm, {nx, ny, nz})
    , plan(alps::spectral::SpectralPlanFactory::create<T, ES, Backend>(grid,
                                                                       0.5,
                                                                       3.0))
    , data1("data1", alps::create_local_layout(grid))
    , data2("data2", alps::create_local_layout(grid))
  {
    auto data1_host = Kokkos::create_mirror_view(data1);
    for (int i = 0; i < data1.extent_int(0); ++i) {
      auto x              = 2 * Kokkos::numbers::pi_v<T> / nx * i;
      data1_host(i, 0, 0) = Kokkos::cos(4 * x);
    }
    Kokkos::deep_copy(data1, data1_host);
    Kokkos::deep_copy(data2, data1);
  }

  void bench_ddx(alps::bench::Bench* bench, std::string base_name)
  {
    bench->name(base_name + " ddx");

    alps::bench::detail::IterationLogic iterationLogic(*bench);
    auto& pc = alps::bench::detail::performanceCounters();

    while (auto n = iterationLogic.numIters()) {
      pc.beginMeasure();
      alps::bench::Clock::time_point const before = alps::bench::Clock::now();
      while (n-- > 0) {
        plan->do_ddx(data1, data2, ES());
      }
      ES().fence();
      alps::bench::Clock::time_point const after = alps::bench::Clock::now();
      pc.endMeasure();
      pc.updateResults(iterationLogic.numIters());
      iterationLogic.add(after - before, pc);
    }
    iterationLogic.moveResultTo(bench->mResults);
  }

  void bench_ddy(alps::bench::Bench* bench, std::string base_name)
  {
    bench->name(base_name + " ddy");

    alps::bench::detail::IterationLogic iterationLogic(*bench);
    auto& pc = alps::bench::detail::performanceCounters();

    while (auto n = iterationLogic.numIters()) {
      pc.beginMeasure();
      alps::bench::Clock::time_point const before = alps::bench::Clock::now();
      while (n-- > 0) {
        plan->do_ddy(data1,
                     data2,
                     alps::spectral::SpectralPostOp::AssignAfterTranspose,
                     ES());
      }
      ES().fence();
      alps::bench::Clock::time_point const after = alps::bench::Clock::now();
      pc.endMeasure();
      pc.updateResults(iterationLogic.numIters());
      iterationLogic.add(after - before, pc);
    }
    iterationLogic.moveResultTo(bench->mResults);
  }
};

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

template<typename T>
void run_benchmarks(alps::bench::Bench&          b,
                    mpipp::communicator const&   world,
                    int                          np1,
                    int                          np2,
                    int                          nx,
                    int                          ny,
                    int                          nz,
                    std::set<std::string> const& backends,
                    std::set<std::string> const& ops)
{
  auto const type_name = []() -> std::string {
    if constexpr (std::is_same_v<T, float>) {
      return "float";
    } else {
      return "double";
    }
  }();

#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
  if (backends.count("vkfft")) {
    BenchmarkCase<T, Kokkos::DefaultExecutionSpace, alps::fft::VKFFT> benchCase(
      world, np1, np2, nx, ny, nz);
    auto const name = type_name + " "
                    + std::string(Kokkos::DefaultExecutionSpace::name())
                    + " VkFFT";
    if (ops.count("ddx") != 0u) {
      benchCase.bench_ddx(&b, name);
    }
    if (ops.count("ddy") != 0u) {
      benchCase.bench_ddy(&b, name);
    }
  }
#endif
#if defined(KOKKOS_ENABLE_CUDA)
  if (backends.count("cufft")) {
    BenchmarkCase<T, Kokkos::DefaultExecutionSpace, alps::fft::CUFFT> benchCase(
      world, np1, np2, nx, ny, nz);
    auto const name = type_name + " "
                    + std::string(Kokkos::DefaultExecutionSpace::name())
                    + " cuFFT";
    if (ops.count("ddx") != 0u) {
      benchCase.bench_ddx(&b, name);
    }
    if (ops.count("ddy") != 0u) {
      benchCase.bench_ddy(&b, name);
    }
  }
#endif
  if (backends.count("fftw")) {
    BenchmarkCase<T, Kokkos::DefaultHostExecutionSpace, alps::fft::FFTW>
               benchCase(world, np1, np2, nx, ny, nz);
    auto const name = type_name + " "
                    + std::string(Kokkos::DefaultHostExecutionSpace::name())
                    + " FFTW";
    if (ops.count("ddx") != 0u) {
      benchCase.bench_ddx(&b, name);
    }
    if (ops.count("ddy") != 0u) {
      benchCase.bench_ddy(&b, name);
    }
  }
}

int alps_main(int argc, char* argv[])
{
  alps::CLIOptions app_options("Benchmark");
  app_options.add_option_with_default<int>(
    "np1", "Number of processes in x-direction", "1");
  app_options.add_option_with_default<int>(
    "np2", "Number of processes in y-direction", "1");
  app_options.add_option<int>("nx", "Number of grid points in x-direction");
  app_options.add_option<int>("ny", "Number of grid points in y-direction");
  app_options.add_option<int>("nz", "Number of grid points in z-direction");
  app_options.add_option_with_default<std::string>(
    "real-type",
    "Comma-separated list of real types to benchmark: float,double "
    "(aliases: single,fp32,f32,fp64,f64)",
    "double");
  app_options.add_option_with_default<std::string>(
    "backend",
    "Comma-separated list of FFT backends to benchmark: vkfft,cufft,fftw. "
    "Backends not available in the current build are silently skipped.",
    "vkfft,cufft,fftw");
  app_options.add_option_with_default<std::string>(
    "ops",
    "Comma-separated list of spectral operations to benchmark: ddx,ddy",
    "ddx,ddy");
  app_options.add_option_with_default<int>(
    "min-epoch-iterations",
    "Minimum number of iterations per benchmark epoch (default: 3)",
    "3");
  app_options.parse(argc, argv);

  const int np1{app_options.get<int>("np1")};
  const int np2{app_options.get<int>("np2")};
  const int nx{app_options.get<int>("nx")};
  const int ny{app_options.get<int>("ny")};
  const int nz{app_options.get<int>("nz")};

  static std::set<std::string> const valid_real_types{
    "float", "double", "single", "fp32", "f32", "fp64", "f64"};
  static std::set<std::string> const valid_backends{"vkfft", "cufft", "fftw"};
  static std::set<std::string> const valid_ops{"ddx", "ddy"};

  const std::set<std::string> real_types =
    parse_csv_option(app_options.get<std::string>("real-type"),
                     valid_real_types,
                     "--real-type",
                     "float, double (aliases: single, fp32, f32, fp64, f64)");
  const std::set<std::string> backends =
    parse_csv_option(app_options.get<std::string>("backend"),
                     valid_backends,
                     "--backend",
                     "vkfft, cufft, fftw");
  const std::set<std::string> ops = parse_csv_option(
    app_options.get<std::string>("ops"), valid_ops, "--ops", "ddx, ddy");

  alps::RuntimeManager::instance().init_runtimes();

  using std::chrono_literals::operator""us;
  alps::bench::Bench b;
  b.title("Spectral derivative")
    .timeUnit(1us, "us")
    .warmup(2)
    .minEpochIterations(app_options.get<int>("min-epoch-iterations"));

  bool const run_float =
    (real_types.count("float") != 0u) || (real_types.count("single") != 0u)
    || (real_types.count("fp32") != 0u) || (real_types.count("f32") != 0u);
  bool const run_double = (real_types.count("double") != 0u)
                       || (real_types.count("fp64") != 0u)
                       || (real_types.count("f64") != 0u);

  if (run_float) {
    run_benchmarks<float>(
      b, mpipp::COMM_WORLD(), np1, np2, nx, ny, nz, backends, ops);
  }

  if (run_double) {
    run_benchmarks<double>(
      b, mpipp::COMM_WORLD(), np1, np2, nx, ny, nz, backends, ops);
  }

  return 0;
}
