#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_tostring.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_templated.hpp>

#include <common/container/view_utils.h>
#include <common/kokkos_abstraction/exec_policy.h>
#include <common/kokkos_abstraction/pool_space.h>
#include <decomp/mdcomm.h>
#include <decomp/pencil_plan.h>
#include <mpipp/collectives.h>
#include <solvers/mesh/curvilinear_mesh.h>
#include <solvers/ns/curvilinear_common/diffusion_cn_zeta_coeff.h>
#include <solvers/ns/curvilinear_common/diffusion_cn_zeta_eqn.h>
#include <spectral/spectral.h>

#include <fmt/format.h>

namespace {
using namespace alps;
using namespace alps::solver;

struct MeshFixture
{
  MPIComm3D              comm;
  PencilPlan             decomp;
  spectral::SpectralGrid spectral_grid;
  BottomWaveMesh         mesh;

  explicit MeshFixture(mpipp::communicator comm_)
    : comm(comm_, {1, 1, comm_.size()}, {1, 1, 0})
    , decomp(comm, {6, 4, 17})
    , spectral_grid(decomp, 1.0, 1.0)
    , mesh(spectral_grid, 1.0)
  {
    std::vector<Real> z_node(17, 0);
    for (int k = 0; k < 15; ++k) {
      Real const s = Real(k) / Real(15);
      z_node[k]    = s * s;
    }
    z_node[15] = Real(1);
    mesh.load(z_node);

    auto      exr_h  = Kokkos::create_mirror_view(mesh.exr);
    auto      eyr_h  = Kokkos::create_mirror_view(mesh.eyr);
    auto      invj_h = Kokkos::create_mirror_view(mesh.invJ);
    int const nx     = mesh.extent(0);
    int const ny     = mesh.extent(1);
    for (int j = 0; j < mesh.extent(1); ++j) {
      for (int i = 0; i < mesh.extent(0); ++i) {
        Real const xi = nx > 1 ? Real(i) / Real(nx - 1) - Real(0.5) : Real(0);
        Real const yj = ny > 1 ? Real(j) / Real(ny - 1) - Real(0.5) : Real(0);
        Real const variation =
          Real(0.06) * xi + Real(0.04) * yj + Real(0.02) * xi * yj;
        exr_h(i, j)  = Real(0.0);
        eyr_h(i, j)  = Real(0.0);
        invj_h(i, j) = Real(1.0) + variation;
      }
    }
    Kokkos::deep_copy(mesh.exr, exr_h);
    Kokkos::deep_copy(mesh.eyr, eyr_h);
    Kokkos::deep_copy(mesh.invJ, invj_h);
  }
};

KOKKOS_FUNCTION Real reference_solution_value(int i, int j, int global_k)
{
  double const base = 0.7 + 0.05 * i + 0.09 * j + 1.3 * global_k
                    - 0.8 * global_k * global_k - 0.02 * i * j;
  return (Real)base;
}

void build_manufactured_system(MDView<Real***> const& ref,
                               MDView<Real***> const& rhs,
                               DiffusionCNBCType      bottom_bc,
                               DiffusionCNBCType      top_bc,
                               BottomWaveMesh const&  mesh,
                               Real                   alpha)
{
  int const z_offset  = mesh.grid.offset(2);
  int const global_nz = mesh.grid.global_extent(2);
  auto      coeff     = DiffusionCNZetaCoeffProvider<CenterPt, BottomWaveMesh>(
    mesh, bottom_bc, top_bc, alpha);

  Kokkos::parallel_for(
    LoopPolicy<3>({0, 0, 0}, local_ends(rhs)),
    KOKKOS_LAMBDA(int i, int j, int k) {
      int const  global_k = z_offset + k;
      auto const c        = coeff.coeff(i, j, k);
      Real const u        = reference_solution_value(i, j, global_k);
      Real       v        = c.diag * u;
      if (global_k > 0) {
        v += c.lower * reference_solution_value(i, j, global_k - 1);
      }
      if (global_k + 1 < global_nz) {
        v += c.upper * reference_solution_value(i, j, global_k + 1);
      }
      ref(i, j, k) = u;
      rhs(i, j, k) = v;
    });
  Kokkos::fence();
}

std::pair<double, double> squared_error(MDView<Real***> const& x,
                                        MDView<Real***> const& x_ref)
{
  Kokkos::complex<double> err = 0;
  Kokkos::parallel_reduce(
    LoopPolicy<3>({0, 0, 0}, local_ends(x)),
    KOKKOS_LAMBDA(int i, int j, int k, Kokkos::complex<double>& val) {
      auto ref_val = double(x_ref(i, j, k));
      auto diff    = double(x(i, j, k) - x_ref(i, j, k));
      val += Kokkos::complex{diff * diff, ref_val * ref_val};
    },
    err);
  Kokkos::fence();
  return {err.real(), err.imag()};
}

DiffusionCNZetaSolverOptions
make_strategy_options(DiffusionCNZetaAlgorithm algorithm)
{
  DiffusionCNZetaSolverOptions options;
  options.algorithm = algorithm;
  options.jacobi_max_iters = 400;
  options.jacobi_abs_tol   = 1e-9;
  options.jacobi_rel_tol   = 1e-11;
  return options;
}

std::pair<double, double> run_strategy_case(mpipp::communicator const& comm,
                                            DiffusionCNBCType        bottom_bc,
                                            DiffusionCNBCType        top_bc,
                                            DiffusionCNZetaAlgorithm algorithm)
{
  MeshFixture fixture(comm);
  auto&       mesh = fixture.mesh;

  auto rhs = MDView<Real***, default_memory_pool>(
    Kokkos::view_alloc("rhs", Kokkos::WithoutInitializing),
    mesh.extent(0),
    mesh.extent(1),
    mesh.extent(2));
  auto ref_solution = MDView<Real***, default_memory_pool>(
    Kokkos::view_alloc("reference", Kokkos::WithoutInitializing), rhs.layout());
  auto test_solution =
    MDView<Real***, default_memory_pool>("test", rhs.layout());

  constexpr double alpha = 0.15;
  build_manufactured_system(
    ref_solution, rhs, bottom_bc, top_bc, mesh, (Real)alpha);

  DiffusionCNZetaEqn<CenterPt> eqn(
    alpha, top_bc, bottom_bc, make_strategy_options(algorithm));

  auto stream = Kokkos::DefaultExecutionSpace();
  eqn.solve(test_solution, rhs, mesh, stream);

  auto [local_error, local_norm] = squared_error(test_solution, ref_solution);
  double global_error{}, global_norm{};
  mpipp::allreduce(local_error, global_error, mpipp::plus<double>(), comm);
  mpipp::allreduce(local_norm, global_norm, mpipp::plus<double>(), comm);
  global_error = std::sqrt(global_error / mesh.global_extent(2));
  global_norm  = std::sqrt(global_norm / mesh.global_extent(2));
  return {global_error, global_norm};
}
} // namespace

CATCH_REGISTER_ENUM(alps::solver::DiffusionCNZetaAlgorithm,
                    alps::solver::DiffusionCNZetaAlgorithm::TDMA,
                    alps::solver::DiffusionCNZetaAlgorithm::PCR,
                    alps::solver::DiffusionCNZetaAlgorithm::Jacobi)

CATCH_REGISTER_ENUM(alps::solver::DiffusionCNBCType,
                    alps::solver::DiffusionCNBCType::Dirichlet,
                    alps::solver::DiffusionCNBCType::Neumann)

struct AbsRelToleranceMatcher : Catch::Matchers::MatcherGenericBase
{
  double norm{0};
  double abs_tol{1e-12};
  double rel_tol{1e-10};

  AbsRelToleranceMatcher(double norm_, double abs_tol_, double rel_tol_)
    : norm(norm_)
    , abs_tol(abs_tol_)
    , rel_tol(rel_tol_)
  {}

  bool match(double value) const
  {
    return std::abs(value) <= abs_tol || std::abs(value) / norm <= rel_tol;
  }

  std::string describe() const override
  {
    return fmt::format(
      "is within absolute tolerance {} or relative tolerance {} of {}",
      abs_tol,
      rel_tol,
      norm);
  }
};

TEST_CASE("DiffusionCNZetaEqn solves manufactured CN-zeta system on one rank",
          "[curvilinear][mpi][single-rank]")
{
  using alps::solver::DiffusionCNBCType;
  using alps::solver::DiffusionCNZetaAlgorithm;

  auto comm = mpipp::COMM_WORLD();
  if (comm.size() != 1) {
    SKIP("test is designed to run with 1 rank. Skipping.");
  }

  auto algorithm = GENERATE(DiffusionCNZetaAlgorithm::TDMA,
                            DiffusionCNZetaAlgorithm::PCR,
                            DiffusionCNZetaAlgorithm::Jacobi);
  auto bottom_bc =
    GENERATE(DiffusionCNBCType::Dirichlet, DiffusionCNBCType::Neumann);
  auto top_bc =
    GENERATE(DiffusionCNBCType::Dirichlet, DiffusionCNBCType::Neumann);

  CAPTURE(algorithm, bottom_bc, top_bc);

  auto [error, norm] = run_strategy_case(comm, bottom_bc, top_bc, algorithm);

  constexpr double abs_tol = 1e-12;
  constexpr double rel_tol = std::is_same_v<Real, double> ? 1e-10 : 1e-6;
  REQUIRE_THAT(error, AbsRelToleranceMatcher(norm, abs_tol, rel_tol));
}

TEST_CASE("DiffusionCNZetaEqn solves manufactured CN-zeta system on multiple ranks",
          "[curvilinear][mpi][multi-rank]")
{
  using alps::solver::DiffusionCNBCType;
  using alps::solver::DiffusionCNZetaAlgorithm;

  auto comm = mpipp::COMM_WORLD();
  if (comm.size() < 2) {
    SKIP("test is designed to run with at least 2 ranks. Skipping.");
  }

  auto algorithm = GENERATE(DiffusionCNZetaAlgorithm::TDMA,
                            DiffusionCNZetaAlgorithm::PCR,
                            DiffusionCNZetaAlgorithm::Jacobi);
  auto bottom_bc =
    GENERATE(DiffusionCNBCType::Dirichlet, DiffusionCNBCType::Neumann);
  auto top_bc =
    GENERATE(DiffusionCNBCType::Dirichlet, DiffusionCNBCType::Neumann);

  CAPTURE(algorithm, bottom_bc, top_bc);

  auto [error, norm] = run_strategy_case(comm, bottom_bc, top_bc, algorithm);

  constexpr double abs_tol = 1e-12;
  constexpr double rel_tol = std::is_same_v<Real, double> ? 1e-10 : 1e-6;
  REQUIRE_THAT(error, AbsRelToleranceMatcher(norm, abs_tol, rel_tol));
}
