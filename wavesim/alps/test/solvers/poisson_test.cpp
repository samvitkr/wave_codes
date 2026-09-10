#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <common/kokkos_abstraction/exec_policy.h>
#include <common/math.h>
#include <common/real_type.h>
#include <decomp/mdcomm.h>
#include <solvers/poisson/tridiagonal_solver.h>

#include <spdlog/spdlog.h>

#include "poisson_test_data.h"

TEST_CASE("Poisson solver (single LU)", "[poisson]")
{
  auto world = mpipp::COMM_WORLD();

  using namespace alps;
  using namespace alps::solver;

  std::vector<int> dims{1, 1, 1};
  std::vector<int> periodic{1, 1, 0};
  std::vector<int> grid{1, 1, 11};

  MPIComm3D comm(world, dims, periodic);
  REQUIRE(comm.size() == 1);

  auto peqn{create_tridiagonal_solver<Real>(
    4, 8, grid[2], comm.comm, TridiagonalAlgorithm::SeqLU)};

  MDView<Real***> d, du, dl, b, x;
  MDView<Real***> lu_d, lu_du, lu_dl;

  get_solution(d, du, dl, b, x);
  peqn->setup(d, dl, du);

  SECTION("LU decomposition")
  {
    get_LU_result(lu_d, lu_du, lu_dl);

    Real err_d = 0;
    Kokkos::parallel_reduce(
      LoopPolicy<3>({0, 0, 0}, {4, 8, grid[2]}),
      KOKKOS_LAMBDA(int i, int j, int k, Real& err) {
        err += alps::square(d(i, j, k) - lu_d(i, j, k));
      },
      err_d);
    Kokkos::fence();
    CHECK(err_d == Catch::Approx(0).margin(1e-12));

    Real err_du = 0;
    Kokkos::parallel_reduce(
      LoopPolicy<3>({0, 0, 0}, {4, 8, grid[2]}),
      KOKKOS_LAMBDA(int i, int j, int k, Real& err) {
        err += alps::square(du(i, j, k) - lu_du(i, j, k));
      },
      err_du);
    Kokkos::fence();
    CHECK(err_du == Catch::Approx(0).margin(1e-12));

    Real err_dl = 0;
    Kokkos::parallel_reduce(
      LoopPolicy<3>({0, 0, 0}, {4, 8, grid[2]}),
      KOKKOS_LAMBDA(int i, int j, int k, Real& err) {
        err += alps::square(dl(i, j, k) - lu_dl(i, j, k));
      },
      err_dl);
    Kokkos::fence();
    CHECK(err_dl == Catch::Approx(0).margin(1e-12));
  }

  SECTION("Solution")
  {
    peqn->solve(b, d, dl, du);

    Real err_x = 0;
    Kokkos::parallel_reduce(
      LoopPolicy<3>({0, 0, 0}, {4, 8, grid[2]}),
      KOKKOS_LAMBDA(int i, int j, int k, Real& err) {
        err += alps::square(x(i, j, k) - b(i, j, k));
      },
      err_x);
    Kokkos::fence();
    CHECK(err_x == Catch::Approx(0).margin(1e-9));
  }
}

#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
TEST_CASE("Poisson solver (cusparse)", "[poisson]")
{
  auto world = mpipp::COMM_WORLD();

  using namespace alps;
  using namespace alps::solver;

  std::vector<int> dims{1, 1, 1};
  std::vector<int> periodic{1, 1, 0};
  std::vector<int> grid{1, 1, 11};

  MPIComm3D comm(world, dims, periodic);
  REQUIRE(comm.size() == 1);

  auto peqn{create_tridiagonal_solver<Real>(
    4, 8, grid[2], comm, TridiagonalAlgorithm::CuSparse)};

  HaloView<Real***> d, du, dl, b, x;
  HaloView<Real***> lu_d, lu_du, lu_dl;

  get_solution(d, du, dl, b, x);
  peqn->setup(d, dl, du);

  peqn->solve(b, d, dl, du);

  Real err_x = 0;
  Kokkos::parallel_reduce(
    LoopPolicy<3>({0, 0, 0}, {4, 8, grid[2]}),
    KOKKOS_LAMBDA(int i, int j, int k, Real& err) {
      err += alps::square(x(i, j, k) - b(i, j, k));
    },
    err_x);
  Kokkos::fence();
  CHECK(err_x == Catch::Approx(0).margin(1e-9));
}
#endif
