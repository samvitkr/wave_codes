#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <common/container/view_utils.h>
#include <common/kokkos_abstraction/exec_policy.h>
#include <common/math.h>
#include <common/real_type.h>
#include <decomp/pencil_plan.h>
#include <solvers/poisson/tridiagonal_solver.h>

#include <spdlog/spdlog.h>

#include "poisson_test_data.h"

TEST_CASE("Poisson solver (LU MPI)", "[poisson]")
{
  auto world = mpipp::COMM_WORLD();

  using namespace alps;
  using namespace alps::solver;

  std::vector<int> dims{1, 2, 2};
  std::vector<int> periodic{1, 1, 0};
  std::vector<int> grid{4, 16, 11};

  MPIComm3D  comm(world, dims, periodic);
  PencilPlan decomp(comm, grid);
  REQUIRE(comm.size() == 4);

  auto peqn{create_tridiagonal_solver<Real>(
    4, 8, decomp.extent(2), comm.axis_comm[2], TridiagonalAlgorithm::MpiLU)};

  MDView<Real***> d_g, du_g, dl_g, b_g, x_g;
  get_solution(d_g, du_g, dl_g, b_g, x_g);

  MDView<Real***> d("", create_local_layout(decomp));
  MDView<Real***> du("", d.layout());
  MDView<Real***> dl("", d.layout());
  MDView<Real***> b("", d.layout());
  auto            x_offset = decomp.offsets();
  Kokkos::parallel_for(
    LoopPolicy<3>({0, 0, 0}, local_ends(d)),
    KOKKOS_LAMBDA(int i, int j, int k) {
      d(i, j, k)  = d_g(i, j, k + x_offset[2]);
      du(i, j, k) = du_g(i, j, k + x_offset[2]);
      dl(i, j, k) = dl_g(i, j, k + x_offset[2]);
      b(i, j, k)  = b_g(i, j, k + x_offset[2]);
    });

  peqn->setup(d, dl, du);

  peqn->solve(b, d, dl, du);

  Real err_x = 0;
  Kokkos::parallel_reduce(
    LoopPolicy<3>({0, 0, 0}, local_ends(d)),
    KOKKOS_LAMBDA(int i, int j, int k, Real& err) {
      err += alps::square(x_g(i, j, k + x_offset[2]) - b(i, j, k));
    },
    err_x);
  Kokkos::fence();
  CHECK(err_x == Catch::Approx(0).margin(1e-9));
}

TEST_CASE("Poisson solver (Wang algorithm)", "[poisson]")
{
  auto world = mpipp::COMM_WORLD();

  using namespace alps;
  using namespace alps::solver;

  std::vector<int> dims{1, 2, 2};
  std::vector<int> periodic{1, 1, 0};
  std::vector<int> grid{4, 16, 11};

  MPIComm3D  comm(world, dims, periodic);
  PencilPlan decomp(comm, grid);
  REQUIRE(comm.size() == 4);

  auto peqn{create_tridiagonal_solver<Real>(
    4, 8, decomp.extent(2), comm.axis_comm[2], TridiagonalAlgorithm::Wang)};

  MDView<Real***> d_g, du_g, dl_g, b_g, x_g;
  get_solution(d_g, du_g, dl_g, b_g, x_g);

  MDView<Real***> d("", create_local_layout(decomp));
  MDView<Real***> du("", d.layout());
  MDView<Real***> dl("", d.layout());
  MDView<Real***> b("", d.layout());
  auto            x_offset = decomp.offsets();
  Kokkos::parallel_for(
    LoopPolicy<3>({0, 0, 0}, local_ends(d)),
    KOKKOS_LAMBDA(int i, int j, int k) {
      d(i, j, k)  = d_g(i, j, k + x_offset[2]);
      du(i, j, k) = du_g(i, j, k + x_offset[2]);
      dl(i, j, k) = dl_g(i, j, k + x_offset[2]);
      b(i, j, k)  = b_g(i, j, k + x_offset[2]);
    });

  peqn->setup(d, dl, du);

  // {
  //     auto& solver = dynamic_cast<TridiagonalWang&>(peqn.backend());
  //     auto& d_host = solver.d_gather;
  //     // auto d_host = create_mirror_view(solver.t_);
  //     // deep_copy(d_host, solver.t_);
  //     for (auto k = d_host.begin(2); k < d_host.end(2); ++k)
  //     {
  //       spdlog::info("d_gather {} {} {}", decomp.my_rank(), k, d_host(0, 0,
  //       k));
  //     }
  // }

  peqn->solve(b, d, dl, du);
  // {
  //     auto d_host = create_mirror_view(b);
  //     deep_copy(d_host, b);
  //     for (auto k = d_host.begin(2); k < d_host.end(2); ++k)
  //     {
  //       spdlog::info("x {} {} {}", decomp.my_rank(), k, d_host(0, 0, k));
  //     }
  // }

  Real err_x = 0;
  Kokkos::parallel_reduce(
    LoopPolicy<3>({0, 0, 0}, local_ends(d)),
    KOKKOS_LAMBDA(int i, int j, int k, Real& err) {
      err += alps::square(x_g(i, j, k + x_offset[2]) - b(i, j, k));
    },
    err_x);
  Kokkos::fence();
  CHECK(err_x == Catch::Approx(0).margin(1e-9));
}
