#include <catch2/catch_test_macros.hpp>

#include <decomp/mdcomm.h>
#include <decomp/pencil_plan.h>

#include <type_traits>
#include <vector>

TEST_CASE(" PencilPlan constructor ", "[PencilPlan]")
{
  auto world = mpipp::COMM_WORLD();
  REQUIRE(world.size() == 4);

  std::vector<int> dims      = {1, 2, 2};
  std::vector<int> periodic  = {1, 1, 0};
  std::vector<int> grid_size = {48, 36, 33};

  alps::MPIComm3D  md_comm(world, dims, periodic);
  alps::PencilPlan plan(md_comm, grid_size);

  CHECK(std::is_copy_constructible_v<alps::PencilPlan>);
  CHECK(std::is_move_constructible_v<alps::PencilPlan>);

  SECTION("Global extents")
  {
    REQUIRE(plan.global_extent(0) == 48);
    REQUIRE(plan.global_extent(1) == 36);
    REQUIRE(plan.global_extent(2) == 33);

    REQUIRE(plan.global_extent(0, alps::Pencil::Y) == 36);
    REQUIRE(plan.global_extent(1, alps::Pencil::Y) == 48);
    REQUIRE(plan.global_extent(2, alps::Pencil::Y) == 33);
  }
}

TEST_CASE(" PencilPlan constructor exceptions ", "[PencilPlan]")
{
  auto world = mpipp::COMM_WORLD();
  REQUIRE(world.size() == 4);

  std::vector<int> dims      = {2, 1, 2};
  std::vector<int> periodic  = {0, 0, 1};
  std::vector<int> grid_size = {48, 36};

  alps::MPIComm3D md_comm(world, dims, periodic);
  REQUIRE_THROWS(alps::PencilPlan(md_comm, grid_size));
}
