#include <catch2/catch_test_macros.hpp>

#include <decomp/block_partition.h>
#include <decomp/mdcomm.h>

#include <type_traits>
#include <vector>

TEST_CASE(" BlockPartition constructor ", "[BlockPartition]")
{
  auto world = mpipp::COMM_WORLD();
  REQUIRE(world.size() == 4);

  std::vector<int> dims      = {2, 1, 2};
  std::vector<int> periodic  = {0, 0, 1};
  std::vector<int> grid_size = {48, 36, 33};

  alps::MPIComm3D      md_comm(world, dims, periodic);
  alps::BlockPartition partition(md_comm, grid_size);

  CHECK(std::is_copy_constructible_v<alps::BlockPartition>);

  SECTION("Global extents")
  {
    REQUIRE(partition.global_extents[0] == grid_size[0]);
    REQUIRE(partition.global_extents[1] == grid_size[1]);
    REQUIRE(partition.global_extents[2] == grid_size[2]);
  }

  SECTION("Local extents")
  {
    // check the rank-coordinate mapping
    auto rank = md_comm.rank();

    switch (rank) {
      case 0:
        REQUIRE(partition.extents[0] == 24);
        REQUIRE(partition.extents[1] == 36);
        REQUIRE(partition.extents[2] == 17);
        break;

      case 1:
        REQUIRE(partition.extents[0] == 24);
        REQUIRE(partition.extents[1] == 36);
        REQUIRE(partition.extents[2] == 17);
        break;

      case 2:
        REQUIRE(partition.extents[0] == 24);
        REQUIRE(partition.extents[1] == 36);
        REQUIRE(partition.extents[2] == 16);
        break;

      case 3:
        REQUIRE(partition.extents[0] == 24);
        REQUIRE(partition.extents[1] == 36);
        REQUIRE(partition.extents[2] == 16);
        break;
    }
  }

  SECTION("Offsets")
  {
    // check the rank-coordinate mapping
    auto rank = md_comm.rank();

    switch (rank) {
      case 0:
        REQUIRE(partition.offsets[0] == 0);
        REQUIRE(partition.offsets[1] == 0);
        REQUIRE(partition.offsets[2] == 0);
        break;

      case 1:
        REQUIRE(partition.offsets[0] == 24);
        REQUIRE(partition.offsets[1] == 0);
        REQUIRE(partition.offsets[2] == 0);
        break;

      case 2:
        REQUIRE(partition.offsets[0] == 0);
        REQUIRE(partition.offsets[1] == 0);
        REQUIRE(partition.offsets[2] == 17);
        break;

      case 3:
        REQUIRE(partition.offsets[0] == 24);
        REQUIRE(partition.offsets[1] == 0);
        REQUIRE(partition.offsets[2] == 17);
        break;
    }
  }
}

TEST_CASE(" BlockPartition constructor exceptions ", "[BlockPartition]")
{
  auto world = mpipp::COMM_WORLD();
  REQUIRE(world.size() == 4);

  std::vector<int> dims      = {2, 1, 2};
  std::vector<int> periodic  = {0, 0, 1};
  std::vector<int> grid_size = {48, 36};

  alps::MPIComm3D md_comm(world, dims, periodic);
  REQUIRE_THROWS(alps::BlockPartition(md_comm, grid_size));
}
