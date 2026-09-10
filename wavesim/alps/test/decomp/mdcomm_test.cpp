#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <decomp/mdcomm.h>
#include <vector>

TEST_CASE(" MDComm constructor with comm, dims, and periodic ", "[MDComm]")
{
  auto world = mpipp::COMM_WORLD();
  REQUIRE(world.size() == 4);

  SECTION("Cartesian topology constructor")
  {
    std::vector<int> dims     = {2, 1, 2};
    std::vector<int> periodic = {0, 0, 1};

    alps::MPIComm3D md_comm(world, dims, periodic);

    // check the dimensions and periodic information of the processor grid
    for (int dim = 0; dim < dims.size(); ++dim) {
      REQUIRE(md_comm.dims[dim] == dims[dim]);
    }

    // check the rank-coordinate mapping
    auto rank   = md_comm.rank();
    auto coords = md_comm.to_coordinates(rank);
    switch (rank) {
      case 0:
        REQUIRE(coords[0] == 0);
        REQUIRE(coords[1] == 0);
        REQUIRE(coords[2] == 0);
        break;

      case 1:
        REQUIRE(coords[0] == 1);
        REQUIRE(coords[1] == 0);
        REQUIRE(coords[2] == 0);
        break;

      case 2:
        REQUIRE(coords[0] == 0);
        REQUIRE(coords[1] == 0);
        REQUIRE(coords[2] == 1);
        break;

      case 3:
        REQUIRE(coords[0] == 1);
        REQUIRE(coords[1] == 0);
        REQUIRE(coords[2] == 1);
        break;
    }
  }

  SECTION("Cartesian topology constructor incorrect argument 1")
  {
    std::vector<int> dims     = {2, 1};
    std::vector<int> periodic = {0, 0, 1};

    REQUIRE_THROWS_WITH(alps::MPIComm3D(world, dims, periodic),
                        Catch::Matchers::ContainsSubstring("incorrect number "
                                                           "of axes"));
  }

  SECTION("Cartesian topology constructor incorrect argument 2")
  {
    std::vector<int> dims     = {2, 1, 3};
    std::vector<int> periodic = {0, 0, 1};

    REQUIRE_THROWS_WITH(alps::MPIComm3D(world, dims, periodic),
                        Catch::Matchers::ContainsSubstring("cannot build a "
                                                           "grid"));
  }
}

TEST_CASE(" Neighbors ", "[MDComm]")
{
  auto world = mpipp::COMM_WORLD();
  REQUIRE(world.size() == 4);

  std::vector<int> dims     = {2, 2, 1};
  std::vector<int> periodic = {1, 0, 0};

  alps::MPIComm3D md_comm(world, dims, periodic);

  switch (md_comm.rank()) {
    case 0:
      REQUIRE(md_comm.next_proc_on_axis[0] == 1);
      REQUIRE(md_comm.next_proc_on_axis[1] == 2);
      REQUIRE(md_comm.next_proc_on_axis[2] == MPI_PROC_NULL);
      REQUIRE(md_comm.prev_proc_on_axis[0] == 1);
      REQUIRE(md_comm.prev_proc_on_axis[1] == MPI_PROC_NULL);
      REQUIRE(md_comm.prev_proc_on_axis[2] == MPI_PROC_NULL);
      break;

    case 1:
      REQUIRE(md_comm.next_proc_on_axis[0] == 0);
      REQUIRE(md_comm.next_proc_on_axis[1] == 3);
      REQUIRE(md_comm.next_proc_on_axis[2] == MPI_PROC_NULL);
      REQUIRE(md_comm.prev_proc_on_axis[0] == 0);
      REQUIRE(md_comm.prev_proc_on_axis[1] == MPI_PROC_NULL);
      REQUIRE(md_comm.prev_proc_on_axis[2] == MPI_PROC_NULL);
      break;

    case 2:
      REQUIRE(md_comm.next_proc_on_axis[0] == 3);
      REQUIRE(md_comm.next_proc_on_axis[1] == MPI_PROC_NULL);
      REQUIRE(md_comm.next_proc_on_axis[2] == MPI_PROC_NULL);
      REQUIRE(md_comm.prev_proc_on_axis[0] == 3);
      REQUIRE(md_comm.prev_proc_on_axis[1] == 0);
      REQUIRE(md_comm.prev_proc_on_axis[2] == MPI_PROC_NULL);
      break;

    case 3:
      REQUIRE(md_comm.next_proc_on_axis[0] == 2);
      REQUIRE(md_comm.next_proc_on_axis[1] == MPI_PROC_NULL);
      REQUIRE(md_comm.next_proc_on_axis[2] == MPI_PROC_NULL);
      REQUIRE(md_comm.prev_proc_on_axis[0] == 2);
      REQUIRE(md_comm.prev_proc_on_axis[1] == 1);
      REQUIRE(md_comm.prev_proc_on_axis[2] == MPI_PROC_NULL);
      break;
  }
}

TEST_CASE(" Axis communicator ", "[MDComm]")
{
  auto world = mpipp::COMM_WORLD();
  REQUIRE(world.size() == 4);

  std::vector<int> dims     = {2, 2, 1};
  std::vector<int> periodic = {1, 0, 0};

  alps::MPIComm3D md_comm(world, dims, periodic);

  for (int dim = 0; dim < dims.size(); ++dim) {
    REQUIRE(md_comm.axis_comm[dim].size() == dims[dim]);
    REQUIRE(md_comm.axis_comm[dim].rank() == md_comm.coords[dim]);
  }
}
