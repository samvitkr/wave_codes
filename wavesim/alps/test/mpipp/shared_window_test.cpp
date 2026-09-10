#include <catch2/catch_test_macros.hpp>

#include <mpipp/comm.h>
#include <mpipp/environment.h>
#include <mpipp/shared_window.h>

TEST_CASE("shared_window operations", "[mpipp][win][shared]")
{
  auto world       = mpipp::COMM_WORLD();
  auto shared_comm = world.split_shared();

  SECTION("allocate shared window (untyped)")
  {
    auto win =
      mpipp::shared_window<void>(100 * sizeof(int), sizeof(int), shared_comm);

    REQUIRE(win.is_valid());
    REQUIRE(win.base_address() != nullptr);
    REQUIRE(win.size_bytes() == 100 * sizeof(int));
    REQUIRE(win.disp_unit() == sizeof(int));
  }

  SECTION("allocate shared window (typed)")
  {
    auto win = mpipp::shared_window<double>(50, shared_comm);

    REQUIRE(win.is_valid());
    REQUIRE(win.base_address() != nullptr);
    REQUIRE(win.size_bytes() == 50 * sizeof(double));
    REQUIRE(win.disp_unit() == sizeof(double));

    // Test local_memory access
    auto local_mem = win.local_memory();
    REQUIRE(local_mem.size() == 50);
    REQUIRE(local_mem.data() == win.base_address());
  }

  SECTION("shared_query")
  {
    auto win = mpipp::shared_window<int>(20, shared_comm);

    // Query own memory
    auto own_mem = win.shared_query(shared_comm.rank());
    REQUIRE(own_mem.size() == 20);
    REQUIRE(own_mem.data() == win.base_address());

    // Query other ranks' memory (if available)
    if (shared_comm.size() > 1) {
      int  other_rank = (shared_comm.rank() + 1) % shared_comm.size();
      auto other_mem  = win.shared_query(other_rank);
      REQUIRE(other_mem.size() == 20);
      REQUIRE(other_mem.data() != nullptr);
    }
  }
}
