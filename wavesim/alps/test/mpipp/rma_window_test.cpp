#include <catch2/catch_test_macros.hpp>

#include <mpipp/comm.h>
#include <mpipp/environment.h>
#include <mpipp/rma_window.h>

#include <vector>

TEST_CASE("rma_window basic operations", "[mpipp][win]")
{
  auto world = mpipp::COMM_WORLD();

  SECTION("create window with MPI_Win_create")
  {
    std::vector<int> data(100, world.rank());
    auto             win = mpipp::rma_window<void>(
      data.data(), data.size() * sizeof(int), sizeof(int), world);

    REQUIRE(win.is_valid());
    REQUIRE(win.base_address() == data.data());
    REQUIRE(win.size_bytes() == data.size() * sizeof(int));
    REQUIRE(win.disp_unit() == sizeof(int));
    REQUIRE(win.use_count() == 1);
  }

  SECTION("create window with span (untyped)")
  {
    std::vector<double>  data(50, 3.14);
    nonstd::span<double> span_data{data.data(), data.size()};
    auto                 win = mpipp::rma_window<void>(span_data, world);

    REQUIRE(win.is_valid());
    REQUIRE(win.base_address() == data.data());
    REQUIRE(win.size_bytes() == data.size() * sizeof(double));
    REQUIRE(win.disp_unit() == sizeof(double));
  }

  SECTION("create window with span (typed)")
  {
    std::vector<float>  data(75, 2.71f);
    nonstd::span<float> span_data{data.data(), data.size()};
    auto                win = mpipp::rma_window<float>(span_data, world);

    REQUIRE(win.is_valid());
    REQUIRE(win.base_address() == data.data());
    REQUIRE(win.size_bytes() == data.size() * sizeof(float));
    REQUIRE(win.disp_unit() == sizeof(float));

    // Test local_memory access
    auto local_mem = win.local_memory();
    REQUIRE(local_mem.size() == data.size());
    REQUIRE(local_mem.data() == data.data());
  }

  SECTION("allocate window (typed)")
  {
    auto win = mpipp::rma_window<double>(50, world);

    REQUIRE(win.is_valid());
    REQUIRE(win.base_address() != nullptr);
    REQUIRE(win.size_bytes() == 50 * sizeof(double));
    REQUIRE(win.disp_unit() == sizeof(double));

    // Test local_memory access
    auto local_mem = win.local_memory();
    REQUIRE(local_mem.size() == 50);
    REQUIRE(local_mem.data() == win.base_address());
  }

  SECTION("window copy semantics")
  {
    auto win1 = mpipp::rma_window<int>(10, world);
    REQUIRE(win1.use_count() == 1);

    auto win2 = win1;
    REQUIRE(win1.use_count() == 2);
    REQUIRE(win2.use_count() == 2);
    REQUIRE(win1.raw_handle() == win2.raw_handle());

    auto win3 = win2;
    REQUIRE(win1.use_count() == 3);
    REQUIRE(win2.use_count() == 3);
    REQUIRE(win3.use_count() == 3);
  }

  SECTION("window move semantics")
  {
    auto win1 = mpipp::rma_window<int>(10, world);
    REQUIRE(win1.use_count() == 1);

    auto win2 = std::move(win1);
    REQUIRE(!win1.is_valid());
    REQUIRE(win2.is_valid());
    REQUIRE(win2.use_count() == 1);
  }

  SECTION("window reset")
  {
    auto win = mpipp::rma_window<int>(10, world);
    REQUIRE(win.is_valid());

    win.reset();
    REQUIRE(!win.is_valid());
    REQUIRE(win.use_count() == 0);
  }
}

TEST_CASE("RMA operations", "[mpipp][win][rma]")
{
  auto world = mpipp::COMM_WORLD();

  if (world.size() < 2) {
    SKIP("RMA operations test requires at least 2 MPI ranks");
  }

  SECTION("put and get operations")
  {
    auto win       = mpipp::rma_window<int>(10, world);
    auto local_mem = win.local_memory();

    // Initialize local memory
    for (std::size_t i = 0; i < local_mem.size(); ++i) {
      local_mem[i] = world.rank() * 100 + static_cast<int>(i);
    }

    // Fence to ensure initialization is complete
    win.fence(0);

    // Put operation: rank 0 writes to rank 1
    if (world.rank() == 0) {
      int value = 999;
      mpipp::put(value, 1, 0, win);
    }

    // Fence to complete RMA operations
    win.fence(0);

    // Verify: rank 1 should have received the value
    if (world.rank() == 1) {
      REQUIRE(local_mem[0] == 999);
    }

    // Get operation: rank 1 reads from rank 0
    if (world.rank() == 1) {
      int value = 0;
      mpipp::get(value, 0, 1, win);
      win.fence(0);
      REQUIRE(value == 1); // rank 0's local_mem[1]
    } else {
      win.fence(0);
    }
  }

  SECTION("accumulate operation")
  {
    auto win       = mpipp::rma_window<int>(10, world);
    auto local_mem = win.local_memory();

    // Initialize local memory
    for (std::size_t i = 0; i < local_mem.size(); ++i) {
      local_mem[i] = static_cast<int>(i);
    }

    win.fence(0);

    // All ranks accumulate to rank 0
    int value = world.rank();
    mpipp::accumulate(value, 0, 0, mpipp::sum<int>{}, win);

    win.fence(0);

    // Verify: rank 0 should have the sum
    if (world.rank() == 0) {
      int expected_sum = 0;
      for (int r = 0; r < world.size(); ++r) {
        expected_sum += r;
      }
      REQUIRE(local_mem[0] == expected_sum);
    }
  }

  SECTION("lock/unlock operations")
  {
    auto win       = mpipp::rma_window<int>(10, world);
    auto local_mem = win.local_memory();

    // Initialize local memory
    for (std::size_t i = 0; i < local_mem.size(); ++i) {
      local_mem[i] = world.rank() * 100 + static_cast<int>(i);
    }

    world.barrier();

    // Use passive target synchronization
    if (world.rank() == 0) {
      win.lock(mpipp::lock_type::exclusive, 1, 0);
      int value = 777;
      mpipp::put(value, 1, 0, win);
      win.unlock(1);
    }

    world.barrier();

    // Verify: rank 1 should have received the value
    if (world.rank() == 1) {
      REQUIRE(local_mem[0] == 777);
    }
  }
}

TEST_CASE("Nonblocking RMA operations", "[mpipp][win][rma][nonblocking]")
{
  auto world = mpipp::COMM_WORLD();

  if (world.size() < 2) {
    SKIP("Nonblocking RMA operations test requires at least 2 MPI ranks");
  }
}

TEST_CASE("Atomic RMA operations", "[mpipp][win][rma][atomic]")
{
  auto world = mpipp::COMM_WORLD();

  if (world.size() < 2) {
    SKIP("Atomic RMA operations test requires at least 2 MPI ranks");
  }

  SECTION("fetch_and_op")
  {
    auto win       = mpipp::rma_window<int>(10, world);
    auto local_mem = win.local_memory();

    // Initialize local memory
    for (std::size_t i = 0; i < local_mem.size(); ++i) {
      local_mem[i] = 100;
    }

    win.fence(0);

    // Fetch and add: rank 0 atomically adds to rank 1
    if (world.rank() == 0) {
      int origin = 50;
      int result = 0;
      mpipp::fetch_and_op(origin, result, 1, 0, mpipp::sum<int>{}, win);
      REQUIRE(result == 100); // Original value
    }

    win.fence(0);

    // Verify: rank 1 should have the updated value
    if (world.rank() == 1) {
      REQUIRE(local_mem[0] == 150);
    }
  }

  SECTION("compare_and_swap")
  {
    auto win       = mpipp::rma_window<int>(10, world);
    auto local_mem = win.local_memory();

    // Initialize local memory
    for (std::size_t i = 0; i < local_mem.size(); ++i) {
      local_mem[i] = 200;
    }

    win.fence(0);

    // Compare and swap: rank 0 swaps value on rank 1
    if (world.rank() == 0) {
      int origin  = 300;
      int compare = 200;
      int result  = 0;
      mpipp::compare_and_swap(origin, compare, result, 1, 0, win);
      REQUIRE(result == 200); // Original value
    }

    win.fence(0);

    // Verify: rank 1 should have the new value
    if (world.rank() == 1) {
      REQUIRE(local_mem[0] == 300);
    }
  }

  SECTION("replace and no_op")
  {
    auto win       = mpipp::rma_window<int>(10, world);
    auto local_mem = win.local_memory();

    for (std::size_t i = 0; i < local_mem.size(); ++i) {
      local_mem[i] = 100;
    }

    win.fence(0);

    if (world.rank() == 0) {
      int val_replace = 500;
      mpipp::accumulate(val_replace, 1, 0, mpipp::replace<int>{}, win);

      int val_no_op = 999;
      int res_no_op = 0;
      mpipp::get_accumulate(
        val_no_op, res_no_op, 1, 1, mpipp::no_op<int>{}, win);
    }

    win.fence(0);

    if (world.rank() == 1) {
      REQUIRE(local_mem[0] == 500);
      REQUIRE(local_mem[1] == 100);
    }
  }
}
