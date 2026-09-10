#include "result_verifier.h"
#include <common/container/view_types.h>
#include <transpose/transposer_mpi_all2all.h>
#include <transpose/transposer_mpi_p2p.h>
#include <transpose/transposer_mpi_p2pshm.h>

#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>

#include <array>
#include <stdexcept>

#define TYPE_LIST_OPENMP                                  \
  (TransposerMPIAll2All<float, Kokkos::OpenMP>),          \
    (TransposerMPIAll2All<double, Kokkos::OpenMP>),       \
    (TransposerMPIPoint2Point<float, Kokkos::OpenMP>),    \
    (TransposerMPIPoint2Point<double, Kokkos::OpenMP>),   \
    (TransposerMPIPoint2PointSHM<float, Kokkos::OpenMP>), \
    (TransposerMPIPoint2PointSHM<double, Kokkos::OpenMP>)
#define TYPE_LIST_DEV                                                 \
  (TransposerMPIAll2All<float, Kokkos::DefaultExecutionSpace>),       \
    (TransposerMPIAll2All<double, Kokkos::DefaultExecutionSpace>),    \
    (TransposerMPIPoint2Point<float, Kokkos::DefaultExecutionSpace>), \
    (TransposerMPIPoint2Point<double, Kokkos::DefaultExecutionSpace>)
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
#define TYPE_LIST TYPE_LIST_OPENMP, TYPE_LIST_DEV
#else
#define TYPE_LIST TYPE_LIST_OPENMP
#endif

template<typename T>
struct ExtractTemplateParams; // Undefined primary template

// Partial specialization to extract the template parameters for Type<T1, T2>
template<template<typename, typename> class Type, typename T1, typename T2>
struct ExtractTemplateParams<Type<T1, T2>>
{
  using first_type  = T1;
  using second_type = T2;
};

using alps::transpose::TransposerMPIAll2All;
using alps::transpose::TransposerMPIPoint2Point;
using alps::transpose::TransposerMPIPoint2PointSHM;
using alps::transpose::TransposerOptions;

TEMPLATE_TEST_CASE("transpose correctness",
                   "[transpose][distributed]",
                   TYPE_LIST)
{
  using T  = typename ExtractTemplateParams<TestType>::first_type;
  using ES = typename ExtractTemplateParams<TestType>::second_type;

  auto comm = mpipp::COMM_WORLD();
  if (comm.size() != 4) {
    SKIP("Requires 4 MPI processes");
  }

  std::array<int, 3> grid_size{48, 96, 43};
  auto               options = TransposerOptions();
  options.tune_tile_sizes    = false; // Autotune is disabled for testing
  TestType transposer(comm, grid_size[0], grid_size[1], grid_size[2], options);

  TransposeResultVerifier<T, ES> verifier(
    grid_size[0], grid_size[1], grid_size[2], comm.size());
  verifier.initialize_reference();

  auto input = verifier.get_input_view(comm.rank());

  SECTION("TransposeOpAssign")
  {
    using alps::transpose::TransposeOpAssign;

    auto output =
      alps::MDView<T***, typename ES::memory_space>("transposed array",
                                                    grid_size[1],
                                                    grid_size[0] / comm.size(),
                                                    grid_size[2]);

    comm.barrier();
    transposer.execute(output, input, grid_size[2], TransposeOpAssign(), ES());
    ES().fence();
    comm.barrier();
    CHECK(verifier.compare(output, comm.rank(), TransposeOpAssign()));
  }

  SECTION("TransposeOpAdd")
  {
    using alps::transpose::TransposeOpAdd;

    auto output =
      alps::MDView<T***, typename ES::memory_space>("transposed array",
                                                    grid_size[1],
                                                    grid_size[0] / comm.size(),
                                                    grid_size[2]);

    auto constexpr add_value = static_cast<T>(47);
    Kokkos::deep_copy(ES(), output, add_value);
    ES().fence();

    comm.barrier();
    transposer.execute(output, input, grid_size[2], TransposeOpAdd(), ES());
    ES().fence();
    comm.barrier();
    CHECK(verifier.compare(output, comm.rank(), TransposeOpAdd(), add_value));
  }

  SECTION("Partial transpose - half the z-planes")
  {
    using alps::transpose::TransposeOpAssign;
    auto const half_nz = grid_size[2] / 2;

    auto output = alps::MDView<T***, typename ES::memory_space>(
      "partial output", grid_size[1], grid_size[0] / comm.size(), half_nz);

    comm.barrier();
    transposer.execute(output, input, half_nz, TransposeOpAssign(), ES());
    ES().fence();
    comm.barrier();

    CHECK(verifier.compare(output, comm.rank(), TransposeOpAssign()));
  }

  SECTION("Single-plane transpose")
  {
    using alps::transpose::TransposeOpAssign;
    auto constexpr single_nz = 1;

    auto output = alps::MDView<T***, typename ES::memory_space>(
      "single output", grid_size[1], grid_size[0] / comm.size(), single_nz);

    comm.barrier();
    transposer.execute(output, input, single_nz, TransposeOpAssign(), ES());
    ES().fence();
    comm.barrier();

    CHECK(verifier.compare(output, comm.rank(), TransposeOpAssign()));
  }
}

TEMPLATE_TEST_CASE("transpose throws on non-divisible grid sizes",
                   "[transpose][distributed]",
                   TYPE_LIST)
{
  using T  = typename ExtractTemplateParams<TestType>::first_type;
  using ES = typename ExtractTemplateParams<TestType>::second_type;

  auto comm = mpipp::COMM_WORLD();
  if (comm.size() != 4) {
    SKIP("Requires 4 MPI processes");
  }

  auto options            = TransposerOptions();
  options.tune_tile_sizes = false;

  // n0=47 is not divisible by 4 processes
  CHECK_THROWS_AS(TestType(comm, 47, 96, 43, options), std::invalid_argument);

  // n1=97 is not divisible by 4 processes
  CHECK_THROWS_AS(TestType(comm, 48, 97, 43, options), std::invalid_argument);

  // Both 47 and 97 are not divisible by 4 processes
  CHECK_THROWS_AS(TestType(comm, 47, 97, 43, options), std::invalid_argument);
}

TEMPLATE_PRODUCT_TEST_CASE("transpose input parameter validation",
                           "[transpose][distributed]",
                           (TransposerMPIAll2All,
                            TransposerMPIPoint2Point,
                            TransposerMPIPoint2PointSHM),
                           ((float, Kokkos::OpenMP)))
{
  using T  = typename ExtractTemplateParams<TestType>::first_type;
  using ES = typename ExtractTemplateParams<TestType>::second_type;

  auto comm = mpipp::COMM_WORLD();
  if (comm.size() != 4) {
    SKIP("Requires 4 MPI processes");
  }

  std::array<int, 3> grid_size{48, 96, 43};
  auto               options = TransposerOptions();
  options.tune_tile_sizes    = false;
  TestType transposer(comm, grid_size[0], grid_size[1], grid_size[2], options);

  using alps::transpose::TransposeOpAssign;

  auto input = alps::MDView<T***, typename ES::memory_space>(
    "input", grid_size[0], grid_size[1] / comm.size(), grid_size[2]);
  auto output = alps::MDView<T***, typename ES::memory_space>(
    "output", grid_size[1], grid_size[0] / comm.size(), grid_size[2]);

  // Invalid howmany value
  CHECK_THROWS_AS(
    transposer.execute(output, input, -1, TransposeOpAssign(), ES()),
    std::invalid_argument);

  CHECK_THROWS_AS(transposer.execute(
                    output, input, grid_size[2] + 1, TransposeOpAssign(), ES()),
                  std::invalid_argument);

  // Invalid input
  auto smaller_array = alps::MDView<T***, typename ES::memory_space>(
    Kokkos::view_alloc("smaller input", Kokkos::WithoutInitializing),
    input.extent(0) - 1,
    input.extent(1),
    input.extent(2));
  CHECK_THROWS_AS(
    transposer.execute(
      output, smaller_array, grid_size[2], TransposeOpAssign(), ES()),
    std::invalid_argument);

  Kokkos::resize(Kokkos::view_alloc(Kokkos::WithoutInitializing),
                 smaller_array,
                 input.extent(0),
                 input.extent(1) - 1,
                 input.extent(2));
  CHECK_THROWS_AS(
    transposer.execute(
      output, smaller_array, grid_size[2], TransposeOpAssign(), ES()),
    std::invalid_argument);

  // Invalid output
  Kokkos::resize(Kokkos::view_alloc(Kokkos::WithoutInitializing),
                 smaller_array,
                 output.extent(0) - 1,
                 output.extent(1),
                 output.extent(2));
  CHECK_THROWS_AS(
    transposer.execute(
      smaller_array, input, grid_size[2], TransposeOpAssign(), ES()),
    std::invalid_argument);

  Kokkos::resize(Kokkos::view_alloc(Kokkos::WithoutInitializing),
                 smaller_array,
                 output.extent(0),
                 output.extent(1) - 1,
                 output.extent(2));
  CHECK_THROWS_AS(
    transposer.execute(
      smaller_array, input, grid_size[2], TransposeOpAssign(), ES()),
    std::invalid_argument);
}

TEMPLATE_TEST_CASE("transpose with embedded arrays",
                   "[transpose][distributed]",
                   TYPE_LIST)
{
  using T  = typename ExtractTemplateParams<TestType>::first_type;
  using ES = typename ExtractTemplateParams<TestType>::second_type;

  auto comm = mpipp::COMM_WORLD();
  if (comm.size() != 4) {
    SKIP("Requires 4 MPI processes");
  }

  Kokkos::Array<int, 3> grid_size{48, 96, 43};
  auto                  options = TransposerOptions();
  options.tune_tile_sizes       = false; // Autotune is disabled for testing
  TestType transposer(comm, grid_size[0], grid_size[1], grid_size[2], options);

  TransposeResultVerifier<T, ES> verifier(
    grid_size[0], grid_size[1], grid_size[2], comm.size());
  verifier.initialize_reference();

  // Create larger arrays with padding in first two dimensions
  int const pad_i = 8;
  int const pad_j = 12;

  // Local sizes for distributed transpose
  int const local_ny_input  = grid_size[1] / comm.size();
  int const local_nx_output = grid_size[0] / comm.size();

  // Input array: larger in first two dimensions
  auto input_large = alps::MDView<T***, typename ES::memory_space>(
    "input large", grid_size[0] + pad_i, local_ny_input + pad_j, grid_size[2]);

  // Output array: larger in first two dimensions (swapped for transpose)
  auto output_large =
    alps::MDView<T***, typename ES::memory_space>("output large",
                                                  grid_size[1] + pad_j,
                                                  local_nx_output + pad_i,
                                                  grid_size[2]);

  // Initialize the large input array with a pattern
  Kokkos::deep_copy(input_large, static_cast<T>(-1));

  // Copy the actual input data into the embedded region
  auto input = verifier.get_input_view(comm.rank());
  Kokkos::deep_copy(ES(),
                    subview(input_large,
                            std::make_pair(0, grid_size[0]),
                            std::make_pair(0, local_ny_input),
                            std::make_pair(0, grid_size[2])),
                    input);
  ES().fence();

  // Initialize output with a marker value to verify only embedded region is
  // modified
  T const marker = static_cast<T>(-999);
  Kokkos::deep_copy(ES(), output_large, marker);
  ES().fence();

  SECTION("Embedded input and output arrays")
  {
    using alps::transpose::TransposeOpAssign;
    comm.barrier();
    transposer.execute(
      output_large, input_large, grid_size[2], TransposeOpAssign(), ES());
    ES().fence();
    comm.barrier();

    // Verify the transposed result in the embedded region
    CHECK(verifier.compare(output_large, comm.rank(), TransposeOpAssign()));

    // Verify that padding regions in output are unchanged
    bool padding_unchanged = true;
    Kokkos::parallel_reduce(
      "check output padding",
      alps::LoopPolicy<3, ES>(
        ES(), alps::begins(output_large), alps::ends(output_large)),
      KOKKOS_LAMBDA(int i, int j, int k, bool& unchanged) {
        // Check if this is in the padding region
        if (i >= grid_size[1] || j >= local_nx_output) {
          if (output_large(i, j, k) != marker) {
            unchanged = false;
          }
        }
      },
      Kokkos::LAnd<bool>(padding_unchanged));
    ES().fence();
    CHECK(padding_unchanged);
  }
}
