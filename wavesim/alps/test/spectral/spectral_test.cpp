#include <catch2/catch_approx.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#define _USE_MATH_DEFINES

#include <Kokkos_Random.hpp>
#include <cmath>
#include <common/container/view_types.h>
#include <common/container/view_utils.h>
#include <common/kokkos_abstraction/exec_policy.h>
#include <limits>
#include <spectral/spectral_base.h>

#if defined(KOKKOS_ENABLE_CUDA)
#define SIG_LIST                                  \
  (float, Kokkos::OpenMP, alps::fft::FFTW, 1),    \
    (double, Kokkos::OpenMP, alps::fft::FFTW, 1), \
    (float, Kokkos::Cuda, alps::fft::VKFFT, 1),   \
    (double, Kokkos::Cuda, alps::fft::VKFFT, 1),  \
    (float, Kokkos::Cuda, alps::fft::CUFFT, 1),   \
    (double, Kokkos::Cuda, alps::fft::CUFFT, 1)
#elif defined(KOKKOS_ENABLE_HIP)
#define SIG_LIST                                  \
  (float, Kokkos::OpenMP, alps::fft::FFTW, 1),    \
    (double, Kokkos::OpenMP, alps::fft::FFTW, 1), \
    (float, Kokkos::HIP, alps::fft::VKFFT, 1),    \
    (double, Kokkos::HIP, alps::fft::VKFFT, 1)
#else
#define SIG_LIST                               \
  (float, Kokkos::OpenMP, alps::fft::FFTW, 1), \
    (double, Kokkos::OpenMP, alps::fft::FFTW, 1)
#endif

TEMPLATE_TEST_CASE_SIG(
  "Spectral construction",
  "[FFT]",
  ((typename T, typename ES, typename Backend, int V), T, ES, Backend, V),
  SIG_LIST)
{
  auto world = mpipp::COMM_WORLD();
  REQUIRE(world.size() == 4);

  std::vector<int> dims{1, 2, 2};
  std::vector<int> periodic{1, 1, 0};
  std::vector<int> grid_size{32, 16, 8};

  using complex = Kokkos::complex<T>;
  alps::MPIComm3D  md_comm(world, dims, periodic);
  alps::PencilPlan decomp(md_comm, grid_size);

  T const pex   = 0.5;
  T const pey   = 3.0;
  auto spectral = alps::spectral::SpectralPlanFactory::create<T, ES, Backend>(
    decomp, pex, pey);

  using MDRPolicy3 = alps::LoopPolicy<3, ES>;

  SECTION("x derivative")
  {
    auto x_offset1 = decomp.offsets()[1];

    alps::MDView<T***, ES> input("original", alps::create_local_layout(decomp));
    alps::MDView<T***, ES> result("result", input.layout());
    alps::MDView<T***, ES> output("output", input.layout());

    auto       nx = grid_size[0];
    auto       ny = grid_size[1];
    MDRPolicy3 range_policy(alps::begins(input), alps::ends(input));
    Kokkos::parallel_for(
      range_policy, KOKKOS_LAMBDA(int i, int j, int k) {
        T x            = i * 2 * M_PI / pex / nx;
        T y            = (x_offset1 + j) * 2 * M_PI / pey / ny;
        input(i, j, k) = k
                       * (std::sin(3 * pex * x + 2 * pey * y)
                          + 0.2 * std::cos(13 * pex * x + 3.0));
        result(i, j, k) = k
                        * (3 * pex * std::cos(3 * pex * x + 2 * pey * y)
                           - 0.2 * 13 * pex * std::sin(13 * pex * x + 3.0));
      });

    spectral->do_ddx(output, input, ES());

    double err = 0;
    Kokkos::parallel_reduce(
      range_policy,
      KOKKOS_LAMBDA(int i, int j, int k, double& s) {
        s += std::abs(output(i, j, k) - result(i, j, k));
      },
      err);
    double norm_result = 0;
    Kokkos::parallel_reduce(
      range_policy,
      KOKKOS_LAMBDA(int i, int j, int k, double& s) {
        s += std::abs(result(i, j, k));
      },
      norm_result);
    Kokkos::fence();

    REQUIRE(err / norm_result
            == Catch::Approx(0).margin(std::numeric_limits<T>::epsilon()
                                       * std::log(grid_size[0]) * 6));
  }

  SECTION("y derivative with transpose")
  {
    auto x_offset0 = decomp.offset(0);
    auto x_offset1 = decomp.offset(1);

    alps::MDView<T***, ES> input("original", alps::create_local_layout(decomp));
    alps::MDView<T***, ES> result("result", alps::create_local_layout(decomp));
    alps::MDView<T***, ES> output("output", alps::create_local_layout(decomp));

    auto       nx = grid_size[0];
    auto       ny = grid_size[1];
    MDRPolicy3 range_policy(alps::begins(input), alps::ends(input));
    Kokkos::parallel_for(
      range_policy, KOKKOS_LAMBDA(int i, int j, int k) {
        T x            = (x_offset0 + i) * 2 * M_PI / pex / nx;
        T y            = (x_offset1 + j) * 2 * M_PI / pey / ny;
        input(i, j, k) = k
                       * (std::sin(3 * pex * x - 2 * pey * y)
                          + 0.2 * std::cos(pey * y + 3.0));
        result(i, j, k) = k
                        * (-2 * pey * std::cos(3 * pex * x - 2 * pey * y)
                           - 0.2 * pey * std::sin(pey * y + 3.0));
      });

    spectral->do_ddy(output,
                     input,
                     alps::spectral::SpectralPostOp::AssignAfterTranspose,
                     ES());

    double err = 0;
    Kokkos::parallel_reduce(
      range_policy,
      KOKKOS_LAMBDA(int i, int j, int k, double& s) {
        s += std::abs(output(i, j, k) - result(i, j, k));
      },
      err);
    double norm_result = 0;
    Kokkos::parallel_reduce(
      range_policy,
      KOKKOS_LAMBDA(int i, int j, int k, double& s) {
        s += std::abs(result(i, j, k));
      },
      norm_result);
    Kokkos::fence();

    REQUIRE(err / norm_result
            == Catch::Approx(0).margin(std::numeric_limits<T>::epsilon()
                                       * std::log(grid_size[1]) * 6));
  }
}
