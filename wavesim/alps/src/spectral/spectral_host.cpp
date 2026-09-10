#include "spectral_host.h"

#include <common/base/logging.h>
#include <common/kokkos_abstraction/pool_space.h>
#include <decomp/pencil_plan.h>
#include <decomp/pencil_transpose.h>

#include <Kokkos_Core.hpp>
#include <omp.h>

#include <array>
#include <string_view>
#include <utility>

namespace alps::spectral {

using Kokkos::WithoutInitializing;
template<class DataType, class... Props>
using buffer_t = Kokkos::View<DataType,
                              Kokkos::LayoutLeft,
                              ::alps::memory_pool<Kokkos::HostSpace>,
                              Props...>;

namespace {
/** Return the smallest multiple of a factor no less than n
 */
template<typename iType1, typename iType2>
constexpr int nextMultiple(iType1 n, iType2 factor)
{
  return static_cast<int>(n + (factor - n % factor) % factor);
}

constexpr std::size_t ComplexBufferAlignment = Kokkos::Impl::MEMORY_ALIGNMENT;
template<class RealT>
constexpr std::size_t ComplexBufferMultiples =
  ComplexBufferAlignment / sizeof(std::complex<RealT>);

template<typename RealT, typename iType>
constexpr auto nextMultipleComplex(iType n)
{
  return nextMultiple(n, ComplexBufferMultiples<RealT>);
}

template<class T>
void omp_normalize(T* ptr_cmplx, int n);

template<class T>
void omp_pad_and_normalize(T* ptr_cmplx, int n, int kc);
} // namespace

template<class T>
SpectralPlan<T, Kokkos::OpenMP, fft::FFTW>::SpectralPlan(
  PencilPlan const& partition,
  T                 kx0,
  T                 ky0)
  : base_t{partition, kx0, ky0, fft::FFTW()}
  , plan_r2c_x(nullptr, &fft::FFTWInterface<T>::destroy_plan)
  , plan_r2c_y(nullptr, &fft::FFTWInterface<T>::destroy_plan)
  , plan_c2r_x(nullptr, &fft::FFTWInterface<T>::destroy_plan)
  , plan_c2r_y(nullptr, &fft::FFTWInterface<T>::destroy_plan)
{
  auto create_1d_plans = [](plan_t& r2c, plan_t& c2r, int const dim) {
    buffer_t<T**>        real_array("tmp real", dim, 2);
    buffer_t<ComplexT**> complex_array(
      "tmp complex", nextMultipleComplex<T>(dim / 2 + 1), 2);
    auto* ptr_real  = &(real_array(0, 1));
    auto* ptr_cmplx = (fftw_complex_t*)(&(complex_array(0, 1))); // NOLINT

    std::array<int, 1> dims = {dim};

    // plan creation uses real_array(0, 1) to account for batch alignment
    // R2C transform is set to preserve the input
    r2c.reset(fft::FFTWInterface<T>::plan_r2c(
      1, dims.data(), ptr_real, ptr_cmplx, FFTW_PRESERVE_INPUT));
    c2r.reset(
      fft::FFTWInterface<T>::plan_c2r(1, dims.data(), ptr_cmplx, ptr_real, 0));
  };
  create_1d_plans(plan_r2c_x, plan_c2r_x, partition.global_extent(0));
  create_1d_plans(plan_r2c_y, plan_c2r_y, partition.global_extent(1));

  tune_transpose_and_transform_chunk_size();

  using namespace std::string_view_literals;
  logger_->debug("Created a SpectralPlan ({}, {}, {})",
                 std::is_same_v<T, double> ? "double"sv : "float"sv,
                 ExecSpace::name(),
                 "FFTW");
}

template<class T>
SpectralPlan<T, Kokkos::OpenMP, fft::FFTW>::~SpectralPlan()
{
  using namespace std::string_view_literals;
  logger_->trace("Destroy SpectralPlan ({}, {}, {})",
                 std::is_same_v<T, double> ? "double"sv : "float"sv,
                 ExecSpace::name(),
                 "FFTW");
}

template<class T>
Kokkos::LayoutLeft
SpectralPlan<T, Kokkos::OpenMP, fft::FFTW>::get_r2c_xy_output_layout() const
{
  return Kokkos::LayoutLeft(
    nextMultipleComplex<T>(pencil.extent(0, Pencil::Y) / 2 + 1) * 2,
    pencil.extent(1, Pencil::Y),
    pencil.extent(2, Pencil::Y));
}

namespace {
template<class T>
void mul_ik(Kokkos::complex<T>* values, int transform_len, T k0)
{
  static_assert(std::is_floating_point_v<T>);
  const auto len_half      = transform_len / 2;
  const auto normalized_k0 = k0 / transform_len;
  auto*      cur           = values;
  for (int i = 0; i < len_half; ++i) {
    auto v = i * normalized_k0;
    *cur   = Kokkos::complex<T>(cur->imag() * (-v), cur->real() * v);
    ++cur;
  }
  values[len_half] = 0;
}

template<class T>
void mul_ik2(Kokkos::complex<T>* values, int transform_len, T k0)
{
  static_assert(std::is_floating_point_v<T>);
  const auto len_half = transform_len / 2;
  auto*      cur      = values;
  for (int i = 0; i < len_half; ++i) {
    auto v = -(i * k0) * (i * k0) / transform_len;
    *cur   = Kokkos::complex<T>(cur->real() * v, cur->imag() * v);
    ++cur;
  }
  values[len_half] = 0;
}
} // namespace

// function template for computing derivatives w.r.t. x
template<class T>
template<class Op>
void SpectralPlan<T, Kokkos::OpenMP, fft::FFTW>::do_x_derivatives(
  view_arg_t const&       output,
  const_view_arg_t const& input,
  Op                      derivative_func,
  ExecSpace const& /*space*/) const
{
  const auto nz = input.extent_int(2);
  const auto ny = input.extent_int(1);
  const auto nx = pencil.extent(0);

  const std::size_t stride    = nextMultipleComplex<T>(nx / 2 + 1);
  const auto        n_threads = omp_get_max_threads();

  auto  mem_space      = default_host_memory_pool();
  auto* complex_buffer = static_cast<ComplexT*>(
    mem_space.allocate(n_threads * stride * sizeof(ComplexT)));

#pragma omp parallel
  {
    auto  thread_id = omp_get_thread_num();
    auto* ptr_cmplx = complex_buffer + thread_id * stride;
#pragma omp for collapse(2)
    for (int k = 0; k < nz; ++k) {
      for (int j = 0; j < ny; ++j) {
        interface::exec_r2c(plan_r2c_x.get(),
                            const_cast<T*>(&(input(0, j, k))),
                            (fftw_complex_t*)ptr_cmplx);

        derivative_func(ptr_cmplx, nx, pex); // multiply ik

        interface::exec_c2r(
          plan_c2r_x.get(), (fftw_complex_t*)ptr_cmplx, &(output(0, j, k)));
      }
    }
  }

  mem_space.deallocate(complex_buffer, n_threads * stride * sizeof(ComplexT));
}

template<class T>
void SpectralPlan<T, Kokkos::OpenMP, fft::FFTW>::do_ddx(
  view_arg_t const&       output,
  const_view_arg_t const& input,
  ExecSpace const&        space) const
{
  logger_->trace("Calculate ddx on {} and write to {} ({}x{}x{})",
                 input.label(),
                 output.label(),
                 input.extent(0),
                 input.extent(1),
                 input.extent(2));

  do_x_derivatives(output, input, mul_ik<T>, space);
}

template<class T>
void SpectralPlan<T, Kokkos::OpenMP, fft::FFTW>::do_d2dx2(
  view_arg_t const&       output,
  const_view_arg_t const& input,
  ExecSpace const&        space) const
{
  logger_->trace("Calculate d2dx2 on {} and write to {} ({}x{}x{})",
                 input.label(),
                 output.label(),
                 input.extent(0),
                 input.extent(1),
                 input.extent(2));

  do_x_derivatives(output, input, mul_ik2<T>, space);
}

// function template for computing derivatives w.r.t. y
template<class T>
template<class Op>
void SpectralPlan<T, Kokkos::OpenMP, fft::FFTW>::do_y_derivatives(
  view_arg_t const&       output,
  const_view_arg_t const& input,
  Op                      derivative_func,
  SpectralPostOp          output_process,
  ExecSpace const&        space) const
{
  using Kokkos::ALL;
  using Kokkos::make_pair;
  using Kokkos::subview;
  using std::min;

  const auto nz = input.extent_int(2);
  const auto ny = pencil.extent(0, Pencil::Y);

  const std::size_t stride    = nextMultipleComplex<T>(ny / 2 + 1);
  const auto        n_threads = omp_get_max_threads();

  auto        mem_space      = default_host_memory_pool();
  auto* const complex_buffer = static_cast<ComplexT*>(
    mem_space.allocate(n_threads * stride * sizeof(ComplexT)));

  // Perform transpose and transform in chunks in the z-direction
  // (improving cache locality)
  const int nz_tile = (pencil.comm().dims[1] > 1) ? nz
                    : (transpose_and_transform_chunk_size_z < 1)
                      ? nz
                      : transpose_and_transform_chunk_size_z;

  buffer_t<T***> buf_y(Kokkos::view_alloc("transpose", WithoutInitializing),
                       ny,
                       pencil.extent(1, Pencil::Y),
                       nz_tile);
  for (int z_begin = 0; z_begin < nz; z_begin += nz_tile) {
    const auto z_end = min(z_begin + nz_tile, nz);
    transpose_xy(buf_y,
                 subview(input, ALL, ALL, make_pair(z_begin, z_end)),
                 pencil,
                 space);

#pragma omp parallel
    {
      auto  thread_id = omp_get_thread_num();
      auto* ptr_cmplx = complex_buffer + thread_id * stride;
#pragma omp for collapse(2)
      for (int k = 0; k < buf_y.extent_int(2); ++k) {
        for (int j = 0; j < buf_y.extent_int(1); ++j) {
          auto* ptr_input = &(buf_y(0, j, k));
          interface::exec_r2c(
            plan_r2c_y.get(), ptr_input, (fftw_complex_t*)ptr_cmplx);

          derivative_func(ptr_cmplx, ny, pey); // multiply ik

          interface::exec_c2r(
            plan_c2r_y.get(), (fftw_complex_t*)ptr_cmplx, ptr_input);
        }
      }
    }

    if (output_process == SpectralPostOp::AssignAfterTranspose) {
      transpose_yx(subview(output, ALL, ALL, make_pair(z_begin, z_end)),
                   buf_y,
                   pencil,
                   space);
    } else if (output_process == SpectralPostOp::AddAfterTranspose) {
      transpose_yx_and_add(subview(output, ALL, ALL, make_pair(z_begin, z_end)),
                           buf_y,
                           pencil,
                           space);
    }
  }

  mem_space.deallocate(complex_buffer, n_threads * stride * sizeof(ComplexT));
}

template<class T>
void SpectralPlan<T, Kokkos::OpenMP, fft::FFTW>::do_ddy(
  view_arg_t const&       output,
  const_view_arg_t const& input,
  SpectralPostOp          output_process,
  ExecSpace const&        space) const
{
  logger_->trace(
    "Calculate ddy on {} and {} to {} ({}x{}x{})",
    input.label(),
    output_process == SpectralPostOp::AssignAfterTranspose ? "write" : "add",
    output.label(),
    input.extent(0),
    input.extent(1),
    input.extent(2));

  do_y_derivatives(output, input, mul_ik<T>, output_process, space);
}

template<class T>
void SpectralPlan<T, Kokkos::OpenMP, fft::FFTW>::do_d2dy2(
  view_arg_t const&       output,
  const_view_arg_t const& input,
  SpectralPostOp          output_process,
  ExecSpace const&        space) const
{
  logger_->trace(
    "Calculate d2dy2 on {} and {} to {} ({}x{}x{})",
    input.label(),
    output_process == SpectralPostOp::AssignAfterTranspose ? "write" : "add",
    output.label(),
    input.extent(0),
    input.extent(1),
    input.extent(2));

  do_y_derivatives(output, input, mul_ik2<T>, output_process, space);
}

template<class T>
void SpectralPlan<T, Kokkos::OpenMP, fft::FFTW>::do_cutoff_xy(
  view_arg_t const& input,
  int               kc_x,
  int               kc_y,
  ExecSpace const&  space) const
{
  using Kokkos::ALL;
  using Kokkos::subview;
  using index_range = std::pair<int, int>;

  logger_->trace("Cutoff {} ({}x{}x{}) to wavenumber kx={} ky={}",
                 input.label(),
                 input.extent(0),
                 input.extent(1),
                 input.extent(2),
                 kc_x,
                 kc_y);

  const auto nz = input.extent_int(2);
  const auto ny = input.extent_int(1);

  // Perform transpose and transform in chunks in the z-direction
  // (improving cache locality)
  const int nz_tile = (pencil.comm().dims[1] > 1) ? nz
                    : (transpose_and_transform_chunk_size_z < 1)
                      ? nz
                      : transpose_and_transform_chunk_size_z;

  // Allocate temporary buffers for transforms and transpose
  const auto x_stride = nextMultipleComplex<T>(pencil.extent(0) / 2 + 1);
  const auto y_stride =
    nextMultipleComplex<T>(pencil.extent(0, Pencil::Y) / 2 + 1);

  buffer_t<T***> buf_x(Kokkos::view_alloc("buf_x", WithoutInitializing),
                       x_stride * 2,
                       ny,
                       nz_tile);
  buffer_t<T***> buf_transpose(
    Kokkos::view_alloc("transpose", WithoutInitializing),
    pencil.extent(0, Pencil::Y),
    pencil.extent(1, Pencil::Y),
    nz_tile);
  const auto    n_threads = omp_get_max_threads();
  buffer_t<T**> complex_y_buffer(
    Kokkos::view_alloc("", WithoutInitializing), y_stride * 2, n_threads);

  for (int z_begin = 0; z_begin < nz; z_begin += nz_tile) {
    const auto z_end = std::min(z_begin + nz_tile, nz);
    const auto buf_y =
      subview(buf_transpose, ALL, ALL, index_range(0, z_end - z_begin));

    // r2c transform in the x-direction
#pragma omp parallel for collapse(2)
    for (int k = 0; k < z_end - z_begin; ++k) {
      for (int j = 0; j < ny; ++j) {
        auto* ptr_input = &(input(0, j, k + z_begin));
        auto* ptr_cmplx = (fftw_complex_t*)&(buf_x(0, j, k));
        interface::exec_r2c(plan_r2c_x.get(), ptr_input, ptr_cmplx);
      }
    }

    // r2c and c2r in the y-direction
    transpose_xy(buf_y, buf_x, pencil, space);

    const auto ny_physical = pencil.global_extent(1);
#pragma omp parallel
    {
      auto  thread_id = omp_get_thread_num();
      auto* ptr_cmplx = &(complex_y_buffer(0, thread_id));
#pragma omp for collapse(2)
      for (int k = 0; k < buf_y.extent_int(2); ++k) {
        for (int j = 0; j < buf_y.extent_int(1); ++j) {
          auto* ptr_input = &(buf_y(0, j, k));
          interface::exec_r2c(
            plan_r2c_y.get(), ptr_input, (fftw_complex_t*)ptr_cmplx);

          omp_pad_and_normalize(ptr_cmplx, ny_physical, kc_y * 2);

          interface::exec_c2r(
            plan_c2r_y.get(), (fftw_complex_t*)ptr_cmplx, ptr_input);
        }
      }
    }

    transpose_yx(buf_x, buf_y, pencil, space);

    // c2r transform in the x-direction
    const auto nx_physical = pencil.global_extent(0);
#pragma omp parallel for collapse(2)
    for (int k = 0; k < z_end - z_begin; ++k) {
      for (int j = 0; j < ny; ++j) {
        auto* ptr_input = &(input(0, j, k + z_begin));
        auto* ptr_cmplx = &(buf_x(0, j, k));

        omp_pad_and_normalize(ptr_cmplx, nx_physical, kc_x * 2);
        interface::exec_c2r(
          plan_c2r_x.get(), (fftw_complex_t*)ptr_cmplx, ptr_input);
      }
    }
  }
}

template<class T>
void SpectralPlan<T, Kokkos::OpenMP, fft::FFTW>::do_fft_r2c_xy(
  view_arg_t const&       output,
  const_view_arg_t const& input,
  ExecSpace const&        space) const
{
  using Kokkos::ALL;
  using Kokkos::subview;
  using index_range = std::pair<int, int>;

  logger_->trace("Calculate 2D R2C transform on {} ({}x{}x{}) and write to {} "
                 "({}x{}x{})",
                 input.label(),
                 input.extent(0),
                 input.extent(1),
                 input.extent(2),
                 output.label(),
                 output.extent(0),
                 output.extent(1),
                 output.extent(2));

  const auto nz = input.extent_int(2);

  // Perform transpose and transform in chunks in the z-direction
  // (improving cache locality)
  const int nz_tile = (pencil.comm().dims[1] > 1) ? nz
                    : (transpose_and_transform_chunk_size_z < 1)
                      ? nz
                      : transpose_and_transform_chunk_size_z;

  const auto     x_stride = nextMultipleComplex<T>(pencil.extent(0) / 2 + 1);
  buffer_t<T***> buf_x(Kokkos::view_alloc("buf_x", WithoutInitializing),
                       x_stride * 2,
                       pencil.extent(1),
                       nz_tile);
  buffer_t<T***> buf_transpose(
    Kokkos::view_alloc("transpose", WithoutInitializing),
    pencil.extent(0, Pencil::Y),
    pencil.extent(1, Pencil::Y),
    nz_tile);
  for (int z_begin = 0; z_begin < nz; z_begin += nz_tile) {
    const auto z_end = std::min(z_begin + nz_tile, nz);
    const auto buf_y =
      subview(buf_transpose, ALL, ALL, index_range(0, z_end - z_begin));

    // r2c transform in the x-direction
#pragma omp parallel for collapse(2)
    for (int k = 0; k < z_end - z_begin; ++k) {
      for (int j = 0; j < buf_x.extent_int(1); ++j) {
        auto* ptr_input = const_cast<T*>(&(input(0, j, k + z_begin)));
        auto* ptr_cmplx = (fftw_complex_t*)&(buf_x(0, j, k));

        interface::exec_r2c(plan_r2c_x.get(), ptr_input, ptr_cmplx);
      }
    }

    transpose_xy(buf_y, buf_x, pencil, space);

    // r2c transform in the y-direction
#pragma omp parallel for collapse(2)
    for (int k = 0; k < z_end - z_begin; ++k) {
      for (int j = 0; j < buf_y.extent_int(1); ++j) {
        auto* ptr_input = &(buf_y(0, j, k));
        auto* ptr_cmplx = (fftw_complex_t*)&(output(0, j, k + z_begin));

        interface::exec_r2c(plan_r2c_y.get(), ptr_input, ptr_cmplx);
      }
    }
  }
}

template<class T>
void SpectralPlan<T, Kokkos::OpenMP, fft::FFTW>::do_fft_c2r_xy(
  view_arg_t const& output,
  view_arg_t const& input,
  bool              dealias,
  ExecSpace const&  space) const
{
  using Kokkos::ALL;
  using Kokkos::subview;
  using index_range = std::pair<int, int>;

  logger_->trace("Calculate 2D C2R transform on {} ({}x{}x{}) and write to {} "
                 "({}x{}x{})",
                 input.label(),
                 input.extent(0),
                 input.extent(1),
                 input.extent(2),
                 output.label(),
                 output.extent(0),
                 output.extent(1),
                 output.extent(2));

  const auto nz = input.extent_int(2);

  // Perform transpose and transform in chunks in the z-direction
  // (improving cache locality)
  const int nz_tile = (pencil.comm().dims[1] > 1) ? nz
                    : (transpose_and_transform_chunk_size_z < 1)
                      ? nz
                      : transpose_and_transform_chunk_size_z;

  const auto x_stride =
    nextMultiple(pencil.extent(0) / 2 + 1, ComplexBufferMultiples<T>);
  buffer_t<T***> buf_x(Kokkos::view_alloc("buf_x", WithoutInitializing),
                       x_stride * 2,
                       pencil.extent(1),
                       nz_tile);
  buffer_t<T***> buf_transpose(
    Kokkos::view_alloc("transpose", WithoutInitializing),
    pencil.extent(0, Pencil::Y),
    pencil.extent(1, Pencil::Y),
    nz_tile);

  for (int z_begin = 0; z_begin < nz; z_begin += nz_tile) {
    const auto z_end = std::min(z_begin + nz_tile, nz);
    const auto buf_y =
      subview(buf_transpose, ALL, ALL, index_range(0, z_end - z_begin));

    // c2r transform in the y-direction
    const auto ny_physical = pencil.global_extent(1);
#pragma omp parallel for collapse(2)
    for (int k = 0; k < z_end - z_begin; ++k) {
      for (int j = 0; j < buf_y.extent_int(1); ++j) {
        auto* ptr_output = &(buf_y(0, j, k));
        auto* ptr_cmplx  = &(input(0, j, k + z_begin));

        if (dealias) {
          const auto kc_y = (ny_physical / 3) * 2; // cutoff k
          omp_pad_and_normalize(ptr_cmplx, ny_physical, kc_y);
        } else {
          omp_normalize(ptr_cmplx, ny_physical);
        }

        interface::exec_c2r(
          plan_c2r_y.get(), (fftw_complex_t*)ptr_cmplx, ptr_output);
      }
    }

    transpose_yx(buf_x, buf_y, pencil, space);

    const auto nx_physical = pencil.global_extent(0);
#pragma omp parallel for collapse(2)
    for (int k = 0; k < z_end - z_begin; ++k) {
      for (int j = 0; j < buf_x.extent_int(1); ++j) {
        auto* ptr_output = &(output(0, j, k + z_begin));
        auto* ptr_cmplx  = &(buf_x(0, j, k));

        if (dealias) {
          const auto kc_x = (nx_physical / 3) * 2; // cutoff k
          omp_pad_and_normalize(ptr_cmplx, nx_physical, kc_x);
        } else {
          omp_normalize(ptr_cmplx, nx_physical);
        }

        interface::exec_c2r(
          plan_c2r_x.get(), (fftw_complex_t*)ptr_cmplx, ptr_output);
      }
    }
  }
}

namespace {
template<class T>
void omp_normalize(T* ptr_cmplx, int n)
{
  static_assert(std::is_floating_point_v<T>);
#pragma omp simd
  for (int i = 0; i < n + 2; ++i) {
    ptr_cmplx[i] /= n;
  }
}

template<class T>
void omp_pad_and_normalize(T* ptr_cmplx, int n, int kc)
{
  static_assert(std::is_floating_point_v<T>);
#pragma omp simd
  for (int i = 0; i < kc; ++i) {
    ptr_cmplx[i] /= n;
  }
#pragma omp simd
  for (int i = kc; i < n + 2; ++i) {
    ptr_cmplx[i] = 0;
  }
}
} // namespace

template class SpectralPlan<float, Kokkos::OpenMP, fft::FFTW>;
template class SpectralPlan<double, Kokkos::OpenMP, fft::FFTW>;

template<>
std::unique_ptr<SpectralPlanBase<float, Kokkos::OpenMP>>
SpectralPlanFactory::create<float, Kokkos::OpenMP, fft::FFTW>(
  PencilPlan const& partition,
  float             kx0,
  float             ky0,
  SpectralOptions const& /*options*/)
{
  return std::make_unique<SpectralPlan<float, Kokkos::OpenMP, fft::FFTW>>(
    partition, kx0, ky0);
}
template<>
std::unique_ptr<SpectralPlanBase<double, Kokkos::OpenMP>>
SpectralPlanFactory::create<double, Kokkos::OpenMP, fft::FFTW>(
  PencilPlan const& partition,
  double            kx0,
  double            ky0,
  SpectralOptions const& /*options*/)
{
  return std::make_unique<SpectralPlan<double, Kokkos::OpenMP, fft::FFTW>>(
    partition, kx0, ky0);
}

} // namespace alps::spectral
