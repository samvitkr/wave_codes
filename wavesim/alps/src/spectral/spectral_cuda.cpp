#include "spectral_cuda.h"

#include <common/base/logging.h>
#include <common/base/macros.h>
#include <common/device/devices.h>
#include <common/kokkos_abstraction/exec_policy.h>
#include <common/kokkos_abstraction/pool_space.h>
#include <decomp/pencil_plan.h>
#include <decomp/pencil_transpose.h>

#include <Kokkos_Core.hpp>

#include <string_view>
#include <utility>

namespace alps::spectral {

template<class DataType, class... Props>
using buffer_t = Kokkos::View<DataType,
                              Kokkos::LayoutLeft,
                              ::alps::memory_pool<Kokkos::CudaSpace>,
                              Props...>;

namespace {
/** Return the smallest multiple of a factor no less than n
 */
template<typename iType1, typename iType2>
constexpr int nextMultiple(iType1 n, iType2 factor)
{
  return static_cast<int>(n + (factor - n % factor) % factor);
}

constexpr std::size_t ComplexBufferAlignment = 32u;
template<typename RealT>
constexpr std::size_t ComplexBufferMultiples =
  ComplexBufferAlignment / sizeof(std::complex<RealT>);

struct MulIkOrder1
{};
struct MulIkOrder2
{};
struct ScaleOnly
{};

template<typename Tag, typename ComplexT>
__global__ void mul_ik_kernel(
  Kokkos::View<ComplexT***, Kokkos::LayoutLeft, Kokkos::CudaSpace> data,
  unsigned                                                         extent2,
  typename ComplexT::value_type                                    k0,
  unsigned                                                         n,
  unsigned                                                         kc)
{
  auto const i       = blockIdx.x * blockDim.x + threadIdx.x;
  auto const j       = blockIdx.y * blockDim.y + threadIdx.y;
  auto const k       = blockIdx.z * blockDim.z + threadIdx.z;
  auto const extent1 = data.extent(1);

  if (j < extent1 && k < extent2) {
    if (i < kc) {
      if constexpr (std::is_same_v<Tag, MulIkOrder1>) {
        auto ik = i * (k0 / n);
        data(i, j, k) *= ComplexT(0, ik);
      } else if constexpr (std::is_same_v<Tag, MulIkOrder2>) {
        data(i, j, k) *= -(i * k0) * (i * k0) / n;
      } else if constexpr (std::is_same_v<Tag, ScaleOnly>) {
        data(i, j, k) /= n;
      } else {
        static_assert(!std::is_same_v<Tag, Tag>, "Unknown tag");
      }
    } else if (i < n / 2 + 1) {
      data(i, j, k) = 0;
    }
  }
  (void)k0; // suppress unused parameter warning
}

template<typename Tag, typename ComplexT>
void execute_mul_ik(
  Kokkos::View<ComplexT***, Kokkos::LayoutLeft, Kokkos::CudaSpace> data,
  int                                                              nz,
  int                                                              n,
  typename ComplexT::value_type                                    k0,
  Kokkos::DefaultExecutionSpace                                    stream)
{
  auto extent0 = static_cast<unsigned>(n / 2 + 1);
  auto extent1 = static_cast<unsigned>(data.extent(1));

  dim3 block(32, 8, 1);
  dim3 grid(
    (extent0 + block.x - 1) / block.x, (extent1 + block.y - 1) / block.y, nz);
  mul_ik_kernel<Tag, ComplexT><<<grid, block, 0, stream.cuda_stream()>>>(
    data, (unsigned)nz, k0, (unsigned)n, (unsigned)n / 2);
  ALPS_CHECK_LAST_DEVICE_ERROR();
}

template<typename ComplexT>
void execute_cutoff(
  Kokkos::View<ComplexT***, Kokkos::LayoutLeft, Kokkos::CudaSpace> data,
  int                                                              nz,
  int                                                              n,
  int                                                              kc,
  Kokkos::DefaultExecutionSpace                                    stream)
{
  auto extent0 = static_cast<unsigned>(n / 2 + 1);
  auto extent1 = static_cast<unsigned>(data.extent(1));

  dim3 block(32, 16, 1);
  dim3 grid(
    (extent0 + block.x - 1) / block.x, (extent1 + block.y - 1) / block.y, nz);
  mul_ik_kernel<ScaleOnly, ComplexT><<<grid, block, 0, stream.cuda_stream()>>>(
    data, (unsigned)nz, 1, (unsigned)n, (unsigned)kc);
  ALPS_CHECK_LAST_DEVICE_ERROR();
}
} // namespace

template<class T>
SpectralPlan<T, Kokkos::Cuda, fft::CUFFT>::SpectralPlan(
  PencilPlan const& partition,
  T                 kx0,
  T                 ky0)
  : base_t{partition, kx0, ky0, fft::CUFFT()}
  , buf_x_spectral(buffer_t<ComplexT***>(
      Kokkos::view_alloc("x dft buffer", Kokkos::WithoutInitializing),
      nextMultiple(pencil.global_extent(0) / 2 + 1, ComplexBufferMultiples<T>),
      pencil.extent(1),
      pencil.extent(2)))
  , buf_y_spectral(buffer_t<ComplexT***>(
      Kokkos::view_alloc("y dft buffer", Kokkos::WithoutInitializing),
      nextMultiple(pencil.global_extent(1) / 2 + 1, ComplexBufferMultiples<T>),
      pencil.extent(1, Pencil::Y),
      pencil.extent(2, Pencil::Y)))
  , buf_y_physical(buffer_t<T***>(
      Kokkos::view_alloc("y reshape buffer", Kokkos::WithoutInitializing),
      create_local_layout(pencil, Pencil::Y)))
{
  setup_plan(pencil.extent(2), Pencil::X, fft::R2C());
  setup_plan(pencil.extent(2), Pencil::Y, fft::R2C());
  setup_plan(pencil.extent(2), Pencil::X, fft::C2R());
  setup_plan(pencil.extent(2), Pencil::Y, fft::C2R());

  logger_->debug("Created a SpectralPlan ({}, {}, {})",
                 std::is_same_v<T, double>  ? "double"
                 : std::is_same_v<T, float> ? "float"
                                            : "unknown type",
                 ExecSpace::name(),
                 "CUFFT");
}

template<class T>
SpectralPlan<T, Kokkos::Cuda, fft::CUFFT>::~SpectralPlan()
{
  logger_->trace("Destroy SpectralPlan ({}, {}, {})",
                 std::is_same_v<T, double>  ? "double"
                 : std::is_same_v<T, float> ? "float"
                                            : "unknown type",
                 ExecSpace::name(),
                 "CUFFT");
}

template<class T>
Kokkos::LayoutLeft
SpectralPlan<T, Kokkos::Cuda, fft::CUFFT>::get_r2c_xy_output_layout() const
{
  const auto& buf = buf_y_spectral;
  return Kokkos::LayoutLeft(buf.extent(0) * 2, buf.extent(1), buf.extent(2));
}

namespace {
template<class T, class Transform>
fft::CUFFTPlan<T, Transform>
create_batch_1d_plan(int n, int howmany, int i_dist, int o_dist)
{
  cufftHandle plan{};

  // Create a cufft handle and set the work area as user managed.
  using cufft_result_id_t = std::underlying_type_t<cufftResult_t>;
  if (auto cufft_result = cufftCreate(&plan); cufft_result != CUFFT_SUCCESS) {
    auto msg = fmt::format("cufftCreate failed with error {}",
                           static_cast<cufft_result_id_t>(cufft_result));
    throw std::runtime_error(msg);
  }
  cufftSetAutoAllocation(plan, 0);

  // Create the cufft plan
  using namespace std::string_view_literals;
  cufftType        type{CUFFT_R2C};
  std::string_view type_name{"R2C"sv};
  if (std::is_same_v<T, double>) {
    type      = std::is_same_v<Transform, fft::R2C> ? CUFFT_D2Z : CUFFT_Z2D;
    type_name = std::is_same_v<Transform, fft::R2C> ? "D2Z"sv : "Z2D"sv;
  } else {
    type      = std::is_same_v<Transform, fft::R2C> ? CUFFT_R2C : CUFFT_C2R;
    type_name = std::is_same_v<Transform, fft::R2C> ? "R2C"sv : "C2R"sv;
  }
  size_t work_size{};

  auto cufft_result = cufftMakePlanMany(plan,
                                        1,
                                        &n,
                                        &i_dist,
                                        1,
                                        i_dist,
                                        &o_dist,
                                        1,
                                        o_dist,
                                        type,
                                        howmany,
                                        &work_size);
  if (cufft_result != CUFFT_SUCCESS) {
    auto msg = fmt::format("cufftMakePlanMany failed with error {}",
                           static_cast<cufft_result_id_t>(cufft_result));
    throw std::runtime_error(msg);
  }

  auto logger = get_logger("fft");
  logger->debug("Created a CUFFT {} plan (id={}) of size {} with batch {}",
                type_name,
                plan,
                n,
                howmany);

  return fft::CUFFTPlan<T, Transform>(plan);
}
} // namespace

template<class T>
void SpectralPlan<T, Kokkos::Cuda, fft::CUFFT>::setup_plan(
  int    nz,
  Pencil direction,
  fft::R2C /*tag*/) const
{
  if (direction == ::alps::Pencil::X) {
    auto n               = pencil.global_extent(0);
    auto batch           = pencil.extent(1, Pencil::X) * nz;
    auto spectral_stride = (int)buf_x_spectral.stride(1);

    plans_x_r2c.try_emplace(
      nz, create_batch_1d_plan<T, fft::R2C>(n, batch, n, spectral_stride));
  } else if (direction == ::alps::Pencil::Y) {
    auto n               = pencil.global_extent(1);
    auto batch           = buf_y_physical.extent_int(1) * nz;
    auto spectral_stride = (int)buf_y_spectral.stride(1);

    plans_y_r2c.try_emplace(
      nz, create_batch_1d_plan<T, fft::R2C>(n, batch, n, spectral_stride));
  }
}

template<class T>
void SpectralPlan<T, Kokkos::Cuda, fft::CUFFT>::setup_plan(
  int    nz,
  Pencil direction,
  fft::C2R /*tag*/) const
{
  if (direction == ::alps::Pencil::X) {
    auto n               = pencil.global_extent(0);
    auto batch           = pencil.extent(1, Pencil::X) * nz;
    auto spectral_stride = (int)buf_x_spectral.stride(1);

    plans_x_c2r.try_emplace(
      nz, create_batch_1d_plan<T, fft::C2R>(n, batch, spectral_stride, n));
  } else if (direction == ::alps::Pencil::Y) {
    auto n               = pencil.global_extent(1);
    auto batch           = buf_y_physical.extent_int(1) * nz;
    auto spectral_stride = (int)buf_y_spectral.stride(1);

    plans_y_c2r.try_emplace(
      nz, create_batch_1d_plan<T, fft::C2R>(n, batch, spectral_stride, n));
  }
}

template<class T>
template<Pencil Direction>
fft::CUFFTPlan<T, fft::R2C> const&
SpectralPlan<T, Kokkos::Cuda, fft::CUFFT>::get_plan(int nz,
                                                    fft::R2C /*tag*/) const
{
  auto& plans = Direction == ::alps::Pencil::X ? plans_x_r2c : plans_y_r2c;
  auto  itr   = plans.find(nz);
  if (itr != plans.cend()) return itr->second;
  setup_plan(nz, Direction, fft::R2C());
  return plans.at(nz);
}

template<class T>
template<Pencil Direction>
fft::CUFFTPlan<T, fft::C2R> const&
SpectralPlan<T, Kokkos::Cuda, fft::CUFFT>::get_plan(int nz,
                                                    fft::C2R /*tag*/) const
{
  auto& plans = Direction == ::alps::Pencil::X ? plans_x_c2r : plans_y_c2r;
  auto  itr   = plans.find(nz);
  if (itr != plans.cend()) return itr->second;
  setup_plan(nz, Direction, fft::C2R());
  return plans.at(nz);
}

// A general function template for do_ddx and do_d2dx2 depending on Tag type
template<class T>
template<class Tag>
void SpectralPlan<T, Kokkos::Cuda, fft::CUFFT>::do_x_derivatives(
  view_arg_t const&       output,
  const_view_arg_t const& input,
  Tag const& /*tag*/,
  Kokkos::Cuda const& space) const
{
  const auto  nz            = input.extent_int(2);
  const auto  nx            = pencil.global_extent(0);
  const auto& forward_plan  = get_plan<::alps::Pencil::X>(nz, fft::R2C());
  const auto& backward_plan = get_plan<::alps::Pencil::X>(nz, fft::C2R());
  const auto& buf           = buf_x_spectral;

  // r2c transform
  forward_plan.run(input.data(), buf.data(), space.cuda_stream());

  // multiply ik
  Kokkos::Profiling::pushRegion("mul_ik");
  execute_mul_ik<Tag, ComplexT>(buf, nz, nx, pex, space);
  Kokkos::Profiling::popRegion();

  // c2r transform
  backward_plan.run(buf.data(), output.data(), space.cuda_stream());
}

template<class T>
void SpectralPlan<T, Kokkos::Cuda, fft::CUFFT>::do_ddx(
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

  do_x_derivatives(output, input, MulIkOrder1(), space);
}

template<class T>
void SpectralPlan<T, Kokkos::Cuda, fft::CUFFT>::do_d2dx2(
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

  do_x_derivatives(output, input, MulIkOrder2(), space);
}

// A general function template for do_ddy and do_d2dy2 depending on Tag type
template<class T>
template<class Tag>
void SpectralPlan<T, Kokkos::Cuda, fft::CUFFT>::do_y_derivatives(
  view_arg_t const&       output,
  const_view_arg_t const& input,
  SpectralPostOp          output_process,
  Tag const& /*tag*/,
  Kokkos::Cuda const& space) const
{
  const auto  nz            = input.extent_int(2);
  const auto  ny            = pencil.global_extent(1);
  const auto& forward_plan  = get_plan<::alps::Pencil::Y>(nz, fft::R2C());
  const auto& backward_plan = get_plan<::alps::Pencil::Y>(nz, fft::C2R());
  const auto& transpose_buf = buf_y_physical;
  const auto& spectral_buf  = buf_y_spectral;

  transpose_xy(transpose_buf, input, pencil, space);

  // r2c transform
  forward_plan.run(
    transpose_buf.data(), spectral_buf.data(), space.cuda_stream());

  // multiply ik
  Kokkos::Profiling::pushRegion("mul_ik");
  execute_mul_ik<Tag, ComplexT>(spectral_buf, nz, ny, pey, space);
  Kokkos::Profiling::popRegion();

  // c2r transform
  backward_plan.run(
    spectral_buf.data(), transpose_buf.data(), space.cuda_stream());

  if (output_process == SpectralPostOp::AssignAfterTranspose) {
    transpose_yx(output, transpose_buf, pencil, space);
  } else if (output_process == SpectralPostOp::AddAfterTranspose) {
    transpose_yx_and_add(output, transpose_buf, pencil, space);
  }
}

template<class T>
void SpectralPlan<T, Kokkos::Cuda, fft::CUFFT>::do_ddy(
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

  do_y_derivatives(output, input, output_process, MulIkOrder1(), space);
}

template<class T>
void SpectralPlan<T, Kokkos::Cuda, fft::CUFFT>::do_d2dy2(
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

  do_y_derivatives(output, input, output_process, MulIkOrder2(), space);
}

template<class T>
void SpectralPlan<T, Kokkos::Cuda, fft::CUFFT>::do_cutoff_xy(
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

  const auto  nx_physical = pencil.global_extent(0);
  const auto  ny_physical = pencil.global_extent(1);
  const auto  nz          = input.extent_int(2);
  const auto& x_forward   = get_plan<::alps::Pencil::X>(nz, fft::R2C());
  const auto& y_forward   = get_plan<::alps::Pencil::Y>(nz, fft::R2C());
  const auto& x_backward  = get_plan<::alps::Pencil::X>(nz, fft::C2R());
  const auto& y_backward  = get_plan<::alps::Pencil::Y>(nz, fft::C2R());

  auto buf_y = subview(buf_y_physical, ALL, ALL, index_range(0, nz));
  auto buf_x_spectral_real_alias = buffer_t<T***>((T*)buf_x_spectral.data(),
                                                  buf_x_spectral.extent(0) * 2,
                                                  buf_x_spectral.extent(1),
                                                  nz);

  // r2c transform in the x-direction
  x_forward.run(input.data(), buf_x_spectral.data(), space.cuda_stream());

  // r2c transform in the y-direction
  transpose_xy(buf_y, buf_x_spectral_real_alias, pencil, space);

  y_forward.run(buf_y.data(), buf_y_spectral.data(), space.cuda_stream());

  // cutoff in the y-direction
  Kokkos::Profiling::pushRegion("cut_y");
  execute_cutoff(buf_y_spectral, nz, ny_physical, kc_y, space);
  Kokkos::Profiling::popRegion();

  // c2r transform in the y-direction
  y_backward.run(buf_y_spectral.data(), buf_y.data(), space.cuda_stream());

  transpose_yx(buf_x_spectral_real_alias, buf_y, pencil, space);

  // c2r transform in the x-direction
  Kokkos::Profiling::pushRegion("cut_x");
  execute_cutoff(buf_x_spectral, nz, nx_physical, kc_x, space);
  Kokkos::Profiling::popRegion();

  x_backward.run(buf_x_spectral.data(), input.data(), space.cuda_stream());
}

template<class T>
void SpectralPlan<T, Kokkos::Cuda, fft::CUFFT>::do_fft_r2c_xy(
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

  auto        nz        = input.extent_int(2);
  const auto& x_forward = get_plan<::alps::Pencil::X>(nz, fft::R2C());
  const auto& y_forward = get_plan<::alps::Pencil::Y>(nz, fft::R2C());

  auto buf_y = subview(buf_y_physical, ALL, ALL, index_range(0, nz));
  auto buf_x_spectral_real_alias = buffer_t<T***>((T*)buf_x_spectral.data(),
                                                  buf_x_spectral.extent(0) * 2,
                                                  buf_x_spectral.extent(1),
                                                  nz);

  // r2c transform in the x-direction
  x_forward.run(input.data(), buf_x_spectral.data(), space.cuda_stream());

  // r2c transform in the y-direction
  transpose_xy(buf_y, buf_x_spectral_real_alias, pencil, space);

  y_forward.run(buf_y.data(), output.data(), space.cuda_stream());
}

template<class T>
void SpectralPlan<T, Kokkos::Cuda, fft::CUFFT>::do_fft_c2r_xy(
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

  const auto  nx_physical = pencil.global_extent(0);
  const auto  ny_physical = pencil.global_extent(1);
  const auto  nx_half     = nx_physical / 2;
  const auto  ny_half     = ny_physical / 2;
  auto        nz          = input.extent_int(2);
  const auto& x_plan      = get_plan<::alps::Pencil::X>(nz, fft::C2R());
  const auto& y_plan      = get_plan<::alps::Pencil::Y>(nz, fft::C2R());

  auto buf_y = subview(buf_y_physical, ALL, ALL, index_range(0, nz));
  auto buf_x_spectral_real_alias = buffer_t<T***>((T*)buf_x_spectral.data(),
                                                  buf_x_spectral.extent(0) * 2,
                                                  buf_x_spectral.extent(1),
                                                  nz);
  Kokkos::View<ComplexT***,
               typename view_arg_t::array_layout,
               typename view_arg_t::memory_space> const
    input_cmplx(
      (ComplexT*)input.data(), input.extent(0) / 2, input.extent(1), nz);

  // cutoff in the y-direction
  const auto kc_y = dealias ? ny_half * 2 / 3 : ny_half;
  Kokkos::Profiling::pushRegion("scale_y");
  execute_cutoff(input_cmplx, nz, ny_physical, kc_y, space);
  Kokkos::Profiling::popRegion();

  y_plan.run(input_cmplx.data(), buf_y.data(), space.cuda_stream());

  transpose_yx(buf_x_spectral_real_alias, buf_y, pencil, space);

  const auto kc_x = dealias ? nx_half * 2 / 3 : nx_half;
  Kokkos::Profiling::pushRegion("scale_x");
  execute_cutoff(buf_x_spectral, nz, nx_physical, kc_x, space);
  Kokkos::Profiling::popRegion();

  x_plan.run(buf_x_spectral.data(), output.data(), space.cuda_stream());
}

template class SpectralPlan<float, Kokkos::Cuda, fft::CUFFT>;
template class SpectralPlan<double, Kokkos::Cuda, fft::CUFFT>;

template<>
std::unique_ptr<SpectralPlanBase<float, Kokkos::Cuda>>
SpectralPlanFactory::create<float, Kokkos::Cuda, fft::CUFFT>(
  PencilPlan const& partition,
  float             kx0,
  float             ky0,
  SpectralOptions const& /*options*/)
{
  return std::make_unique<SpectralPlan<float, Kokkos::Cuda, fft::CUFFT>>(
    partition, kx0, ky0);
}
template<>
std::unique_ptr<SpectralPlanBase<double, Kokkos::Cuda>>
SpectralPlanFactory::create<double, Kokkos::Cuda, fft::CUFFT>(
  PencilPlan const& partition,
  double            kx0,
  double            ky0,
  SpectralOptions const& /*options*/)
{
  return std::make_unique<SpectralPlan<double, Kokkos::Cuda, fft::CUFFT>>(
    partition, kx0, ky0);
}

} // namespace alps::spectral
