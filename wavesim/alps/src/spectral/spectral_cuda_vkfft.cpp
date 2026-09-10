#include "spectral_cuda_vkfft.h"

#include <common/base/logging.h>
#include <common/kokkos_abstraction/pool_space.h>
#include <decomp/pencil_plan.h>
#include <decomp/pencil_transpose.h>
#include <fft/vkfft.h>
#include <spectral/detail/vkfft_app_utils.h>

#include <Kokkos_Core.hpp>
#if defined(KOKKOS_ENABLE_CUDA)
#include <cuda.h>
#elif defined(KOKKOS_ENABLE_HIP)
#include <hip/hip_runtime.h>
#endif

#include <stdexcept>
#include <utility>

namespace alps::spectral {

template<class DataType, class... Props>
using buffer_t =
  Kokkos::View<DataType,
               Kokkos::LayoutLeft,
               ::alps::memory_pool<Kokkos::DefaultExecutionSpace::memory_space>,
               Props...>;

namespace {
/** Return the smallest multiple of a factor no less than n
 */
template<typename iType1, typename iType2>
constexpr int nextMultiple(iType1 n, iType2 factor)
{
  return static_cast<int>(n + (factor - n % factor) % factor);
}

bool check_radix(int n, std::vector<int> radixes)
{
  if (n < 0) return false;
  for (const auto radix : radixes) {
    if (radix > 1 && radix <= n) {
      while (n % radix == 0) {
        n /= radix;
      }
    }
  }
  return n == 1;
}

constexpr std::size_t ComplexBufferAlignment = 32u;
template<typename RealT>
constexpr std::size_t ComplexBufferMultiples =
  ComplexBufferAlignment / sizeof(std::complex<RealT>);

template<typename RealT, typename iType>
constexpr auto nextMultipleComplex(iType n)
{
  return nextMultiple(n, ComplexBufferMultiples<RealT>);
}

template<class T>
Kokkos::View<Kokkos::complex<T>* [2],
             Kokkos::LayoutLeft,
             Kokkos::DefaultExecutionSpace::memory_space>
prepare_derivative_kernel(T k0, int n);

template<class T>
__global__ void
prepare_cutoff_kernel(Kokkos::complex<T>* kernel, int n, int cutoff_idx);

bool check_vkfft_plan(VkFFTApplication const* app,
                      VkFFTPlan const*        plan,
                      Logger                  logger);

template<typename Transform, typename stream_t>
void exec(VkFFTApplication* app,
          VkFFTLaunchParams launchParams,
          Transform const& /*tag*/,
          stream_t stream,
          Logger   logger)
{
  // Overwrite the execution stream and save the original stream in the app
  // configuration for cleanup
  auto old_num_streams = std::exchange(app->configuration.num_streams, 1);
  auto old_stream      = std::exchange(app->configuration.stream, &stream);

  int const inverse = std::is_same_v<Transform, fft::C2R> ? 1 : -1;
  if (auto result = VkFFTAppend(app, inverse, &launchParams);
      result != VKFFT_SUCCESS) {
    logger->error(
      "VKFFT execution error code {}",
      static_cast<std::underlying_type_t<decltype(result)>>(result));
  }

  // Restore the original stream configuration in the app for proper cleanup
  app->configuration.num_streams = old_num_streams;
  app->configuration.stream      = old_stream;
}
} // namespace

template<class T>
SpectralPlan<T, Kokkos::DefaultExecutionSpace, fft::VKFFT>::SpectralPlan(
  PencilPlan const&      partition,
  T                      kx0,
  T                      ky0,
  SpectralOptions const& options)
  : base_t{partition, kx0, ky0, fft::VKFFT()}
  , buf_spectral(buffer_t<ComplexT**>(
      Kokkos::view_alloc("dft buffer", Kokkos::WithoutInitializing),
      std::max(nextMultipleComplex<T>(pencil.extent(0) / 2 + 1)
                 * pencil.extent(1),
               nextMultipleComplex<T>(pencil.extent(0, Pencil::Y) / 2 + 1)
                 * pencil.extent(1, Pencil::Y)),
      pencil.extent(2))) // the larger size of the x and y transform results
  , buf_y_physical(buffer_t<T***>(
      Kokkos::view_alloc("y reshape buffer", Kokkos::WithoutInitializing),
      create_local_layout(pencil, Pencil::Y)))
  , kernel_cutoff_x(buffer_t<ComplexT*>(
      Kokkos::view_alloc("cutoff x kernel", Kokkos::WithoutInitializing),
      pencil.global_extent(0)))
  , kernel_cutoff_y(buffer_t<ComplexT*>(
      Kokkos::view_alloc("cutoff y kernel", Kokkos::WithoutInitializing),
      pencil.global_extent(1)))
  , tuning_config_(detail::make_vkfft_tuning_config(options.vkfft_tuning))
{
  if (!check_radix(pencil.global_extent(0), {2, 3, 5, 7, 11, 13})
      || !check_radix(pencil.global_extent(1), {2, 3, 5, 7, 11, 13})) {
    if (pencil.comm().rank() == 0) {
      throw std::invalid_argument(
        fmt::format("The Nx ({}) and Ny ({}) should be decomposed to radixes "
                    "of (2, 3, 5, 7, 11, 13) for vkFFT.",
                    pencil.global_extent(0),
                    pencil.global_extent(1)));
    }
  }

  kernel_ddx = prepare_derivative_kernel(kx0, pencil.global_extent(0));
  kernel_ddy = prepare_derivative_kernel(ky0, pencil.global_extent(1));
  Kokkos::deep_copy(kernel_cutoff_x, Kokkos::complex<T>(0, 0));
  Kokkos::deep_copy(kernel_cutoff_y, Kokkos::complex<T>(0, 0));

  logger_->debug("Created a SpectralPlan ({}, {}, {})",
                 std::is_same_v<T, double> ? "double" : "float",
                 ExecSpace::name(),
                 "VKFFT");
}

template<class T>
SpectralPlan<T, Kokkos::DefaultExecutionSpace, fft::VKFFT>::~SpectralPlan()
{
  logger_->trace("Destroy SpectralPlan ({}, {}, {})",
                 std::is_same_v<T, double> ? "double" : "float",
                 ExecSpace::name(),
                 "VKFFT");
}

template<class T>
Kokkos::LayoutLeft SpectralPlan<T, Kokkos::DefaultExecutionSpace, fft::VKFFT>::
  get_r2c_xy_output_layout() const
{
  return Kokkos::LayoutLeft(
    nextMultipleComplex<T>(pencil.extent(0, Pencil::Y) / 2 + 1) * 2,
    pencil.extent(1, Pencil::Y),
    pencil.extent(2, Pencil::Y));
}

template<class T>
template<::alps::Pencil Direction>
VkFFTConfiguration SpectralPlan<T, Kokkos::DefaultExecutionSpace, fft::VKFFT>::
  get_plan_base_config() const
{
  static_assert(Direction == ::alps::Pencil::X
                || Direction == ::alps::Pencil::Y);

  VkFFTConfiguration config{};
  config.FFTdim           = 1;
  config.performR2C       = 1;
  config.doublePrecision  = std::is_same_v<T, double> ? 1 : 0;
  config.normalize        = 1;
  config.isInputFormatted = 1; // Read input from another array and write
  // output to another array when performing R2C

  if constexpr (Direction == ::alps::Pencil::X) {
    // transform configuration specific to x direction
    config.size[0]              = pencil.global_extent(0);
    config.inputBufferStride[0] = pencil.global_extent(0);
    config.bufferStride[0] = nextMultipleComplex<T>(config.size[0] / 2 + 1);
  } else {
    // transform configuration specific to y direction
    config.size[0]              = pencil.global_extent(1);
    config.inputBufferStride[0] = pencil.global_extent(1);
    config.bufferStride[0] = nextMultipleComplex<T>(config.size[0] / 2 + 1);
  }
  return config;
}

template<class T>
template<::alps::Pencil Direction>
const auto&
SpectralPlan<T, Kokkos::DefaultExecutionSpace, fft::VKFFT>::get_plan(
  int nz) const
{
  static_assert(Direction == ::alps::Pencil::X
                || Direction == ::alps::Pencil::Y);

  auto itr = plans.find(
    {pencil.extent(0, Direction), nz * pencil.extent(1, Direction), 0});
  if (itr == plans.cend()) {
    itr = emplace_plan<Direction>(nz);
  }
  return itr->second.plan;
}

template<class T>
template<::alps::Pencil Direction>
auto SpectralPlan<T, Kokkos::DefaultExecutionSpace, fft::VKFFT>::emplace_plan(
  int nz) const
{
  static_assert(Direction == ::alps::Pencil::X
                || Direction == ::alps::Pencil::Y);

  VKFFTPlanWithStorage new_plan{VkApp{nullptr, &deleteVkFFT}};
#if defined(KOKKOS_ENABLE_CUDA)
  cuDeviceGet(&(new_plan.device), Kokkos::Cuda().cuda_device());
#elif defined(KOKKOS_ENABLE_HIP)
  KOKKOS_IMPL_HIP_SAFE_CALL(
    hipDeviceGet(&(new_plan.device), Kokkos::HIP().hip_device()));
#endif

  VkFFTConfiguration config;
  if constexpr (Direction == ::alps::Pencil::X) {
    // transform in the x direction
    config               = get_plan_base_config<::alps::Pencil::X>();
    config.numberBatches = (uint64_t)nz * pencil.extent(1, ::alps::Pencil::X);
    new_plan.buffer_size = buf_spectral.span() * sizeof(ComplexT);
    new_plan.buffer      = buf_spectral.data();
  } else {
    // transform in the y direction
    config               = get_plan_base_config<::alps::Pencil::Y>();
    config.numberBatches = (uint64_t)nz * pencil.extent(1, ::alps::Pencil::Y);
    new_plan.buffer_size = buf_spectral.span() * sizeof(ComplexT);
    new_plan.buffer      = buf_spectral.data();
  }
  config.device                     = &(new_plan.device);
  config.bufferSize                 = &(new_plan.buffer_size);
  config.buffer                     = &(new_plan.buffer);
  config.inverseReturnToInputBuffer = 1; // Return output to input buffer when
                                         // performing C2R
  if (tuning_config_.enabled) {
    auto tune_input = buffer_t<T***>(
      Kokkos::view_alloc("vkfft tune input", Kokkos::WithoutInitializing),
      pencil.extent(0, Direction),
      pencil.extent(1, Direction),
      nz);
    detail::VkFFTExecStream tune_stream = nullptr;

    new_plan.plan = detail::tune_scheduler_params(
      [&]() -> VkFFTConfiguration { return config; },
      [&](VkFFTApplication* app, detail::VkFFTExecStream stream) {
        VkFFTLaunchParams params{};
        void*             input_ptr = tune_input.data();
        params.inputBuffer          = &input_ptr;

        exec(app, params, fft::R2C{}, stream, logger_);
        exec(app, params, fft::C2R{}, stream, logger_);
      },
      tune_stream,
      tuning_config_,
      logger_);
  } else {
    new_plan.plan = detail::make_vkfft_app(config);
  }
  if (!check_vkfft_plan(
        new_plan.plan.get(), new_plan.plan->localFFTPlan_inverse, logger_)) {
    throw std::runtime_error("Transform size too large for VkFFT");
  }

  logger_->debug("Created a VkFFT plan for size {}x{}.",
                 config.size[0],
                 config.numberBatches);
  return plans
    .try_emplace(
      {pencil.extent(0, Direction), nz * pencil.extent(1, Direction), 0},
      std::move(new_plan))
    .first;
}

template<class T>
template<::alps::Pencil Direction>
const auto&
SpectralPlan<T, Kokkos::DefaultExecutionSpace, fft::VKFFT>::get_cutoff_plan(
  int nz,
  int cutoff_idx) const
{
  static_assert(Direction == ::alps::Pencil::X
                || Direction == ::alps::Pencil::Y);

  auto itr = cutoff_plans.find({pencil.extent(0, Direction),
                                nz * pencil.extent(1, Direction),
                                cutoff_idx});
  if (itr == cutoff_plans.cend()) {
    itr = emplace_cutoff_plan<Direction>(nz, cutoff_idx);
  }
  return itr->second.plan;
}

template<class T>
template<::alps::Pencil Direction>
auto SpectralPlan<T, Kokkos::DefaultExecutionSpace, fft::VKFFT>::
  emplace_cutoff_plan(int nz, int cutoff_idx) const
{
  static_assert(Direction == ::alps::Pencil::X
                || Direction == ::alps::Pencil::Y);

  VKFFTPlanWithStorage new_plan{VkApp{nullptr, &deleteVkFFT}};
#if defined(KOKKOS_ENABLE_CUDA)
  cuDeviceGet(&(new_plan.device), Kokkos::Cuda().cuda_device());
#elif defined(KOKKOS_ENABLE_HIP)
  KOKKOS_IMPL_HIP_SAFE_CALL(
    hipDeviceGet(&(new_plan.device), Kokkos::HIP().hip_device()));
#endif

  VkFFTConfiguration config;
  if constexpr (Direction == ::alps::Pencil::X) {
    // transform in the x direction
    config               = get_plan_base_config<::alps::Pencil::X>();
    config.numberBatches = (uint64_t)nz * pencil.extent(1, ::alps::Pencil::X);
    new_plan.buffer_size = buf_spectral.span() * sizeof(ComplexT);
    new_plan.buffer      = buf_spectral.data();
  } else {
    // transform in the y direction
    config               = get_plan_base_config<::alps::Pencil::Y>();
    config.numberBatches = (uint64_t)nz * pencil.extent(1, ::alps::Pencil::Y);
    new_plan.buffer_size = buf_spectral.span() * sizeof(ComplexT);
    new_plan.buffer      = buf_spectral.data();
  }
  config.device                     = &(new_plan.device);
  config.bufferSize                 = &(new_plan.buffer_size);
  config.buffer                     = &(new_plan.buffer);
  config.performZeropadding[0]      = 1;
  config.frequencyZeroPadding       = 1;
  config.fft_zeropad_left[0]        = cutoff_idx;
  config.fft_zeropad_right[0]       = config.size[0] / 2 + 1;
  config.inverseReturnToInputBuffer = 1; // Return output to input buffer when
                                         // performing C2R
  config.makeInversePlanOnly = 1;

  if (tuning_config_.enabled) {
    auto tune_input =
      buffer_t<T***>(Kokkos::view_alloc("vkfft tune cutoff input",
                                        Kokkos::WithoutInitializing),
                     pencil.extent(0, Direction),
                     pencil.extent(1, Direction),
                     nz);
    detail::VkFFTExecStream tune_stream = nullptr;

    new_plan.plan = detail::tune_scheduler_params(
      [&]() -> VkFFTConfiguration { return config; },
      [&](VkFFTApplication* app, detail::VkFFTExecStream stream) {
        VkFFTLaunchParams params{};
        void*             input_ptr  = tune_input.data();
        void*             buffer_ptr = buf_spectral.data();
        params.inputBuffer           = &input_ptr;
        params.buffer                = &buffer_ptr;

        exec(app, params, fft::C2R{}, stream, logger_);
      },
      tune_stream,
      tuning_config_,
      logger_);
  } else {
    new_plan.plan = detail::make_vkfft_app(config);
  }
  if (!check_vkfft_plan(
        new_plan.plan.get(), new_plan.plan->localFFTPlan_inverse, logger_)) {
    throw std::runtime_error("Transform size too large for VkFFT");
  }

  logger_->debug("Created a VkFFT cutoff plan for size {}x{} (cutoff={}).",
                 config.size[0],
                 config.numberBatches,
                 cutoff_idx);
  return cutoff_plans
    .try_emplace({pencil.extent(0, Direction),
                  nz * pencil.extent(1, Direction),
                  cutoff_idx},
                 std::move(new_plan))
    .first;
}

template<class T>
template<::alps::Pencil Direction>
const auto& SpectralPlan<T, Kokkos::DefaultExecutionSpace, fft::VKFFT>::
  get_convolution_plan(int nz) const
{
  static_assert(Direction == ::alps::Pencil::X
                || Direction == ::alps::Pencil::Y);

  auto itr = conv_plans.find(
    {pencil.extent(0, Direction), nz * pencil.extent(1, Direction), 0});
  if (itr == conv_plans.cend()) {
    itr = emplace_convolution_plan<Direction>(nz);
  }
  return itr->second.plan;
}

template<class T>
template<::alps::Pencil Direction>
auto SpectralPlan<T, Kokkos::DefaultExecutionSpace, fft::VKFFT>::
  emplace_convolution_plan(int nz) const
{
  static_assert(Direction == ::alps::Pencil::X
                || Direction == ::alps::Pencil::Y);

  VKFFTPlanWithStorage new_plan{VkApp{nullptr, &deleteVkFFT}};
#if defined(KOKKOS_ENABLE_CUDA)
  cuDeviceGet(&(new_plan.device), Kokkos::Cuda().cuda_device());
#elif defined(KOKKOS_ENABLE_HIP)
  KOKKOS_IMPL_HIP_SAFE_CALL(
    hipDeviceGet(&(new_plan.device), Kokkos::HIP().hip_device()));
#endif

  VkFFTConfiguration config;
  if constexpr (Direction == ::alps::Pencil::X) {
    // transform in the x direction
    config               = get_plan_base_config<::alps::Pencil::X>();
    config.numberBatches = (uint64_t)nz * pencil.extent(1, ::alps::Pencil::X);
    new_plan.buffer_size = buf_spectral.span() * sizeof(ComplexT);
    new_plan.kernel_size = kernel_ddx.stride(1) * sizeof(ComplexT);
    new_plan.buffer      = buf_spectral.data();
  } else {
    // transform in the y direction
    config               = get_plan_base_config<::alps::Pencil::Y>();
    config.numberBatches = (uint64_t)nz * pencil.extent(1, ::alps::Pencil::Y);
    new_plan.buffer_size = buf_spectral.span() * sizeof(ComplexT);
    new_plan.kernel_size = kernel_ddy.stride(1) * sizeof(ComplexT);
    new_plan.buffer      = buf_spectral.data();
  }
  config.device                      = &(new_plan.device);
  config.isOutputFormatted           = true;
  config.outputBufferStride[0]       = config.inputBufferStride[0];
  config.performConvolution          = 1;
  config.numberKernels               = 1;
  config.singleKernelMultipleBatches = 1;
  config.bufferSize                  = &(new_plan.buffer_size);
  config.buffer                      = &(new_plan.buffer);
  config.kernelSize                  = &(new_plan.kernel_size);
  config.kernel                      = nullptr;

  if (tuning_config_.enabled) {
    auto tune_input = buffer_t<T***>(
      Kokkos::view_alloc("vkfft tune conv input", Kokkos::WithoutInitializing),
      pencil.extent(0, Direction),
      pencil.extent(1, Direction),
      nz);
    auto tune_kernel = buffer_t<ComplexT*>(
      Kokkos::view_alloc("vkfft tune kernel", Kokkos::WithoutInitializing),
      new_plan.kernel_size / sizeof(ComplexT));
    Kokkos::deep_copy(tune_kernel, ComplexT(1, 0));
    detail::VkFFTExecStream tune_stream = nullptr;

    new_plan.plan = detail::tune_scheduler_params(
      [&]() -> VkFFTConfiguration { return config; },
      [&](VkFFTApplication* app, detail::VkFFTExecStream stream) {
        VkFFTLaunchParams params{};
        void*             input_ptr  = tune_input.data();
        void*             output_ptr = tune_input.data();
        void*             kernel_ptr = tune_kernel.data();
        params.inputBuffer           = &input_ptr;
        params.outputBuffer          = &output_ptr;
        params.kernel                = &kernel_ptr;

        exec(app, params, fft::R2C{}, stream, logger_);
      },
      tune_stream,
      tuning_config_,
      logger_);
  } else {
    new_plan.plan = detail::make_vkfft_app(config);
  }
  if (!check_vkfft_plan(
        new_plan.plan.get(), new_plan.plan->localFFTPlan_inverse, logger_)) {
    throw std::runtime_error("Transform size too large for VkFFT");
  }

  logger_->debug("Created a VkFFT convolution plan for size {}x{}.",
                 config.size[0],
                 config.numberBatches);
  return conv_plans
    .try_emplace(
      {pencil.extent(0, Direction), nz * pencil.extent(1, Direction), 0},
      std::move(new_plan))
    .first;
}

template<class T>
void SpectralPlan<T, Kokkos::DefaultExecutionSpace, fft::VKFFT>::
  apply_x_convolution(view_arg_t const&       output,
                      const_view_arg_t const& input,
                      ComplexT*               kernel,
                      ExecSpace const&        space) const
{
  int         nz             = input.extent_int(2);
  const auto& transform_plan = get_convolution_plan<::alps::Pencil::X>(nz);

  VkFFTLaunchParams params{};
  auto*             ptr_input  = input.data();
  auto*             ptr_output = output.data();
  params.inputBuffer           = (void**)&ptr_input;
  params.outputBuffer          = (void**)&ptr_output;
  void* ptr_kernel             = kernel;
  params.kernel                = &ptr_kernel;
#if defined(KOKKOS_ENABLE_CUDA)
  exec(transform_plan.get(), params, fft::R2C(), space.cuda_stream(), logger_);
#elif defined(KOKKOS_ENABLE_HIP)
  exec(transform_plan.get(), params, fft::R2C(), space.hip_stream(), logger_);
#endif
}

template<class T>
void SpectralPlan<T, Kokkos::DefaultExecutionSpace, fft::VKFFT>::do_ddx(
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

  apply_x_convolution(output, input, kernel_ddx.data(), space);
}

template<class T>
void SpectralPlan<T, Kokkos::DefaultExecutionSpace, fft::VKFFT>::do_d2dx2(
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

  apply_x_convolution(
    output, input, kernel_ddx.data() + kernel_ddx.stride(1), space);
}

template<class T>
void SpectralPlan<T, Kokkos::DefaultExecutionSpace, fft::VKFFT>::
  apply_y_convolution(view_arg_t const&       output,
                      const_view_arg_t const& input,
                      ComplexT*               kernel,
                      SpectralPostOp          output_process,
                      ExecSpace const&        space) const
{
  int         nz             = input.extent_int(2);
  const auto& transform_plan = get_convolution_plan<::alps::Pencil::Y>(nz);

  transpose_xy(buf_y_physical, input, pencil, space);

  VkFFTLaunchParams params{};
  auto*             ptr_input = buf_y_physical.data();
  params.inputBuffer          = (void**)&ptr_input;
  params.outputBuffer         = (void**)&ptr_input;
  void* ptr_kernel            = kernel;
  params.kernel               = &ptr_kernel;
#if defined(KOKKOS_ENABLE_CUDA)
  exec(transform_plan.get(), params, fft::R2C(), space.cuda_stream(), logger_);
#elif defined(KOKKOS_ENABLE_HIP)
  exec(transform_plan.get(), params, fft::R2C(), space.hip_stream(), logger_);
#endif

  if (output_process == SpectralPostOp::AssignAfterTranspose) {
    transpose_yx(output, buf_y_physical, pencil, space);
  } else if (output_process == SpectralPostOp::AddAfterTranspose) {
    transpose_yx_and_add(output, buf_y_physical, pencil, space);
  }
}

template<class T>
void SpectralPlan<T, Kokkos::DefaultExecutionSpace, fft::VKFFT>::do_ddy(
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

  apply_y_convolution(output, input, kernel_ddy.data(), output_process, space);
}

template<class T>
void SpectralPlan<T, Kokkos::DefaultExecutionSpace, fft::VKFFT>::do_d2dy2(
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

  auto* kernel = kernel_ddy.data() + kernel_ddy.stride(1);
  apply_y_convolution(output, input, kernel, output_process, space);
}

template<class T>
void SpectralPlan<T, Kokkos::DefaultExecutionSpace, fft::VKFFT>::do_cutoff_xy(
  view_arg_t const& input,
  int               kc_x,
  int               kc_y,
  ExecSpace const&  space) const
{
  logger_->trace("Cutoff {} ({}x{}x{}) to wavenumber kx={} ky={}",
                 input.label(),
                 input.extent(0),
                 input.extent(1),
                 input.extent(2),
                 kc_x,
                 kc_y);

  constexpr int BLOCK_SIZE = 256;

  int nx = pencil.global_extent(0);
  int ny = pencil.global_extent(1);
  int nz = input.extent_int(2);

  // temporary storage for low-pass filtered result in the x-direction
  auto buf_x_transformed = buffer_t<T***>(
    (T*)buf_spectral.data(), pencil.extent(0), pencil.extent(1), nz);

  // cutoff in the x-direction as convolution
  // prepare the cutoff kernel
  auto grid_size_x = (nx + BLOCK_SIZE - 1) / BLOCK_SIZE;
#if defined(KOKKOS_ENABLE_CUDA)
  prepare_cutoff_kernel<T><<<grid_size_x, BLOCK_SIZE, 0, space.cuda_stream()>>>(
    kernel_cutoff_x.data(), nx, kc_x);
#elif defined(KOKKOS_ENABLE_HIP)
  prepare_cutoff_kernel<T><<<grid_size_x, BLOCK_SIZE, 0, space.hip_stream()>>>(
    kernel_cutoff_x.data(), nx, kc_x);
#endif

  apply_x_convolution(buf_x_transformed, input, kernel_cutoff_x.data(), space);

  // cutoff in the y-direction as convolution
  // prepare the cutoff kernel
  auto grid_size_y = (ny + BLOCK_SIZE - 1) / BLOCK_SIZE;
#if defined(KOKKOS_ENABLE_CUDA)
  prepare_cutoff_kernel<T><<<grid_size_y, BLOCK_SIZE, 0, space.cuda_stream()>>>(
    kernel_cutoff_y.data(), ny, kc_y);
#elif defined(KOKKOS_ENABLE_HIP)
  prepare_cutoff_kernel<T><<<grid_size_y, BLOCK_SIZE, 0, space.hip_stream()>>>(
    kernel_cutoff_y.data(), ny, kc_y);
#endif

  apply_y_convolution(input,
                      buf_x_transformed,
                      kernel_cutoff_y.data(),
                      SpectralPostOp::AssignAfterTranspose,
                      space);
}

template<class T>
void SpectralPlan<T, Kokkos::DefaultExecutionSpace, fft::VKFFT>::do_fft_r2c_xy(
  view_arg_t const&       output,
  const_view_arg_t const& input,
  ExecSpace const&        space) const
{
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

  int nz = input.extent_int(2);

  auto buf_x_transformed =
    buffer_t<T***>((T*)buf_spectral.data(),
                   nextMultipleComplex<T>(pencil.extent(0) / 2 + 1) * 2,
                   pencil.extent(1),
                   nz);

  // r2c transform in the x-direction
  {
    const auto&       x_plan = get_plan<::alps::Pencil::X>(nz);
    VkFFTLaunchParams params{};
    void* ptr_input = const_cast<void*>(static_cast<const void*>(input.data()));
    void* buffer    = buf_x_transformed.data();
    params.inputBuffer = &ptr_input;
    params.buffer      = &buffer;
#if defined(KOKKOS_ENABLE_CUDA)
    exec(x_plan.get(), params, fft::R2C(), space.cuda_stream(), logger_);
#elif defined(KOKKOS_ENABLE_HIP)
    exec(x_plan.get(), params, fft::R2C(), space.hip_stream(), logger_);
#endif
  }

  // r2c transform in the y-direction
  transpose_xy(buf_y_physical, buf_x_transformed, pencil, space);

  {
    const auto&       y_plan = get_plan<::alps::Pencil::Y>(nz);
    VkFFTLaunchParams params{};
    void*             ptr_output = output.data();
    void*             ptr_input  = buf_y_physical.data();
    params.inputBuffer           = &ptr_input;
    params.buffer                = &ptr_output;
#if defined(KOKKOS_ENABLE_CUDA)
    exec(y_plan.get(), params, fft::R2C(), space.cuda_stream(), logger_);
#elif defined(KOKKOS_ENABLE_HIP)
    exec(y_plan.get(), params, fft::R2C(), space.hip_stream(), logger_);
#endif
  }
}

template<class T>
void SpectralPlan<T, Kokkos::DefaultExecutionSpace, fft::VKFFT>::do_fft_c2r_xy(
  view_arg_t const& output,
  view_arg_t const& input,
  bool              dealias,
  ExecSpace const&  space) const
{
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

  int  nz = input.extent_int(2);
  auto buf_x_transformed =
    buffer_t<T***>((T*)buf_spectral.data(),
                   nextMultipleComplex<T>(pencil.extent(0) / 2 + 1) * 2,
                   pencil.extent(1),
                   nz);

  // c2r transform in the y-direction
  {
    const auto& y_plan =
      (dealias)
        ? get_cutoff_plan<::alps::Pencil::Y>(nz, pencil.global_extent(1) / 3)
        : get_plan<::alps::Pencil::Y>(nz);
    VkFFTLaunchParams params{};
    void*             ptr_input = input.data();
    void*             buffer    = buf_y_physical.data();
    params.inputBuffer          = &buffer;
    params.buffer               = &ptr_input;

#if defined(KOKKOS_ENABLE_CUDA)
    exec(y_plan.get(), params, fft::C2R(), space.cuda_stream(), logger_);
#elif defined(KOKKOS_ENABLE_HIP)
    exec(y_plan.get(), params, fft::C2R(), space.hip_stream(), logger_);
#endif
  }

  transpose_yx(buf_x_transformed, buf_y_physical, pencil, space);

  // c2r transform in the x-direction
  {
    const auto& x_plan =
      (dealias)
        ? get_cutoff_plan<::alps::Pencil::X>(nz, pencil.global_extent(0) / 3)
        : get_plan<::alps::Pencil::X>(nz);
    VkFFTLaunchParams params{};
    void*             ptr_output = output.data();
    void*             buffer     = buf_x_transformed.data();
    params.inputBuffer           = &ptr_output;
    params.buffer                = &buffer;

#if defined(KOKKOS_ENABLE_CUDA)
    exec(x_plan.get(), params, fft::C2R(), space.cuda_stream(), logger_);
#elif defined(KOKKOS_ENABLE_HIP)
    exec(x_plan.get(), params, fft::C2R(), space.hip_stream(), logger_);
#endif
  }
}

namespace {
template<class T>
Kokkos::View<Kokkos::complex<T>* [2],
             Kokkos::LayoutLeft,
             Kokkos::DefaultExecutionSpace::memory_space>
prepare_derivative_kernel(T k0, int n)
{
  using ComplexT = Kokkos::complex<T>;
  buffer_t<ComplexT* [2]> kernel("spectral kernel", n);

  auto kernel_host =
    Kokkos::create_mirror_view(Kokkos::WithoutInitializing, kernel);

  // Set the Fourier coefficients of the kernel
  for (int i = 0; i < n / 2; ++i) {
    kernel_host(i, 0) = ComplexT(0, i * k0);
    for (int j = 1; j < kernel.extent_int(1); ++j) {
      kernel_host(i, j) = kernel_host(i, j - 1) * ComplexT(0, i * k0);
    }
  }
  for (int j = 0; j < kernel.extent_int(1); ++j) {
    kernel_host(n / 2, j) = 0;
  }
  for (int i = n / 2 + 1; i < n; ++i) {
    kernel_host(i, 0) = ComplexT(0, (i - n) * k0);
    for (int j = 1; j < kernel.extent_int(1); ++j) {
      kernel_host(i, j) = kernel_host(i, j - 1) * ComplexT(0, (i - n) * k0);
    }
  }

  Kokkos::deep_copy(kernel, kernel_host);

  return kernel;
}

template<class T>
__global__ void
prepare_cutoff_kernel(Kokkos::complex<T>* kernel, int n, int cutoff_idx)
{
  auto i = int(blockIdx.x * blockDim.x + threadIdx.x);

  // Set the Fourier coefficients of the kernel
  if (i < n) {
    if (i >= cutoff_idx && i <= n - cutoff_idx) {
      kernel[i] = Kokkos::complex<T>(0, 0);
    } else {
      kernel[i] = Kokkos::complex<T>(1, 0);
    }
  }
}

bool check_vkfft_plan(VkFFTApplication const* app,
                      VkFFTPlan const*        plan,
                      Logger                  logger)
{
  if (plan->bigSequenceEvenR2C != 0u) {
    logger->error("Requested FFT size ({}) unsupported by vkFFT.",
                  app->configuration.size[0]);
    return false;
  }
  return true;
}

} // namespace

template class SpectralPlan<float, Kokkos::DefaultExecutionSpace, fft::VKFFT>;
template class SpectralPlan<double, Kokkos::DefaultExecutionSpace, fft::VKFFT>;

template<>
std::unique_ptr<SpectralPlanBase<float, Kokkos::DefaultExecutionSpace>>
SpectralPlanFactory::create<float, Kokkos::DefaultExecutionSpace, fft::VKFFT>(
  PencilPlan const&      partition,
  float                  kx0,
  float                  ky0,
  SpectralOptions const& options)
{
  return std::make_unique<
    SpectralPlan<float, Kokkos::DefaultExecutionSpace, fft::VKFFT>>(
    partition, kx0, ky0, options);
}
template<>
std::unique_ptr<SpectralPlanBase<double, Kokkos::DefaultExecutionSpace>>
SpectralPlanFactory::create<double, Kokkos::DefaultExecutionSpace, fft::VKFFT>(
  PencilPlan const&      partition,
  double                 kx0,
  double                 ky0,
  SpectralOptions const& options)
{
  return std::make_unique<
    SpectralPlan<double, Kokkos::DefaultExecutionSpace, fft::VKFFT>>(
    partition, kx0, ky0, options);
}

} // namespace alps::spectral
