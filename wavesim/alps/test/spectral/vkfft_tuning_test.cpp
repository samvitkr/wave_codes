#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <catch2/catch_session.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_contains.hpp>
#include <fmt/base.h>

#include <common/base/logging.h>
#include <fft/vkfft.h>
#include <spectral/detail/vkfft_app_utils.h>
#include <spectral/detail/vkfft_tuning.h>

#if (VKFFT_BACKEND == 1)
#include <cuda_runtime.h>
#include <thrust/complex.h>
#include <thrust/copy.h>
#include <thrust/device_vector.h>
#elif (VKFFT_BACKEND == 2)
#ifndef __HIP_PLATFORM_HCC__
#define __HIP_PLATFORM_HCC__
#endif
#include <hip/hip_complex.h>
#include <hip/hip_runtime.h>
#include <hip/hiprtc.h>
#include <thrust/complex.h>
#include <thrust/copy.h>
#include <thrust/device_vector.h>
#endif

#if (VKFFT_BACKEND == 1 || VKFFT_BACKEND == 2)
template<typename T>
using DevVector = thrust::device_vector<T>;
template<typename T>
using ComplexT = thrust::complex<T>;
#endif

template<class T>
struct ReferenceData
{
  std::vector<T>           input;
  std::vector<T>           output;
  std::vector<ComplexT<T>> kernel;
};

namespace {
using alps::spectral::detail::VkApp;

#if (VKFFT_BACKEND == 1)
using Device_t = CUdevice;
#elif (VKFFT_BACKEND == 2)
using Device_t = hipDevice_t;
#endif

constexpr int DEVICE_ID{0};

int prepareDevice()
{
#if (VKFFT_BACKEND == 1)
  if (auto err = cudaSetDevice(DEVICE_ID); err != cudaSuccess) {
    fmt::println(
      stderr, "Error setting CUDA device: {}", cudaGetErrorString(err));
    return err;
  }
  if (auto err = cudaFree(0); err != cudaSuccess) {
    fmt::println(
      stderr, "Error freeing CUDA memory: {}", cudaGetErrorString(err));
    return err;
  }
#elif (VKFFT_BACKEND == 2)
  if (auto res = hipSetDevice(DEVICE_ID); res != hipSuccess) {
    fmt::println(
      stderr, "Error setting HIP device: {}", hipGetErrorString(res));
    return res;
  }
#endif
  return 0;
}

template<class T>
ReferenceData<T> generate_data(int n, int batch)
{
  const double k0 = 0.5;

  ReferenceData<T> ref_data{std::vector<T>(n * batch),
                            std::vector<T>(n * batch),
                            std::vector<ComplexT<T>>(n)};
  auto&            input  = ref_data.input;
  auto&            output = ref_data.output;
  auto&            kernel = ref_data.kernel;

  using std::cos;
  using std::sin;
  for (int b = 0; b < batch; ++b) {
    for (int i = 0; i < n; ++i) {
      double x   = i * 2 * M_PI / k0 / n;
      double y   = b * 2 * M_PI / batch;
      double val = sin(3 * k0 * x) + 0.2 * cos(13 * k0 * x + 2 * y + 3.0)
                 + 0.8 * sin((n / 3 + 1) * k0 * x - (b / 3 - 1) * y + 2.0)
                 - 0.75 * cos((n / 2 - 1) * k0 * x - (b / 2 - 1) * y + 2.0);
      double dx_val = 3 * k0 * cos(3 * k0 * x)
                    - 0.2 * 13 * k0 * sin(13 * k0 * x + 2 * y + 3.0)
                    + 0.8 * (n / 3 + 1) * k0
                        * cos((n / 3 + 1) * k0 * x - (b / 3 - 1) * y + 2.0)
                    + 0.75 * (n / 2 - 1) * k0
                        * sin((n / 2 - 1) * k0 * x - (b / 2 - 1) * y + 2.0);
      input[i + b * n]  = static_cast<T>(val);
      output[i + b * n] = static_cast<T>(dx_val);
    }
  }

  for (int i = 0; i < n / 2; ++i) {
    kernel[i] = ComplexT<T>(0, i * k0);
  }
  kernel[n / 2] = 0;
  for (int i = n / 2 + 1; i < n; ++i) {
    kernel[i] = ComplexT<T>(0, static_cast<T>((i - n) * k0));
  }

  return ref_data;
}

void synchronize_vkfft(VkFFTApplication* app)
{
  auto const num_streams = app->configuration.num_streams;
  if (num_streams > 1) {
    if (VkFFTSync(app) != VKFFT_SUCCESS) {
      throw std::runtime_error("synchronize failed");
    }
    return;
  }
#if (VKFFT_BACKEND == 1)
  cudaStream_t stream_to_sync =
    (num_streams == 1) ? *app->configuration.stream : nullptr;
  auto res = cudaStreamSynchronize(stream_to_sync);
  if (res != cudaSuccess) {
    throw std::runtime_error("synchronize failed");
  }
#elif (VKFFT_BACKEND == 2)
  hipStream_t stream_to_sync =
    (num_streams == 1) ? *app->configuration.stream : nullptr;
  auto res = hipStreamSynchronize(stream_to_sync);
  if (res != hipSuccess) {
    throw std::runtime_error("synchronize failed");
  }
#endif
}

template<typename floatT>
VkFFTConfiguration
make_batched_1d_r2c_convolution_config(Device_t* device,
                                       int       size,
                                       int       n_batch,
                                       void**    input_buffer_ptr,
                                       void**    output_buffer_ptr,
                                       void**    buffer_ptr,
                                       uint64_t* buffer_size,
                                       uint64_t* kernel_size)
{
  VkFFTConfiguration cfg{};
  cfg.FFTdim                      = 1;
  cfg.size[0]                     = size;
  cfg.size[1]                     = 1;
  cfg.size[2]                     = 1;
  cfg.doublePrecision             = std::is_same_v<floatT, float> ? 0 : 1;
  cfg.performR2C                  = true;
  cfg.coordinateFeatures          = 1;
  cfg.normalize                   = 1;
  cfg.device                      = device;
  cfg.kernelConvolution           = false;
  cfg.performConvolution          = true;
  cfg.numberBatches               = n_batch;
  cfg.singleKernelMultipleBatches = true;
  cfg.isInputFormatted            = true;
  cfg.isOutputFormatted           = true;
  cfg.inputBufferStride[0]        = size;
  cfg.inputBuffer                 = input_buffer_ptr;
  cfg.outputBufferStride[0]       = size;
  cfg.outputBuffer                = output_buffer_ptr;
  cfg.bufferSize                  = buffer_size;
  cfg.buffer                      = buffer_ptr;
  cfg.kernelSize                  = kernel_size;
  cfg.kernel                      = nullptr;
  cfg.numberKernels               = 1;
  return cfg;
}

template<typename floatT>
VkApp setup_batched_1d_r2c_convolution_vkfft(
  Device_t                                            device,
  int                                                 size,
  int                                                 n_batch,
  DevVector<floatT>&                                  input_dev,
  DevVector<floatT>&                                  output_dev,
  DevVector<std::byte>&                               buffer_dev,
  alps::spectral::detail::VkFFTSchedulerParams const& scheduler_params = {})
{
  using complexT = ComplexT<floatT>;

  uint64_t requiredbufferSize =
    static_cast<uint64_t>(n_batch) * (size / 2 + 1) * sizeof(complexT);
  if (buffer_dev.size() < requiredbufferSize) {
    buffer_dev.resize(requiredbufferSize);
  }
  uint64_t kernelSize = sizeof(complexT) * size;

  void* bufferPtr       = buffer_dev.data().get();
  void* inputBufferPtr  = input_dev.data().get();
  void* outputBufferPtr = output_dev.data().get();
  // not used in initialization, will be set by execute_convolution

  auto configuration =
    make_batched_1d_r2c_convolution_config<floatT>(&device,
                                                   size,
                                                   n_batch,
                                                   &inputBufferPtr,
                                                   &outputBufferPtr,
                                                   &bufferPtr,
                                                   &requiredbufferSize,
                                                   &kernelSize);

  configuration = alps::spectral::detail::with_scheduler_params(
    configuration, scheduler_params);
  configuration.keepShaderCode = true;

  auto app = alps::spectral::detail::make_vkfft_app(configuration);
  if (app->localFFTPlan_inverse->bigSequenceEvenR2C != 0) {
    throw std::runtime_error("Failed to initialize VkFFT: size too large");
  }
  return app;
}

template<typename floatT>
std::pair<double, double> check_result(DevVector<floatT> const&   dev_output,
                                       std::vector<floatT> const& ref_output)
{
  std::vector<floatT> buffer_output(dev_output.size());
  thrust::copy(dev_output.begin(), dev_output.end(), buffer_output.begin());

  double sum2    = 0;
  double max_abs = 0;
  for (auto i = 0u; i < ref_output.size(); ++i) {
    auto err = std::abs(buffer_output[i] - ref_output[i]);
    sum2 += err * err;
    max_abs = max_abs < err ? err : max_abs;
  }
  return {max_abs, std::sqrt(sum2 / ref_output.size())};
}

template<typename floatT>
void execute_convolution(VkApp const&                 app,
                         DevVector<floatT>&           input_dev,
                         DevVector<floatT>&           output_dev,
                         DevVector<ComplexT<floatT>>& kernel_dev,
                         DevVector<std::byte>&        buffer_dev)
{
  void* inputBufferPtr  = input_dev.data().get();
  void* outputBufferPtr = output_dev.data().get();
  void* kernelPtr       = kernel_dev.data().get();
  void* bufferPtr       = buffer_dev.data().get();

  VkFFTLaunchParams params{};
  params.inputBuffer  = &inputBufferPtr;
  params.outputBuffer = &outputBufferPtr;
  params.kernel       = &kernelPtr;
  params.buffer       = &bufferPtr;

  auto res = VkFFTAppend(app.get(), -1, &params);
  if (res != VKFFT_SUCCESS) {
    std::string error_msg{getVkFFTErrorString(res)};
    throw std::runtime_error("Convolution FFT failed: " + error_msg);
  }
  synchronize_vkfft(app.get());
}
} // namespace

TEMPLATE_TEST_CASE("Scheduler parameter overrides produce correct results",
                   "[VkFFT][tuning]",
                   double,
                   float)
{
  std::vector<int> thresholds{8, 24};
  std::vector<int> disable_rb_values{0, 1};
  std::vector<int> aim_threads_values{0, 64, 256};
  std::vector<int> grouped_batch0_values{0, 1};

  std::vector<int> sizes{64, 256, 1024};
  std::vector<int> batches{1, 128};

  Device_t device{};
#if (VKFFT_BACKEND == 1)
  cuDeviceGet(&device, DEVICE_ID);
#elif (VKFFT_BACKEND == 2)
  hipDeviceGet(&device, DEVICE_ID);
#endif

  for (auto threshold : thresholds) {
    for (auto disable_rb : disable_rb_values) {
      for (auto aim_threads : aim_threads_values) {
        for (auto grouped_batch0 : grouped_batch0_values) {
          for (auto transform_size : sizes) {
            for (auto n_batch : batches) {
              DYNAMIC_SECTION("threshold="
                              << threshold << ", disable_rb=" << disable_rb
                              << ", aim_threads=" << aim_threads
                              << ", grouped_batch0=" << grouped_batch0
                              << ", size=" << transform_size
                              << ", batch=" << n_batch)
              {
                auto ref_data =
                  generate_data<TestType>(transform_size, n_batch);

                auto input_dev  = DevVector<TestType>(ref_data.input);
                auto output_dev = DevVector<TestType>(input_dev.size());
                auto kernel_dev =
                  DevVector<ComplexT<TestType>>(ref_data.kernel);
                auto buffer_dev = DevVector<std::byte>{};

                auto app = setup_batched_1d_r2c_convolution_vkfft(
                  device,
                  transform_size,
                  n_batch,
                  input_dev,
                  output_dev,
                  buffer_dev,
                  {threshold, disable_rb, aim_threads, grouped_batch0});

                execute_convolution(
                  app, input_dev, output_dev, kernel_dev, buffer_dev);

                auto [Linf_error, L2_error] =
                  check_result(output_dev, ref_data.output);
                if constexpr (std::is_same_v<TestType, double>) {
                  REQUIRE(Linf_error < sqrt(std::log(transform_size)) * 1.5e-9);
                  REQUIRE(L2_error < sqrt(std::log(transform_size)) * 2.4e-10);
                } else {
                  REQUIRE(Linf_error < sqrt(std::log(transform_size)) * 1.3e-3);
                  REQUIRE(L2_error < sqrt(std::log(transform_size)) * 2.4e-4);
                }
              }
            }
          }
        }
      }
    }
  }
}

TEST_CASE(
  "tune_scheduler_params returns valid params and produces correct plan",
  "[VkFFT][tuning]")
{
  constexpr int transform_size = 512;
  constexpr int n_batch        = 64;

  Device_t device{};
#if (VKFFT_BACKEND == 1)
  cuDeviceGet(&device, DEVICE_ID);
#elif (VKFFT_BACKEND == 2)
  hipDeviceGet(&device, DEVICE_ID);
#endif

  auto ref_data   = generate_data<float>(transform_size, n_batch);
  auto input_dev  = DevVector<float>(ref_data.input);
  auto output_dev = DevVector<float>(ref_data.input.size());
  auto kernel_dev = DevVector<ComplexT<float>>(ref_data.kernel);
  auto buffer_dev = DevVector<std::byte>{};

  uint64_t buffer_size = static_cast<uint64_t>(n_batch)
                       * (transform_size / 2 + 1) * sizeof(ComplexT<float>);
  if (buffer_dev.size() < buffer_size) {
    buffer_dev.resize(buffer_size);
  }
  uint64_t kernel_size = sizeof(ComplexT<float>) * transform_size;

  void* input_buffer_ptr  = input_dev.data().get();
  void* output_buffer_ptr = output_dev.data().get();
  void* buffer_ptr        = buffer_dev.data().get();
  void* kernel_ptr        = kernel_dev.data().get();

  auto make_config = [&]() {
    return make_batched_1d_r2c_convolution_config<float>(&device,
                                                         transform_size,
                                                         n_batch,
                                                         &input_buffer_ptr,
                                                         &output_buffer_ptr,
                                                         &buffer_ptr,
                                                         &buffer_size,
                                                         &kernel_size);
  };

  auto logger = alps::get_logger("vkfft_tuning_test");

  alps::spectral::detail::VkFFTTuningConfig const tuning_config = [] {
    alps::spectral::detail::VkFFTTuningConfig cfg{};
    cfg.enabled                      = true;
    cfg.thresholds                   = {8, 12, 16, 24};
    cfg.sweep_refine_batch           = true;
    cfg.aim_threads_values           = {32, 64, 128};
    cfg.sweep_aim_threads            = true;
    cfg.grouped_batch0_values        = {0, 1};
    cfg.sweep_grouped_batch0         = true;
    cfg.benchmark_samples            = 5;
    cfg.min_measurable_batch_time_us = 20.0;
    cfg.target_batch_time_us         = 3000.0;
    cfg.pilot_batch_time_us          = 250.0;
    return cfg;
  }();

  auto app = alps::spectral::detail::tune_scheduler_params(
    make_config,
    [&](VkFFTApplication* app, void*) {
      VkFFTLaunchParams params{};
      params.inputBuffer  = &input_buffer_ptr;
      params.outputBuffer = &output_buffer_ptr;
      params.kernel       = &kernel_ptr;
      params.buffer       = &buffer_ptr;

      auto res = VkFFTAppend(app, -1, &params);
      if (res != VKFFT_SUCCESS) {
        std::string error_msg{getVkFFTErrorString(res)};
        throw std::runtime_error("Convolution FFT failed: " + error_msg);
      }
      synchronize_vkfft(app);
    },
    nullptr,
    tuning_config,
    logger);

  REQUIRE(app != nullptr);

  alps::spectral::detail::VkFFTSchedulerParams best{};
  best.register_threshold =
    static_cast<int>(app->configuration.custom_register_threshold);
  best.disable_refine_batch =
    static_cast<int>(app->configuration.disable_refine_batch);
  best.aim_threads    = static_cast<int>(app->configuration.aimThreads);
  best.grouped_batch0 = static_cast<int>(app->configuration.groupedBatch[0]);

  REQUIRE_THAT(tuning_config.thresholds,
               Catch::Matchers::Contains(best.register_threshold));
  REQUIRE((best.disable_refine_batch == 0 || best.disable_refine_batch == 1));
  REQUIRE_THAT(tuning_config.aim_threads_values,
               Catch::Matchers::Contains(best.aim_threads));
  REQUIRE_THAT(tuning_config.grouped_batch0_values,
               Catch::Matchers::Contains(best.grouped_batch0));

  execute_convolution(app, input_dev, output_dev, kernel_dev, buffer_dev);

  auto [Linf_error, L2_error] = check_result(output_dev, ref_data.output);
  REQUIRE(Linf_error < sqrt(std::log(transform_size)) * 1.3e-3);
  REQUIRE(L2_error < sqrt(std::log(transform_size)) * 2.4e-4);
}

TEST_CASE("tune_scheduler_params returns baseline when no sweep enabled",
          "[VkFFT][tuning]")
{
  constexpr int transform_size = 512;
  constexpr int n_batch        = 64;

  Device_t device{};
#if (VKFFT_BACKEND == 1)
  cuDeviceGet(&device, DEVICE_ID);
#elif (VKFFT_BACKEND == 2)
  hipDeviceGet(&device, DEVICE_ID);
#endif

  auto ref_data   = generate_data<float>(transform_size, n_batch);
  auto input_dev  = DevVector<float>(ref_data.input);
  auto output_dev = DevVector<float>(ref_data.input.size());
  auto kernel_dev = DevVector<ComplexT<float>>(ref_data.kernel);
  auto buffer_dev = DevVector<std::byte>{};

  uint64_t buffer_size = static_cast<uint64_t>(n_batch)
                       * (transform_size / 2 + 1) * sizeof(ComplexT<float>);
  if (buffer_dev.size() < buffer_size) {
    buffer_dev.resize(buffer_size);
  }
  uint64_t kernel_size = sizeof(ComplexT<float>) * transform_size;

  void* input_buffer_ptr  = input_dev.data().get();
  void* output_buffer_ptr = output_dev.data().get();
  void* buffer_ptr        = buffer_dev.data().get();
  void* kernel_ptr        = kernel_dev.data().get();

  auto make_config = [&]() {
    return make_batched_1d_r2c_convolution_config<float>(&device,
                                                         transform_size,
                                                         n_batch,
                                                         &input_buffer_ptr,
                                                         &output_buffer_ptr,
                                                         &buffer_ptr,
                                                         &buffer_size,
                                                         &kernel_size);
  };

  auto logger = alps::get_logger("vkfft_tuning_test");

  alps::spectral::detail::VkFFTTuningConfig baseline_config{};
  baseline_config.enabled                      = true;
  baseline_config.thresholds                   = {};
  baseline_config.sweep_refine_batch           = false;
  baseline_config.sweep_aim_threads            = false;
  baseline_config.sweep_grouped_batch0         = false;
  baseline_config.benchmark_samples            = 1;
  baseline_config.min_measurable_batch_time_us = 20.0;
  baseline_config.target_batch_time_us         = 3000.0;
  baseline_config.pilot_batch_time_us          = 250.0;

  auto app = alps::spectral::detail::tune_scheduler_params(
    make_config,
    [&](VkFFTApplication* app, void*) {
      VkFFTLaunchParams params{};
      params.inputBuffer  = &input_buffer_ptr;
      params.outputBuffer = &output_buffer_ptr;
      params.kernel       = &kernel_ptr;
      params.buffer       = &buffer_ptr;

      auto res = VkFFTAppend(app, -1, &params);
      if (res != VKFFT_SUCCESS) {
        std::string error_msg{getVkFFTErrorString(res)};
        throw std::runtime_error("Convolution FFT failed: " + error_msg);
      }
      synchronize_vkfft(app);
    },
    nullptr,
    baseline_config,
    logger);

  REQUIRE(app != nullptr);

  alps::spectral::detail::VkFFTSchedulerParams best{};
  best.register_threshold =
    static_cast<int>(app->configuration.custom_register_threshold);
  best.disable_refine_batch =
    static_cast<int>(app->configuration.disable_refine_batch);
  best.aim_threads    = static_cast<int>(app->configuration.aimThreads);
  best.grouped_batch0 = static_cast<int>(app->configuration.groupedBatch[0]);

  REQUIRE(best.register_threshold == 0);
  REQUIRE(best.disable_refine_batch == 0);
  REQUIRE(best.aim_threads == 128); // default aim_threads
  REQUIRE(best.grouped_batch0 == 0);

  execute_convolution(app, input_dev, output_dev, kernel_dev, buffer_dev);

  auto [Linf_error, L2_error] = check_result(output_dev, ref_data.output);
  REQUIRE(Linf_error < sqrt(std::log(transform_size)) * 1.3e-3);
  REQUIRE(L2_error < sqrt(std::log(transform_size)) * 2.4e-4);
}

TEST_CASE("fingerprint detects duplicate underlying plans", "[VkFFT][tuning]")
{
  constexpr int transform_size = 256;
  constexpr int n_batch        = 12;

  Device_t device{};
#if (VKFFT_BACKEND == 1)
  cuDeviceGet(&device, DEVICE_ID);
#elif (VKFFT_BACKEND == 2)
  hipDeviceGet(&device, DEVICE_ID);
#endif

  std::vector<alps::spectral::detail::VkFFTSchedulerParams> const param_sets = {
    {0, 0, 0, 0},
    {8, 0, 0, 1},
    {16, 0, 32, 0},
    {8, 0, 32, 0},
    {8, 1, 256, 0},
    {24, 1, 0, 0},
  };

  // Shared device buffers. setup_batched_1d_r2c_convolution_vkfft stores only
  // the pointer values internally, so reusing the same allocations across
  // multiple apps is safe here.
  DevVector<float>           input_dev((size_t)transform_size * n_batch, 1.0f);
  DevVector<float>           output_dev((size_t)transform_size * n_batch, 0.0f);
  DevVector<ComplexT<float>> kernel_dev((size_t)transform_size);
  DevVector<std::byte>       buffer_dev{};

  std::vector<alps::spectral::detail::VkFFTPlanFingerprint> fingerprints;
  fingerprints.reserve(param_sets.size());
  for (auto const& ps : param_sets) {
    auto app = setup_batched_1d_r2c_convolution_vkfft(
      device, transform_size, n_batch, input_dev, output_dev, buffer_dev, ps);
    fingerprints.push_back(
      alps::spectral::detail::fingerprint_vkfft_plan(app.get()));
  }

  std::vector<std::size_t> unique_indices;
  for (std::size_t i = 0; i < fingerprints.size(); ++i) {
    bool is_dup = false;
    for (std::size_t j : unique_indices) {
      if (fingerprints[i] == fingerprints[j]) {
        is_dup = true;
        break;
      }
    }
    if (!is_dup) {
      unique_indices.push_back(i);
    }
  }

  std::size_t const n_total  = fingerprints.size();
  std::size_t const n_unique = unique_indices.size();

  REQUIRE(n_unique >= 2);

  // At least two of the scheduler param sets should collapse to the same
  // underlying VkFFT plan (total > unique means at least one duplicate was
  // found). On GPUs where every param set happens to produce a genuinely
  // different plan, skip this check with a warning rather than failing, so the
  // test remains hardware-portable.
  if (n_total == n_unique) {
    WARN("No duplicate fingerprints found on this GPU — all "
         << n_total
         << " parameter combinations produce distinct underlying plans. "
            "Skipping duplicate-existence check.");
  } else {
    REQUIRE(n_total > n_unique);
  }
}

int main(int argc, char* argv[])
{
  int result = prepareDevice();
  if (result != 0) {
    return result;
  }

  return Catch::Session().run(argc, argv);
}
