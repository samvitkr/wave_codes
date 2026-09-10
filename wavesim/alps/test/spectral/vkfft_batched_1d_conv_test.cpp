// general parts
#include <algorithm>
#include <cstdio>
#include <memory>
#include <string>
#include <vector>

#include <catch2/catch_session.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <fmt/base.h>

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

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wdouble-promotion"
#ifdef __NVCC_DIAG_PRAGMA_SUPPORT__
#pragma nv_diagnostic push
#pragma nv_diag_suppress 68
#pragma nv_diag_suppress 177
#pragma nv_diag_suppress 550
#else
#pragma diagnostic push
#pragma diag_suppress 68
#pragma diag_suppress 177
#pragma diag_suppress 550
#endif
#include <vkFFT.h>
#ifdef __NVCC_DIAG_PRAGMA_SUPPORT__
#pragma nv_diagnostic pop
#else
#pragma diagnostic pop
#endif
#pragma GCC diagnostic pop

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
using VkApp = std::unique_ptr<VkFFTApplication, void (*)(VkFFTApplication*)>;

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
  // Force context initialization before CUDA 12
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
VkApp setup_batched_1d_r2c_convolution_VkFFT(Device_t              device,
                                             int                   size,
                                             int                   n_batch,
                                             DevVector<floatT>&    input_dev,
                                             DevVector<floatT>&    output_dev,
                                             DevVector<std::byte>& buffer_dev)
{
  using complexT = ComplexT<floatT>;

  VkFFTConfiguration configuration = {};

  // Multidimensional FFT dimensions sizes
  configuration.FFTdim  = 1; // FFT dimension, 1D, 2D or 3D (default 1).
  configuration.size[0] = size;
  configuration.size[1] = 1;
  configuration.size[2] = 1;

  configuration.doublePrecision    = std::is_same_v<floatT, float> ? 0 : 1;
  configuration.performR2C         = true;
  configuration.coordinateFeatures = 1;
  configuration.normalize          = 1; // normalize iFFT

  // After this, configuration file contains pointers to Vulkan objects needed
  // to work with the GPU: VkDevice* device - created device, [uint64_t
  // *bufferSize, VkBuffer *buffer, VkDeviceMemory* bufferDeviceMemory] -
  // allocated GPU memory FFT is performed on. [uint64_t *kernelSize, VkBuffer
  // *kernel, VkDeviceMemory* kernelDeviceMemory] - allocated GPU memory,
  // where kernel for convolution is stored.
  configuration.device = &device;

  configuration.kernelConvolution           = false;
  configuration.performConvolution          = true;
  configuration.numberBatches               = n_batch;
  configuration.singleKernelMultipleBatches = true;
  configuration.isInputFormatted            = true;
  configuration.isOutputFormatted           = true;

  // Allocate separate buffer for the input data.
  uint64_t requiredbufferSize = configuration.numberBatches
                              * (configuration.size[0] / 2 + 1)
                              * sizeof(complexT);
  if (buffer_dev.size() < requiredbufferSize) {
    buffer_dev.resize(requiredbufferSize);
  }
  uint64_t kernelSize = sizeof(complexT) * size;

  void* bufferPtr       = buffer_dev.data().get();
  void* inputBufferPtr  = input_dev.data().get();
  void* outputBufferPtr = output_dev.data().get();

  configuration.inputBufferStride[0]  = configuration.size[0];
  configuration.inputBuffer           = &inputBufferPtr;
  configuration.outputBufferStride[0] = configuration.size[0];
  configuration.outputBuffer          = &outputBufferPtr;
  configuration.bufferSize            = &requiredbufferSize;
  configuration.buffer                = &bufferPtr;
  configuration.kernelSize            = &kernelSize;
  configuration.kernel                = nullptr;
  configuration.numberKernels         = 1;

  // Initialize application responsible for the convolution.
  auto app = VkApp{new VkFFTApplication(), &deleteVkFFT};
  // configuration.keepShaderCode = 1;
  auto resFFT = initializeVkFFT(app.get(), configuration);
  if (resFFT != VKFFT_SUCCESS) {
    std::string error_msg{getVkFFTErrorString(resFFT)};
    throw std::runtime_error("Failed to initialize VkFFT: " + error_msg);
  }
  if (app->localFFTPlan_inverse->bigSequenceEvenR2C != 0) {
    throw std::runtime_error("Failed to initialize VkFFT: size too large");
  }

  return std::move(app);
}

VkFFTResult execute_fft(VkApp const& app, VkFFTLaunchParams launchParams)
{
  auto resFFT = VkFFTAppend(app.get(), -1, &launchParams);
  if (resFFT != VKFFT_SUCCESS) {
    std::string error_msg{getVkFFTErrorString(resFFT)};
    throw std::runtime_error("Forward FFT failed: " + error_msg);
  }
  synchronize_vkfft(app.get());
  return resFFT;
}

template<typename floatT>
std::pair<double, double> check_result(DevVector<floatT> const&   dev_output,
                                       std::vector<floatT> const& ref_output)
{
  // Transfer data from GPU using staging buffer.
  std::vector<floatT> buffer_output(dev_output.size());
  thrust::copy(dev_output.begin(), dev_output.end(), buffer_output.begin());

  // Check result
  double sum2    = 0;
  double max_abs = 0;
  for (auto i = 0u; i < ref_output.size(); ++i) {
    auto err = std::abs(buffer_output[i] - ref_output[i]);
    sum2 += err * err;
    max_abs = max_abs < err ? err : max_abs;
  }

  return {max_abs, std::sqrt(sum2 / ref_output.size())};
}
} // anonymous namespace

TEMPLATE_TEST_CASE("Batched 1D convolution", "[VkFFT]", double, float)
{
  std::vector<int> test_sizes{32,   40,   42,   48,   60,   64,   72,   80,
                              96,   128,  160,  144,  192,  224,  256,  288,
                              320,  384,  420,  512,  576,  640,  768,  896,
                              1024, 1152, 1280, 1536, 1792, 2048, 4096, 5120};
  std::vector<int> test_batches{1, 2, 3, 6, 7, 9, 12};
  std::copy_if(test_sizes.begin(),
               test_sizes.end(),
               std::back_inserter(test_batches),
               [](int const n) { return n <= 1024; });

  Device_t device{};
#if (VKFFT_BACKEND == 1)
  cuDeviceGet(&device, DEVICE_ID);
#elif (VKFFT_BACKEND == 2)
  hipDeviceGet(&device, DEVICE_ID);
#endif

  for (auto transform_size : test_sizes) {
    for (auto n_batch : test_batches) {
      DYNAMIC_SECTION("transform size: " << transform_size << " x " << n_batch)
      {
        auto ref_data = generate_data<TestType>(transform_size, n_batch);

        auto input_dev  = DevVector<TestType>(ref_data.input);
        auto output_dev = DevVector<TestType>(input_dev.size());
        auto kernel_dev = DevVector<ComplexT<TestType>>(ref_data.kernel);
        auto buffer_dev = DevVector<std::byte>{};

        auto app = setup_batched_1d_r2c_convolution_VkFFT(
          device, transform_size, n_batch, input_dev, output_dev, buffer_dev);

        void* inputBufferPtr  = input_dev.data().get();
        void* outputBufferPtr = output_dev.data().get();
        void* kernelPtr       = kernel_dev.data().get();
        void* bufferPtr       = buffer_dev.data().get();
        auto  launchParams    = [&] {
          VkFFTLaunchParams params{};
          params.inputBuffer  = &inputBufferPtr;
          params.outputBuffer = &outputBufferPtr;
          params.kernel       = &kernelPtr;
          params.buffer       = &bufferPtr;
          return params;
        }();

        execute_fft(app, launchParams);

        auto [Linf_error, L2_error] = check_result(output_dev, ref_data.output);
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

int main(int argc, char* argv[])
{
  int result = prepareDevice();
  if (result != 0) {
    return result;
  }

  return Catch::Session().run(argc, argv);
}
