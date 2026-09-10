#include "fftplan_cufft.h"

#include <common/base/logging.h>
#include <common/memory/dynamic_size_pool.h>
#include <common/runtime/manager.h>

#include <cuda.h>
#include <cufft.h>
#include <fmt/ranges.h>

#include <utility>

namespace Kokkos {
class CudaSpace;
} // namespace Kokkos

namespace alps::fft {

template<class T, class TransformType>
struct CUFFTInterface;

template<>
struct CUFFTInterface<float, R2C>
{
  using handle_t = cufftHandle;

  static constexpr auto transform_type = CUFFT_R2C;
  using input_t                        = cufftReal;
  using output_t                       = cufftComplex;
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
  static constexpr auto& exec = cufftExecR2C;
};

template<>
struct CUFFTInterface<double, R2C>
{
  using handle_t = cufftHandle;

  static constexpr auto transform_type = CUFFT_D2Z;
  using input_t                        = cufftDoubleReal;
  using output_t                       = cufftDoubleComplex;
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
  static constexpr auto& exec = cufftExecD2Z;
};

template<>
struct CUFFTInterface<float, C2R>
{
  using handle_t = cufftHandle;

  static constexpr auto transform_type = CUFFT_C2R;
  using input_t                        = cufftComplex;
  using output_t                       = cufftReal;
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
  static constexpr auto& exec = cufftExecC2R;
};

template<>
struct CUFFTInterface<double, C2R>
{
  using handle_t = cufftHandle;

  static constexpr auto transform_type = CUFFT_Z2D;
  using input_t                        = cufftDoubleComplex;
  using output_t                       = cufftDoubleReal;
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
  static constexpr auto& exec = cufftExecZ2D;
};

using cufft_result_id_t = std::underlying_type_t<cufftResult_t>;

template<class T, class TransformType>
FFTPlan<T, TransformType, CUFFT>::FFTPlan(const handle_t cufft_handle)
  : plan_{cufft_handle}
  , initialized_{true}
  , stream_{[](Logger logger) -> cudaStream_t {
    cudaStream_t stream{};
    if (auto result = cudaStreamCreate(&stream); result != cudaSuccess) {
      logger->warn("cudaStreamCreate error: {}", cudaGetErrorString(result));
      return nullptr;
    }
    return stream;
  }(logger_)}
  , work_size_{work_size()}
{
  cufftSetStream(plan_, stream_);

  // Set a work area to the plan to prevent errors
  auto& allocator = RuntimeManager::instance().memory_pool<Kokkos::CudaSpace>();
  auto* work_area = allocator.allocate(128u);
  cufftSetWorkArea(plan_, work_area);
  allocator.deallocate(work_area);
}

template<class T, class TransformType>
FFTPlan<T, TransformType, CUFFT>&
FFTPlan<T, TransformType, CUFFT>::operator=(FFTPlan&& src) noexcept
{
  if (this != &src) {
    // Destroy the old plan to prepare moving the plan from src to this.
    if (initialized_) {
      cuStreamSynchronize(stream_);
      if (auto result = cufftDestroy(plan_); CUFFT_SUCCESS != result) {
        logger_->warn("Failed to destroy a cuFFT plan with error {}",
                      static_cast<cufft_result_id_t>(result));
      }
    }
    plan_        = std::move(src.plan_);
    stream_      = std::move(src.stream_);
    initialized_ = std::move(src.initialized_);

    src.plan_        = {};
    src.initialized_ = false;
    src.stream_      = nullptr;

    if (initialized_) {
      work_size_ = work_size();
    }
  }
  return *this;
}

template<class T, class TransformType>
FFTPlan<T, TransformType, CUFFT>::FFTPlan(FFTPlan&& src) noexcept
  : plan_{std::move(src.plan_)}
  , initialized_{std::move(src.initialized_)}
  , stream_{std::move(src.stream_)}
  , work_size_{work_size()}
{
  src.plan_        = {};
  src.initialized_ = false;
  src.stream_      = nullptr;
}

template<class T, class TransformType>
FFTPlan<T, TransformType, CUFFT>::~FFTPlan()
{
  // Destroy the cufft plan if it is initialized.
  if (initialized_) {
    cuStreamSynchronize(stream_);
    if (auto result = cufftDestroy(plan_); CUFFT_SUCCESS != result) {
      logger_->warn("Failed to destroy a cuFFT plan with error {}",
                    static_cast<cufft_result_id_t>(result));
    }
  }

  // Destroy the stream if not the default stream.
  if (stream_ != nullptr) {
    if (auto result = cudaStreamDestroy(stream_); result != cudaSuccess) {
      logger_->warn("cudaStreamDestroy error: {}", cudaGetErrorString(result));
    }
  }
}

template<class T, class TransformType>
std::size_t FFTPlan<T, TransformType, CUFFT>::work_size() const noexcept
{
  if (!initialized_) {
    return 0;
  }

  std::size_t size{0};
  if (auto result = cufftGetSize(plan_, &size); result != CUFFT_SUCCESS) {
    logger_->error("cufftGetSize failed with error {}",
                   static_cast<cufft_result_id_t>(result));
  }
  return size;
}

namespace {
struct workarea_alloc_data
{
  void*       ptr;
  std::size_t work_size;
  cufftHandle plan;
};

void CUDART_CB allocate_workarea(void* args)
{
  auto* alloc_data = static_cast<workarea_alloc_data*>(args);

  auto& allocator = RuntimeManager::instance().memory_pool<Kokkos::CudaSpace>();
  alloc_data->ptr = allocator.allocate(alloc_data->work_size);
  cufftSetWorkArea(alloc_data->plan, alloc_data->ptr);
}

void CUDART_CB deallocate_workarea(void* args)
{
  auto* alloc_data = static_cast<workarea_alloc_data*>(args);
  auto  ptr_to_alloc_data =
    ::std::unique_ptr<workarea_alloc_data>(alloc_data); // ensure the memory is
                                                        // deallocated
  auto& allocator = RuntimeManager::instance().memory_pool<Kokkos::CudaSpace>();
  allocator.deallocate(alloc_data->ptr);
}
} // namespace

template<class T, class TransformType>
void FFTPlan<T, TransformType, CUFFT>::run(void const*  input,
                                           void*        output,
                                           cudaStream_t stream) const
{
  // Reserve memory in the pool to avoid cudaMalloc being called in cuda
  // callback function
  if (auto& allocator =
        RuntimeManager::instance().memory_pool<Kokkos::CudaSpace>();
      allocator.getLargestAvailableBlock() < work_size_) {
    logger_->debug("The memory pool is not large enough to allocate the "
                   "work area for the FFT plan. {} bytes will be reserved.",
                   work_size_);
    void* tmp = allocator.allocate(work_size_);
    allocator.deallocate(tmp);
  }

  cudaStream_t ext_stream = stream;

  cudaEvent_t event{};
  cudaEventCreateWithFlags(&event, cudaEventDisableTiming);

  // Enqueue the event in the external stream.
  // The execution of the plan is blocked until the external stream
  // reaches the event.
  cudaEventRecord(event, ext_stream);
  cudaStreamWaitEvent(stream_, event);

  // Allocate a work area.
  auto* workload = new workarea_alloc_data{nullptr, work_size_, plan_};
  cudaLaunchHostFunc(stream_, allocate_workarea, (void*)workload);

  using input_t  = typename CUFFTInterface<T, TransformType>::input_t;
  using output_t = typename CUFFTInterface<T, TransformType>::output_t;
  auto result    = CUFFTInterface<T, TransformType>::exec(
    plan_, (input_t*)input, (output_t*)output);

  cudaLaunchHostFunc(stream_, deallocate_workarea, (void*)workload);

  // Enqueue the finish event in the external stream to synchronize the
  // external stream with the fft
  cudaEventRecord(event, stream_);
  cudaStreamWaitEvent(ext_stream, event);

  cudaEventDestroy(event);

  if (result != CUFFT_SUCCESS) {
    logger_->error("CUFFT error {}", static_cast<cufft_result_id_t>(result));
  }
}

// Explicit instantiation
template class FFTPlan<float, R2C, CUFFT>;
template class FFTPlan<float, C2R, CUFFT>;
template class FFTPlan<double, R2C, CUFFT>;
template class FFTPlan<double, C2R, CUFFT>;

} // namespace alps::fft
