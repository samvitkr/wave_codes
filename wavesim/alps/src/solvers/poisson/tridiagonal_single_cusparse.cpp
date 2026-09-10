#include "tridiagonal_single_cusparse.h"

#include <common/base/logging.h>
#include <common/kokkos_abstraction/pool_space.h>

#include <Kokkos_Core.hpp>

#include <stdexcept>

#if defined(KOKKOS_ENABLE_CUDA)
#include <cusparse.h>
#elif defined(KOKKOS_ENABLE_HIP)
#include <rocsparse/rocsparse.h>

#define MY_VERSION_NUMBER(major, minor, patch) \
  (((major) * 10000) + ((minor) * 100) + (patch))
#define MY_ROCSPARSE_VERSION \
  MY_VERSION_NUMBER(         \
    ROCSPARSE_VERSION_MAJOR, ROCSPARSE_VERSION_MINOR, ROCSPARSE_VERSION_PATCH)
#endif

namespace alps::solver {
namespace detail {
template<class T>
struct GTSVHelper;

#if defined(KOKKOS_ENABLE_CUDA)
template<>
struct GTSVHelper<double>
{
  static constexpr auto& interleavedBatchBufferSize =
    cusparseDgtsvInterleavedBatch_bufferSizeExt;
  static constexpr auto& interleavedBatch = cusparseDgtsvInterleavedBatch;
};

template<>
struct GTSVHelper<float>
{
  static constexpr auto& interleavedBatchBufferSize =
    cusparseSgtsvInterleavedBatch_bufferSizeExt;
  static constexpr auto& interleavedBatch = cusparseSgtsvInterleavedBatch;
};

template<typename T>
class GTSVImpl : public GTSVHelper<T>
{
 public:
  std::size_t      buffer_size{0u};
  cusparseHandle_t handle{};
  cudaStream_t     stream{nullptr};

  static constexpr int algorithm{0}; // cuThomas algorithm (c.f. cuSPARSE doc)

  GTSVImpl()
  {
    KOKKOS_IMPL_CUDA_SAFE_CALL(cudaStreamCreate(&stream));
    if (auto result = cusparseCreate(&handle);
        result != CUSPARSE_STATUS_SUCCESS) {
      auto msg = fmt::format("Failed to create cusparse handle: {} ({})",
                             cusparseGetErrorName(result),
                             cusparseGetErrorString(result));
      throw std::runtime_error(msg);
    }
  }

  ~GTSVImpl() noexcept
  {
    cusparseDestroy(handle);
    cudaStreamDestroy(stream);
  }
};
#elif defined(KOKKOS_ENABLE_HIP)
template<>
struct GTSVHelper<double>
{
  static constexpr auto& interleavedBatchBufferSize =
    rocsparse_dgtsv_interleaved_batch_buffer_size;
  static constexpr auto& interleavedBatch = rocsparse_dgtsv_interleaved_batch;
};

template<>
struct GTSVHelper<float>
{
  static constexpr auto& interleavedBatchBufferSize =
    rocsparse_sgtsv_interleaved_batch_buffer_size;
  static constexpr auto& interleavedBatch = rocsparse_sgtsv_interleaved_batch;
};

template<typename T>
class GTSVImpl : public GTSVHelper<T>
{
 public:
  std::size_t      buffer_size{0u};
  rocsparse_handle handle{nullptr};
  hipStream_t      stream{nullptr};

  static constexpr rocsparse_gtsv_interleaved_alg algorithm{
    rocsparse_gtsv_interleaved_alg_thomas};

  GTSVImpl()
  {
    KOKKOS_IMPL_HIP_SAFE_CALL(hipStreamCreate(&stream));
    if (auto result = rocsparse_create_handle(&handle);
        result != rocsparse_status_success) {
#if MY_ROCSPARSE_VERSION >= MY_VERSION_NUMBER(3, 0, 2)
      auto msg = fmt::format("Failed to create rocsparse handle: {} ({})",
                             rocsparse_get_status_name(result),
                             rocsparse_get_status_description(result));
#else
      auto msg = fmt::format(
        "Failed to create rocsparse handle (error {})",
        static_cast<std::underlying_type_t<rocsparse_status>>(result));
#endif
      throw std::runtime_error(msg);
    }
  }

  ~GTSVImpl() noexcept
  {
    rocsparse_destroy_handle(handle);
    handle = nullptr;

    auto result = hipStreamDestroy(stream);
    (void)result;
  }
};
#else
template<typename T>
class GTSVImpl
{
 public:
  std::size_t buffer_size{0u};
};
#endif
} // namespace detail
} // namespace alps::solver

namespace alps {
namespace solver {

template<class DataType, class... Properties>
using MDView = Kokkos::View<DataType, Kokkos::LayoutLeft, Properties...>;

template<typename ValueT>
TridiagonalCuSparse<ValueT>::TridiagonalCuSparse(
  int                        batch_count_1,
  int                        batch_count_2,
  int                        n_eqns,
  const mpipp::communicator& comm)
  : base_t(batch_count_1, batch_count_2, n_eqns, comm)
{
  this->is_coeff_overwritten = true;

#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
  if (this->nproc_ != 1) {
    throw std::runtime_error("The tridiagonal solver supports only single "
                             "process.");
  }
  impl_ = std::make_unique<impl_t>();
#else
  throw std::runtime_error("Unsupported tridiagonal solver");
#endif
}

template<typename ValueT>
TridiagonalCuSparse<ValueT>::~TridiagonalCuSparse() = default;

template<typename ValueT>
void TridiagonalCuSparse<ValueT>::setup_impl(coeff_t const& d,
                                             coeff_t const& dl,
                                             coeff_t const& du)
{
  if (d.extent_int(2) != this->nz_) {
    throw std::invalid_argument("Poisson solver coefficient d size mismatch.");
  }
  if (dl.extent_int(2) != this->nz_) {
    throw std::invalid_argument("Poisson solver coefficient dl size mismatch.");
  }
  if (du.extent_int(2) != this->nz_) {
    throw std::invalid_argument("Poisson solver coefficient du size mismatch.");
  }

  MDView<ValueT***, default_memory_pool> x(
    Kokkos::view_alloc("", Kokkos::WithoutInitializing), d.layout());
  this->batch_ = static_cast<int>(x.stride(2));

#if defined(KOKKOS_ENABLE_CUDA)
  cusparseSetStream(impl_->handle, impl_->stream);
  auto result = impl_t::interleavedBatchBufferSize(impl_->handle,
                                                   impl_t::algorithm,
                                                   this->nz_,
                                                   dl.data(),
                                                   d.data(),
                                                   du.data(),
                                                   x.data(),
                                                   this->batch_,
                                                   &impl_->buffer_size);
  if (result != CUSPARSE_STATUS_SUCCESS) {
    auto msg = fmt::format("cuSparse error {} ({})",
                           cusparseGetErrorName(result),
                           cusparseGetErrorString(result));
    throw std::runtime_error(msg);
  }
  KOKKOS_IMPL_CUDA_SAFE_CALL(cudaStreamSynchronize(impl_->stream));
#elif defined(KOKKOS_ENABLE_HIP)
  rocsparse_set_stream(impl_->handle, impl_->stream);
  auto result = impl_t::interleavedBatchBufferSize(impl_->handle,
                                                   impl_t::algorithm,
                                                   dl.extent_int(2),
                                                   dl.data(),
                                                   d.data(),
                                                   du.data(),
                                                   x.data(),
                                                   this->batch_,
                                                   this->batch_,
                                                   &impl_->buffer_size);
  if (result != rocsparse_status_success) {
#if MY_ROCSPARSE_VERSION >= MY_VERSION_NUMBER(3, 0, 2)
    auto msg = fmt::format("rocSPARSE error {} ({})",
                           rocsparse_get_status_name(result),
                           rocsparse_get_status_description(result));
#else
    auto msg = fmt::format(
      "rocSPARSE error ({})",
      static_cast<std::underlying_type_t<rocsparse_status>>(result));
#endif
    throw std::runtime_error(msg);
  }
  KOKKOS_IMPL_HIP_SAFE_CALL(hipStreamSynchronize(impl_->stream));
#endif
}

template<typename ValueT>
void TridiagonalCuSparse<ValueT>::solve_impl(solution_t const& x,
                                             coeff_t const&    d,
                                             coeff_t const&    dl,
                                             coeff_t const&    du)
{
  auto constexpr* label{"gtsv_buffer"};
  auto const      buffer_size = impl_->buffer_size;

  auto buffer = std::unique_ptr<std::byte, std::function<void(std::byte*)>>(
    (std::byte*)default_memory_pool().allocate(label, buffer_size),
    [=](std::byte* ptr) {
      default_memory_pool().deallocate(label, ptr, buffer_size);
    });

#if defined(KOKKOS_ENABLE_CUDA)
  cusparseSetStream(impl_->handle, impl_->stream);
  auto result = impl_t::interleavedBatch(impl_->handle,
                                         impl_t::algorithm,
                                         dl.extent_int(2),
                                         dl.data(),
                                         d.data(),
                                         du.data(),
                                         x.data(),
                                         this->batch_,
                                         buffer.get());
  if (result != CUSPARSE_STATUS_SUCCESS) {
    auto msg = fmt::format("cuSparse error {} ({})",
                           cusparseGetErrorName(result),
                           cusparseGetErrorString(result));
    throw std::runtime_error(msg);
  }
  KOKKOS_IMPL_CUDA_SAFE_CALL(cudaStreamSynchronize(impl_->stream));
#elif defined(KOKKOS_ENABLE_HIP)
  rocsparse_set_stream(impl_->handle, impl_->stream);
  auto result = impl_t::interleavedBatch(impl_->handle,
                                         impl_t::algorithm,
                                         dl.extent_int(2),
                                         dl.data(),
                                         d.data(),
                                         du.data(),
                                         x.data(),
                                         this->batch_,
                                         this->batch_,
                                         buffer.get());
  if (result != rocsparse_status_success) {
#if MY_ROCSPARSE_VERSION >= MY_VERSION_NUMBER(3, 0, 2)
    auto msg = fmt::format("rocSPARSE error {} ({})",
                           rocsparse_get_status_name(result),
                           rocsparse_get_status_description(result));
#else
    auto msg = fmt::format(
      "rocSPARSE error ({})",
      static_cast<std::underlying_type_t<rocsparse_status>>(result));
#endif
    throw std::runtime_error(msg);
  }
  KOKKOS_IMPL_HIP_SAFE_CALL(hipStreamSynchronize(impl_->stream));
#else
  (void)x;
  (void)d;
  (void)dl;
  (void)du;
#endif
}

template class TridiagonalCuSparse<double>;
template class TridiagonalCuSparse<float>;

} // namespace solver
} // namespace alps
