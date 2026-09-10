#include "transposer_single.h"

#include <common/base/logging.h>
#include <common/device/device_traits.h>

#include <Kokkos_Core.hpp>

namespace alps::transpose {
namespace {

#if defined(KOKKOS_ENABLE_CUDA)
constexpr int TILE_DIM   = 32;
constexpr int BLOCK_ROWS = 16;
#elif defined(KOKKOS_ENABLE_HIP)
constexpr int TILE_DIM   = 64;
constexpr int BLOCK_ROWS = 8;
#endif
/// @brief Functor to transpose a 3D array.
template<typename T, typename Op, bool isMultiples>
__global__ void transpose_cuda_kernel(T* __restrict__ odata,
                                      T const* __restrict__ idata,
                                      unsigned dim0,
                                      unsigned dim1,
                                      unsigned istride1,
                                      unsigned istride2,
                                      unsigned ostride1,
                                      unsigned ostride2)
{
  __shared__ T tile[TILE_DIM][TILE_DIM + 1];

  unsigned int x        = blockIdx.x * TILE_DIM + threadIdx.x;
  unsigned int y        = blockIdx.y * TILE_DIM + threadIdx.y;
  auto         index_in = x + y * istride1 + blockIdx.z * istride2;

  for (int j = 0; j < TILE_DIM; j += BLOCK_ROWS) {
    if (isMultiples || (x < dim0 && (y + j) < dim1)) {
      tile[threadIdx.y + j][threadIdx.x] = idata[index_in + j * istride1];
    }
  }

  x              = blockIdx.y * TILE_DIM + threadIdx.x;
  y              = blockIdx.x * TILE_DIM + threadIdx.y;
  auto index_out = x + y * ostride1 + blockIdx.z * ostride2;

  __syncthreads();

  for (int j = 0; j < TILE_DIM; j += BLOCK_ROWS) {
    if (isMultiples || (x < dim1 && (y + j) < dim0)) {
      Op::apply(odata[index_out + j * ostride1],
                tile[threadIdx.x][threadIdx.y + j]);
    }
  }
}
} // anonymous namespace

template<class T, typename ExecSpace>
TransposerSingle<T, ExecSpace>::TransposerSingle(int n0,
                                                 int n1,
                                                 int /*max_nz_hint*/,
                                                 TransposerOptions /*options*/)
  : base_t(n0, n1, 1)
{
  static_assert(is_cuda_execution_space_v<ExecSpace>
                || is_hip_execution_space_v<ExecSpace>);
}

template<class T, class ExecSpace>
void TransposerSingle<T, ExecSpace>::execute_impl(const OutType&   out,
                                                  const InType&    in,
                                                  int              howmany,
                                                  TransposeOps     op,
                                                  const ExecSpace& space) const
{
  Kokkos::Tools::pushRegion("TransposeSingle");
  auto ceil_div = [](unsigned x) { return (x + TILE_DIM - 1) / TILE_DIM; };

  auto* kernel = std::visit(
    [divisible = (n0np % TILE_DIM == 0) && (n1np % TILE_DIM == 0)](auto&& op_) {
      using Op = std::decay_t<decltype(op_)>;
      return divisible ? (void*)(&transpose_cuda_kernel<T, Op, true>)
                       : (void*)(&transpose_cuda_kernel<T, Op, false>);
    },
    op);
  uint3 dims{unsigned(n0np), unsigned(n1np), unsigned(howmany)};
  dim3  grid_dim(ceil_div(dims.x), ceil_div(dims.y), dims.z);
  dim3  block_dim(TILE_DIM, BLOCK_ROWS, 1);

  auto* odata             = out.data();
  auto* idata             = in.data();
  auto  istride1          = unsigned(in.stride(1));
  auto  istride2          = unsigned(in.stride(2));
  auto  ostride1          = unsigned(out.stride(1));
  auto  ostride2          = unsigned(out.stride(2));
  void* launchArguments[] = {&odata,
                             &idata,
                             &dims.x,
                             &dims.y,
                             &istride1,
                             &istride2,
                             &ostride1,
                             &ostride2};
#if defined(KOKKOS_ENABLE_CUDA)
  KOKKOS_IMPL_CUDA_SAFE_CALL(cudaLaunchKernel(
    kernel, grid_dim, block_dim, launchArguments, 0, space.cuda_stream()));
#elif defined(KOKKOS_ENABLE_HIP)
  KOKKOS_IMPL_HIP_SAFE_CALL(hipLaunchKernel(
    kernel, grid_dim, block_dim, launchArguments, 0, space.hip_stream()));
#endif
  Kokkos::Tools::popRegion();
}

template class TransposerSingle<float, Kokkos::DefaultExecutionSpace>;
template class TransposerSingle<double, Kokkos::DefaultExecutionSpace>;

} // namespace alps::transpose
