#include "transposer_single.h"

#include <common/base/logging.h>
#include <common/kokkos_abstraction/pool_space.h>
#include <common/utils/bench.h>

#include <Kokkos_Core.hpp>

#include <array>
#include <limits>
#include <mutex>

namespace alps::transpose {
namespace {
template<class T>
std::array<int, 2> tune_tile_sizes(int dim0, int dim1, int dim2);

/// @brief Functor to transpose a 3D array.
template<typename Op, typename OutType, typename InType>
void transpose_single_openmp_impl(OutType const&        output,
                                  InType const&         input,
                                  std::array<int, 3>    i_extents,
                                  std::array<int, 2>    tile_dims,
                                  Kokkos::OpenMP const& space)
{
  auto* KOKKOS_RESTRICT       out       = output.data();
  auto const* KOKKOS_RESTRICT in        = input.data();
  std::array<size_t, 2> const o_strides = {output.stride(1), output.stride(2)};
  std::array<size_t, 2> const i_strides = {input.stride(1), input.stride(2)};

  auto* m_instance = space.impl_internal_space_instance();

  // serialize the kernel launch
  std::lock_guard<std::mutex> lock(m_instance->m_instance_mutex);
#pragma omp parallel for schedule(static) collapse(3) \
  num_threads(m_instance->thread_pool_size())
  for (int k = 0; k < i_extents[2]; ++k) {
    for (int jj = 0; jj < i_extents[0]; jj += tile_dims[1]) {
      for (int ii = 0; ii < i_extents[1]; ii += tile_dims[0]) {
        const int j_max = std::min(jj + tile_dims[1], i_extents[0]);
        const int i_max = std::min(ii + tile_dims[0], i_extents[1]);
        for (int j = jj; j < j_max; ++j) {
#pragma omp simd
          for (int i = ii; i < i_max; ++i) {
            Op::apply(out[i + j * o_strides[0] + k * o_strides[1]],
                      in[j + i * i_strides[0] + k * i_strides[1]]);
          }
        }
      }
    }
  }
}
} // anonymous namespace

template<class T, typename ExecSpace>
TransposerSingle<T, ExecSpace>::TransposerSingle(
  int                     n0,
  int                     n1,
  int                     max_nz_hint,
  const TransposerOptions options)
  : base_t(n0, n1, 1)
  , tile_sizes_{options.tile_sizes}
{
  static_assert(std::is_same_v<ExecSpace, Kokkos::OpenMP>);
  if (options.tune_tile_sizes) {
    logger_->trace("Tune TransposerSingle (from {0}x{1}x{2} to "
                   "{1}x{0}x{2}) on {3}",
                   n0,
                   n1,
                   max_nz_hint,
                   ExecSpace().name());

    tile_sizes_ = tune_tile_sizes<T>(n0, n1, max_nz_hint);

    logger_->debug("Best tile size tuned for TransposeSingle from "
                   "({0}x{1}x{2}) to ({1}x{0}x{2}) on {3}: {4} x {5}",
                   n0,
                   n1,
                   max_nz_hint,
                   ExecSpace().name(),
                   tile_sizes_[0],
                   tile_sizes_[1]);
  }
  if (tile_sizes_[0] <= 0 || tile_sizes_[1] <= 0) {
    tile_sizes_ = std::array<int, 2>{128, 128};
    logger_->warn("TransposerSingle tile sizes are not set. Using default "
                  "tile size: {} x {}.",
                  tile_sizes_[0],
                  tile_sizes_[1]);
  }
}

template<class T, class ExecSpace>
void TransposerSingle<T, ExecSpace>::execute_impl(const OutType&   out,
                                                  const InType&    in,
                                                  int              howmany,
                                                  TransposeOps     op,
                                                  const ExecSpace& space) const
{
  static_assert(std::is_same_v<ExecSpace, Kokkos::OpenMP>);
  Kokkos::Tools::pushRegion("TransposeSingle");
  std::visit(
    [&, this](auto&& op_) {
      using Op = std::decay_t<decltype(op_)>;
      transpose_single_openmp_impl<Op>(
        out, in, {n0np, n1np, howmany}, tile_sizes_, space);
    },
    op);
  Kokkos::Tools::popRegion();
}

// Explicit instantiation
template class TransposerSingle<float, Kokkos::OpenMP>;
template class TransposerSingle<double, Kokkos::OpenMP>;

namespace {

template<class T, class MemSpace>
void omp_fill_zero(Kokkos::View<T***, Kokkos::LayoutLeft, MemSpace> const& view)
{
#pragma omp parallel for collapse(2)
  for (int k = 0; k < view.extent_int(2); ++k) {
    for (int j = 0; j < view.extent_int(1); ++j) {
#pragma omp simd
      for (int i = 0; i < view.extent_int(0); ++i) {
        view(i, j, k) = 0;
      }
    }
  }
}

template<class T>
std::array<int, 2> tune_tile_sizes(int dim0, int dim1, int dim2)
{
  using mem_space = ::alps::memory_pool<Kokkos::HostSpace>;
  auto logger     = get_logger("transposer_bench");

  auto       exec = Kokkos::OpenMP();
  const auto max_threads_per_block =
    ::Kokkos::Impl::get_tile_size_properties(exec).max_total_tile_size;
  constexpr auto min_threads_per_block = 32;
  constexpr auto max_tile_dim          = 512;
  constexpr auto min_tile_dim          = 4;

  std::array<int, 2> best_tile{};
  double             best_time = std::numeric_limits<double>::max();

  // Allocate and initialize the workspace
  Kokkos::View<T***, Kokkos::LayoutLeft, mem_space> src(
    Kokkos::view_alloc("", Kokkos::WithoutInitializing), dim0, dim1, dim2);
  Kokkos::View<T***, Kokkos::LayoutLeft, mem_space> dst(
    Kokkos::view_alloc("", Kokkos::WithoutInitializing), dim1, dim0, dim2);
  // Fill with non-zero values to avoid any optimization
  Kokkos::deep_copy(src, T(1));
  Kokkos::deep_copy(dst, T(2));

  // Loop thru all tile dimensions with a constraint on the total tile threads
  for (int tile_1 = min_tile_dim; tile_1 <= max_tile_dim; tile_1 *= 2) {
    for (int tile_0 = min_tile_dim; tile_0 <= max_tile_dim; tile_0 *= 2) {
      if (tile_0 * tile_1 > max_threads_per_block) break;
      if (tile_0 * tile_1 < min_threads_per_block) continue;

      logger->trace(
        "Bench transpose {}x{} with tile {}x{}", dim0, dim1, tile_0, tile_1);

      ::alps::bench::Bench b;
      b.output(nullptr).performanceCounters(false);
      b.run("TransposeSingle", [&] {
        transpose_single_openmp_impl<TransposeOpAssign>(
          dst, src, {dim0, dim1, dim2}, {tile_0, tile_1}, exec);
        exec.fence();
      });
      auto result  = b.results().front();
      auto elapsed = result.average(::alps::bench::Result::Measure::elapsed);

      logger->trace("Bench result {}x{} with tile {}x{}: {:.2f} ms",
                    dim0,
                    dim1,
                    tile_0,
                    tile_1,
                    elapsed * 1000);
      if (elapsed < best_time) {
        best_time = elapsed;
        best_tile = std::array<int, 2>{tile_0, tile_1};
      }
    }
  }

  return best_tile;
}
} // anonymous namespace
} // namespace alps::transpose
