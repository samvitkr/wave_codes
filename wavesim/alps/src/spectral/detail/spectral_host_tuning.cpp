//
// Created by xuananqing on 3/22/23.
//

#include <spectral/spectral_host.h>

#include <common/base/logging.h>
#include <common/kokkos_abstraction/pool_space.h>
#include <common/utils/bench.h>

#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE
#include <OpenMP/Kokkos_OpenMP.hpp>
#undef KOKKOS_IMPL_PUBLIC_INCLUDE
#endif

namespace alps::spectral {

template<class DataType, class... Props>
using buffer_t = Kokkos::View<DataType,
                              Kokkos::LayoutLeft,
                              ::alps::memory_pool<Kokkos::HostSpace>,
                              Props...>;

template<class T>
void SpectralPlan<T, Kokkos::OpenMP, fft::FFTW>::
  tune_transpose_and_transform_chunk_size()
{
  // For transpose across multiple processes using MPI, small chunk sizes do not
  // seem to benefit so the tuning is skipped
  if (pencil.comm().dims[1] > 1) {
    transpose_and_transform_chunk_size_z = 0;
    return;
  }

  auto logger = get_logger("transform_bench");

  double best_time       = std::numeric_limits<double>::max();
  int    best_chunk_size = 0;
  for (int chunk_size = pencil.extent(2); chunk_size > 0; --chunk_size) {
    transpose_and_transform_chunk_size_z = chunk_size;

    buffer_t<T***> input("", create_local_layout(pencil));
    buffer_t<T***> output("", input.layout());

    ::alps::bench::Bench b;
    b.output(nullptr).performanceCounters(false);
    // Calculate an in-place ddy
    b.run("ddy", [&output, &input, space = Kokkos::OpenMP(), this] {
      do_ddy(output, input, SpectralPostOp::AssignAfterTranspose, space);
      space.fence();
    });
    auto elapsed =
      b.results().back().median(::alps::bench::Result::Measure::elapsed);
    if (elapsed < best_time) {
      best_time       = elapsed;
      best_chunk_size = chunk_size;
    }

    logger->trace("Bench y transform of size {}x{}x{} with chunk size {}: t={}",
                  pencil.extent(0),
                  pencil.extent(1),
                  pencil.extent(2),
                  chunk_size,
                  elapsed);
  }

  transpose_and_transform_chunk_size_z = best_chunk_size;

  logger->debug(
    "Tuning result for y transform of size {}x{}x{}: best chunk size {}",
    pencil.extent(0),
    pencil.extent(1),
    pencil.extent(2),
    best_chunk_size);
}

template void SpectralPlan<float, Kokkos::OpenMP, fft::FFTW>::
  tune_transpose_and_transform_chunk_size();
template void SpectralPlan<double, Kokkos::OpenMP, fft::FFTW>::
  tune_transpose_and_transform_chunk_size();
} // namespace alps::spectral
