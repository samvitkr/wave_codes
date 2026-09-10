//
// Created by xuanx004 on 3/23/23.
//

#include "../tridiagonal_wang.h"

#include <common/base/logging.h>
#include <common/kokkos_abstraction/pool_space.h>

#include <Kokkos_Core.hpp>
#include <Kokkos_Timer.hpp>
#include <mpipp/collectives.h>

namespace alps::solver {

namespace {
/// Divide n into n_chunks parts, and return the index range of each part.
/** First n_chunks-1 parts have the same size and the last part may be smaller.
 */
std::vector<std::pair<int, int>> distribute_work(int n, int n_chunks)
{
  int const                        chunk_size = (n + n_chunks - 1) / n_chunks;
  std::vector<std::pair<int, int>> distribution;
  if ((n_chunks - 1) * chunk_size >= n) return distribution;
  distribution.reserve(n_chunks);

  for (int i = 0; i < n_chunks; ++i) {
    distribution.emplace_back(i * chunk_size,
                              std::min(n, (i + 1) * chunk_size));
  }
  return distribution;
}
} // anonymous namespace

template<typename ValueT>
void TridiagonalWang<ValueT>::tune_algorithms(coeff_t const& d,
                                              coeff_t const& dl,
                                              coeff_t const& du)
{
  // Tune the number of chunks when doing the reduce_and_backward step
  double          best_time  = std::numeric_limits<double>::max();
  int             best_chunk = 1;
  ReductionMethod best_method{ReductionMethod::MPI_Alltoall};

  const auto rank = rank_;
  if (rank == 0) {
    this->logger_->debug(
      "Wang reduce_and_backward: Benchmarking and tuning algorithms...");
  }
  Kokkos::View<ValueT***, Kokkos::LayoutLeft, default_memory_pool> x(
    Kokkos::view_alloc("", Kokkos::WithoutInitializing), d.layout());
  static_assert(
    std::is_same_v<typename decltype(x)::array_layout, Kokkos::LayoutLeft>);

  const auto bench_timer = [&, this](std::string const& algorithm_name) {
    solve_impl(x, d, dl, du); // warm-up
    Kokkos::fence();
    mpipp::barrier(comm_);
    int constexpr n_runs = 3;
    double elapsed_time{0};
    for (int i = 0; i < n_runs; ++i) {
      const Kokkos::Timer timer;
      solve_impl(x, d, dl, du);
      Kokkos::fence();
      mpipp::barrier(comm_);
      elapsed_time = std::min(timer.seconds(), elapsed_time);
    }
    double max_used_time{elapsed_time};
    mpipp::allreduce(elapsed_time, max_used_time, mpipp::max<double>(), comm_);
    if (rank == 0) {
      this->logger_->trace("Wang reduce_and_backward using {}: # of "
                           "chunks {}, t={}",
                           algorithm_name,
                           this->ny_chunk_distribution_.size(),
                           max_used_time);
    }
    return max_used_time;
  };

  for (int n_chunks =
         (this->n2_ + minimum_ny_chunk_size - 1) / minimum_ny_chunk_size;
       n_chunks > 0;
       --n_chunks) {
    ny_chunk_distribution_ = distribute_work(this->n2_, n_chunks);
    if (ny_chunk_distribution_.empty()) continue;

    reduction_method = ReductionMethod::MPI_Alltoall;
    if (auto time = bench_timer("Alltoall"); time < best_time) {
      best_time   = time;
      best_chunk  = n_chunks;
      best_method = ReductionMethod::MPI_Alltoall;
    }
  }

  reduction_method       = best_method;
  ny_chunk_distribution_ = distribute_work(this->n2_, best_chunk);
  if (rank == 0) {
    this->logger_->debug(
      "Wang reduce_and_backward: best # of chunks {} using method {}",
      best_chunk,
      int(best_method));
  }
}

template void TridiagonalWang<double>::tune_algorithms(coeff_t const& d,
                                                       coeff_t const& dl,
                                                       coeff_t const& du);
template void TridiagonalWang<float>::tune_algorithms(coeff_t const& d,
                                                      coeff_t const& dl,
                                                      coeff_t const& du);

} // namespace alps::solver
