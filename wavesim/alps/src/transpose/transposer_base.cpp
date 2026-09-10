//
// Created by xuananqing on 3/27/23.
//

#include "transposer_base.h"

#include "transposer_mpi_all2all.h"
#include "transposer_mpi_p2p.h"
#include "transposer_mpi_p2pshm.h"
#include "transposer_single.h"
#include <common/base/logging.h>
#include <common/device/device_traits.h>
#include <common/kokkos_abstraction/pool_space.h>
#include <common/utils/bench.h>

#include <Kokkos_Core.hpp>
#include <fmt/format.h>
#include <mpipp/collectives.h>

#include <limits>
#include <memory>
#include <vector>

namespace alps::transpose {

namespace {
template<class T, typename ExecSpace>
std::unique_ptr<TransposerBase<T, ExecSpace>>
create_transposer_impl(const mpipp::communicator& comm,
                       int                        n0,
                       int                        n1,
                       int                        max_nz,
                       TransposerOptions          options)
{
  using base_t = TransposerBase<T, ExecSpace>;

  if (options.method == TransposeMethod::Single) {
    return std::unique_ptr<base_t>(static_cast<base_t*>(
      new TransposerSingle<T, ExecSpace>(n0, n1, max_nz, std::move(options))));
  }

  if (options.method == TransposeMethod::All2All) {
    return std::unique_ptr<base_t>(
      static_cast<base_t*>(new TransposerMPIAll2All<T, ExecSpace>(
        comm, n0, n1, max_nz, std::move(options))));
  }

  if (options.method == TransposeMethod::Point2Point) {
    return std::unique_ptr<base_t>(
      static_cast<base_t*>(new TransposerMPIPoint2Point<T, ExecSpace>(
        comm, n0, n1, max_nz, std::move(options))));
  }

  if (options.method == TransposeMethod::Point2PointSHM) {
    if constexpr (is_openmp_execution_space_v<ExecSpace>) {
      return std::unique_ptr<base_t>(
        static_cast<base_t*>(new TransposerMPIPoint2PointSHM<T, Kokkos::OpenMP>(
          comm, n0, n1, max_nz, std::move(options))));
    } else {
      throw std::invalid_argument(
        "Point2PointSHM requires OpenMP execution space");
    }
  }

  throw std::runtime_error("Unknown or unimplemented method");
}

template<typename T, typename ExecSpace>
TransposerOptions tune_transpose_method(TransposerOptions const      options,
                                        const mpipp::communicator&   comm,
                                        int                          n0,
                                        int                          n1,
                                        int                          nz,
                                        std::vector<TransposeMethod> methods,
                                        Logger                       logger)
{
  constexpr double MAX_BENCH_TIME       = 5.0; // seconds
  constexpr double DEFAULT_BENCH_EPOCHS = 19.0;
  auto             tuned_options        = options;
  auto             best_time            = std::numeric_limits<double>::max();

  auto const exec = ExecSpace();
  auto const rank = comm.rank();
  if (rank == 0) {
    logger->debug(
      "Benchmarking and tuning transposition methods on {}x{}x{}", n0, n1, nz);
  }

  using view_t = Kokkos::
    View<T***, Kokkos::LayoutLeft, PoolSpace<typename ExecSpace::memory_space>>;
  view_t in("input", n0, n1 / comm.size(), nz);
  view_t out("output", n1, n0 / comm.size(), nz);
  Kokkos::deep_copy(exec, in, 1.0);
  Kokkos::deep_copy(exec, out, 2.0);
  exec.fence();

  for (auto method : methods) {
    if (rank == 0) {
      logger->debug("Benchmarking method: {}",
                    TransposeMethod_traits::to_string_or_throw(method));
    }

    auto opts   = options;
    opts.method = method;
    auto transposer =
      create_transposer_impl<T, ExecSpace>(comm, n0, n1, nz, opts);

    namespace bench = ::alps::bench;
    bench::Bench b;
    b.comm(comm.raw_handle()).output(nullptr).performanceCounters(false);

    // Perform a pilot run to estimate the time per measurement
    exec.fence();
    b.epochs(1).run("pilot", [&] {
      transposer->execute(out, in, nz, TransposeOpAssign(), exec);
      exec.fence();
    });
    auto epochs = [&] {
      auto local = b.results().back().sumProduct(
        bench::Result::Measure::elapsed, bench::Result::Measure::iterations);
      auto global{local};
      mpipp::allreduce(local, global, mpipp::max<decltype(local)>(), comm);
      return static_cast<size_t>(std::max(
        1.0, std::min(DEFAULT_BENCH_EPOCHS, MAX_BENCH_TIME / global + 0.5)));
    }();

    // Perform the real benchmark
    exec.fence();
    b.epochs(epochs).run("Transpose", [&] {
      transposer->execute(out, in, nz, TransposeOpAssign(), exec);
      exec.fence();
    });
    auto elapsed = [&] {
      auto local = b.results().back().median(bench::Result::Measure::elapsed);
      auto global{local};
      mpipp::allreduce(local, global, mpipp::max<decltype(local)>(), comm);
      return global;
    }();

    if (rank == 0) {
      logger->debug("Benchmarking method: {} took t={}s",
                    TransposeMethod_traits::to_string_or_throw(method),
                    elapsed);
    }

    if (elapsed < best_time) {
      best_time            = elapsed;
      tuned_options.method = method;
    }
  }

  if (rank == 0) {
    logger->debug(
      "Best transpose method: {} (t={})",
      TransposeMethod_traits::to_string_or_throw(tuned_options.method),
      best_time);
  }

  return tuned_options;
}
} // anonymous namespace

template<class T, typename ExecSpace>
std::unique_ptr<TransposerBase<T, ExecSpace>>
create_transposer(const mpipp::communicator& comm,
                  int                        n0,
                  int                        n1,
                  int                        max_nz,
                  TransposerOptions          options)
{
  if (options.method == TransposeMethod::Default) {
    if (comm.size() == 1) {
      options.method = TransposeMethod::Single;
    } else {
      options.method = TransposeMethod::All2All;
    }
  }

  if (options.method == TransposeMethod::Autotune) {
    if (comm.size() == 1) {
      options.method = TransposeMethod::Single;
    } else {
      auto candidates = std::vector<TransposeMethod>{
        TransposeMethod::All2All, TransposeMethod::Point2Point};
      if constexpr (is_openmp_execution_space_v<ExecSpace>) {
        candidates.push_back(TransposeMethod::Point2PointSHM);
      }
      options = tune_transpose_method<T, ExecSpace>(
        options, comm, n0, n1, max_nz, candidates, get_logger("transposer"));
    }
  }

  return create_transposer_impl<T, ExecSpace>(comm, n0, n1, max_nz, options);
}

#define INSTANTIATE(T, ExecSpace)                        \
  template std::unique_ptr<TransposerBase<T, ExecSpace>> \
  create_transposer<T, ExecSpace>(                       \
    const mpipp::communicator&, int, int, int, TransposerOptions);

INSTANTIATE(double, Kokkos::DefaultExecutionSpace)
INSTANTIATE(float, Kokkos::DefaultExecutionSpace)

#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
INSTANTIATE(double, Kokkos::DefaultHostExecutionSpace)
INSTANTIATE(float, Kokkos::DefaultHostExecutionSpace)
#endif

#undef INSTANTIATE

template<class T, class ExecSpace>
TransposerBase<T, ExecSpace>::TransposerBase(int n0, int n1, int n_blocks)
  : n0np{n0 / n_blocks}
  , n1np{n1 / n_blocks}
  , np{n_blocks}
  , logger_{get_logger("transposer")}
{
  if (np <= 0) {
    throw std::invalid_argument("Number of process blocks must be positive");
  }
  if (n0 <= 0 || n1 <= 0) {
    throw std::invalid_argument("Grid dimensions must be positive");
  }
  if (n0 % np != 0 || n1 % np != 0) {
    auto err_msg = fmt::format("Number of process blocks {} does not "
                               "divide the grid number {} or {}",
                               np,
                               n0,
                               n1);
    throw std::invalid_argument(err_msg);
  }
}

template<class T, class ExecSpace>
TransposerBase<T, ExecSpace>::~TransposerBase() = default;

template<class T, class ExecSpace>
void TransposerBase<T, ExecSpace>::execute(OutType const&   out,
                                           InType const&    in,
                                           TransposeOps     op,
                                           ExecSpace const& space) const
{
  const auto howmany = Kokkos::min(in.extent_int(2), out.extent_int(2));
  this->execute(out, in, howmany, op, space);
}

template<class T, class ExecSpace>
void TransposerBase<T, ExecSpace>::execute(OutType const&   out,
                                           InType const&    in,
                                           int              howmany,
                                           TransposeOps     op,
                                           ExecSpace const& space) const
{
  if (howmany < 0 || howmany > in.extent_int(2)
      || howmany > out.extent_int(2)) {
    auto err_msg = fmt::format("Invalid number of planes to transpose: {}. "
                               "Must be non-negative and no more than the "
                               "third dimension of input ({}) and output ({})",
                               howmany,
                               in.extent_int(2),
                               out.extent_int(2));
    throw std::invalid_argument(err_msg);
  }
  if (in.extent(0) < (size_t)n0np * np || in.extent(1) < (size_t)n1np) {
    auto err_msg = fmt::format("Input view is too small for the configured "
                               "grid and process blocks. Expected at least "
                               "{}x{} but got {}x{}",
                               n0np * np,
                               n1np,
                               in.extent(0),
                               in.extent(1));
    throw std::invalid_argument(err_msg);
  }
  if (out.extent(0) < (size_t)n1np * np || out.extent(1) < (size_t)n0np) {
    auto err_msg = fmt::format("Output view is too small for the configured "
                               "grid and process blocks. Expected at least "
                               "{}x{} but got {}x{}",
                               n1np * np,
                               n0np,
                               out.extent(0),
                               out.extent(1));
    throw std::invalid_argument(err_msg);
  }

  logger_->trace("Transpose {} ({}x{}x{}) to {} ({}x{}x{}) on {}",
                 in.label(),
                 in.extent(0),
                 in.extent(1),
                 howmany,
                 out.label(),
                 out.extent(0),
                 out.extent(1),
                 howmany,
                 space.name());

  if (howmany == 0) {
    logger_->trace("No planes to transpose, skipping execution");
    return;
  }
  this->execute_impl(out, in, howmany, op, space);
}

template class TransposerBase<float, Kokkos::DefaultExecutionSpace>;
template class TransposerBase<double, Kokkos::DefaultExecutionSpace>;

#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
template class TransposerBase<float, Kokkos::DefaultHostExecutionSpace>;
template class TransposerBase<double, Kokkos::DefaultHostExecutionSpace>;
#endif

} // namespace alps::transpose
