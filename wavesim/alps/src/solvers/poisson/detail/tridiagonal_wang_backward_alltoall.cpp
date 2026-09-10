#include "../tridiagonal_wang.h"

#include <common/async/streams.h>
#include <common/device/device_traits.h>
#include <common/kokkos_abstraction/exec_policy.h>
#include <common/kokkos_abstraction/pool_space.h>

#include <Kokkos_Core.hpp>
#include <mpipp/collectives.h>

namespace alps::solver {

template<class DataType, class... Properties>
using MDView = Kokkos::View<DataType, Kokkos::LayoutLeft, Properties...>;
using Kokkos::ALL;
using Kokkos::subview;

namespace {
struct AlltoallCountsDisplacements
{
  using counts_t = MDView<int**, Kokkos::HostSpace>;
  using displs_t = MDView<int**, Kokkos::HostSpace>;
  counts_t            send_counts, recv_counts;
  displs_t            send_displacements, recv_displacements;
  std::vector<size_t> rt_displacements;

  AlltoallCountsDisplacements(int                                     nproc,
                              int                                     rank,
                              int                                     n1,
                              std::vector<std::pair<int, int>> const& y_chunks)
    : send_counts("send counts", nproc, y_chunks.size())
    , recv_counts("recv counts", send_counts.layout())
    , send_displacements("send displs", send_counts.layout())
    , recv_displacements("recv displs", recv_counts.layout())
    , rt_displacements(y_chunks.size() + 1)
  {
    for (size_t chunk_i = 0; chunk_i < y_chunks.size(); ++chunk_i) {
      auto const& chunk   = y_chunks[chunk_i];
      auto const  y_block = chunk.second - chunk.first;

      // Calculate the send/recv counts
      for (int p = 0; p < nproc; ++p) {
        send_counts(p, chunk_i) = (int)(y_block / nproc) * n1;
      }
      if (int residual = y_block % nproc; residual != 0) {
        for (int p = 0; p < residual; ++p) {
          send_counts(p, chunk_i) += n1;
        }
      }
      for (int p = 0; p < nproc; ++p) {
        recv_counts(p, chunk_i) = send_counts(rank, chunk_i);
      }

      // Calculate the send/recv displacements
      for (int p = 1; p < nproc; ++p) {
        send_displacements(p, chunk_i) =
          send_displacements(p - 1, chunk_i) + send_counts(p - 1, chunk_i);
        recv_displacements(p, chunk_i) =
          recv_displacements(p - 1, chunk_i) + recv_counts(p - 1, chunk_i);
      }

      rt_displacements[chunk_i + 1] =
        rt_displacements[chunk_i] + recv_counts(0, chunk_i) * (nproc + 1);
    }
  }

  template<typename T, size_t U>
  void exchange_forward(nonstd::span<T, U>         sendbuf,
                        T*                         recvbuf,
                        int                        chunk_i,
                        bool                       is_uniform,
                        mpipp::communicator const& comm) const
  {
    auto const nproc = send_counts.extent(0);
    if (is_uniform) {
      mpipp::alltoall(sendbuf, send_counts(0, chunk_i), recvbuf, comm);
    } else {
      mpipp::alltoallv(sendbuf,
                       nonstd::span(&send_counts(0, chunk_i), nproc),
                       nonstd::span(&send_displacements(0, chunk_i), nproc),
                       recvbuf,
                       nonstd::span(&recv_counts(0, chunk_i), nproc),
                       nonstd::span(&recv_displacements(0, chunk_i), nproc),
                       comm);
    }
  }

  template<typename T, size_t U>
  void exchange_backward(nonstd::span<T, U>         sendbuf,
                         T*                         recvbuf,
                         int                        chunk_i,
                         bool                       is_uniform,
                         mpipp::communicator const& comm) const
  {
    auto const nproc = send_counts.extent(0);
    if (is_uniform) {
      mpipp::alltoall(sendbuf, recv_counts(0, chunk_i), recvbuf, comm);
    } else {
      mpipp::alltoallv(sendbuf,
                       nonstd::span(&recv_counts(0, chunk_i), nproc),
                       nonstd::span(&recv_displacements(0, chunk_i), nproc),
                       recvbuf,
                       nonstd::span(&send_counts(0, chunk_i), nproc),
                       nonstd::span(&send_displacements(0, chunk_i), nproc),
                       comm);
    }
  }
};

auto constexpr default_tile_size = []() -> Kokkos::Array<std::int64_t, 2> {
  using Device = Kokkos::DefaultExecutionSpace;
  if constexpr (is_cuda_execution_space_v<Device>) return {32, 4};
  if constexpr (is_hip_execution_space_v<Device>) return {64, 8};
  return {0, 0};
}();
} // namespace

template<typename ValueT>
void TridiagonalWang<ValueT>::reduce_and_backward_alltoall(
  solution_t                       x,
  coeff_t                          d,
  std::vector<std::pair<int, int>> ny_chunk_distribution,
  std::true_type /*mpi_can_access_default_space*/) const
{
  using nonstd::span;
  const auto                               n1 = this->n1_;
  MDView<ReduceX_t**, default_memory_pool> rtemp_rank(
    Kokkos::view_alloc("r_rank", Kokkos::WithoutInitializing), n1, this->n2_);

  mpipp::barrier(comm_);

  AlltoallCountsDisplacements counts_displacements(
    nproc_, rank_, xn_dev.extent_int(0), ny_chunk_distribution);
  auto const& [send_counts,
               recv_counts,
               send_displacements,
               recv_displacements,
               rt_displacements] = counts_displacements;

  const MDView<ReduceX_t*, default_memory_pool> rt_storage(
    Kokkos::view_alloc(Kokkos::WithoutInitializing, "rt"),
    rt_displacements.back());

  for (size_t chunk_i = 0; chunk_i < ny_chunk_distribution.size(); ++chunk_i) {
    auto const& chunk   = ny_chunk_distribution[chunk_i];
    auto const  y_start = chunk.first;
    auto const  y_end   = chunk.second;
    auto const  y_block = y_end - y_start;
    auto const  y_offset =
      y_start + send_displacements(rank_, chunk_i) / xn_dev.extent_int(0);
    auto const is_uniform = (y_block % nproc_ == 0);

    const MDView<ReduceX_t***> rt(
      Kokkos::view_wrap(rt_storage.data() + rt_displacements[chunk_i]),
      xn_dev.extent(0),
      recv_counts(0, chunk_i) / xn_dev.extent(0),
      nproc_ + 1);

    auto stream = get_next_stream(stream_pool_);
    auto policy2d_ =
      LoopPolicy<2>(stream, {0, 0}, {n1, rt.extent_int(1)}, default_tile_size);

    counts_displacements.exchange_forward(
      span(xn_dev.data() + y_start * xn_dev.stride(1),
           y_block * xn_dev.stride(1)),
      rt.data() + rt.stride(2),
      chunk_i,
      is_uniform,
      comm_);

    Kokkos::Tools::pushRegion("Wang solve reduced system");
    auto const  nproc = this->nproc_;
    auto const& f_g   = f_gather;
    Kokkos::parallel_for(
      "solve reduced system 1", policy2d_, KOKKOS_LAMBDA(int i, int j) {
        for (int p = 1; p < nproc; ++p) {
          rt(i, j, p + 1) =
            ReduceX_t(ValueT(rt(i, j, p + 1))
                      - f_g(i, j + y_offset, p) * ValueT(rt(i, j, p)));
        }
      });

    stream.fence();
    counts_displacements.exchange_backward(
      span(rt.data(), recv_counts(0, chunk_i) * nproc_),
      rtemp_rank.data() + rtemp_rank.stride(1) * y_start,
      chunk_i,
      is_uniform,
      comm_);

    auto const& g_g = g_gather;
    Kokkos::parallel_for(
      "solve reduced system 2", policy2d_, KOKKOS_LAMBDA(int i, int j) {
        for (int p = nproc - 2; p >= 0; --p) {
          rt(i, j, p + 1) =
            ReduceX_t(ValueT(rt(i, j, p + 1))
                      - g_g(i, j + y_offset, p + 1) * ValueT(rt(i, j, p + 2)));
        }
      });
    stream.fence();
    Kokkos::Tools::popRegion();

    counts_displacements.exchange_backward(
      span(rt.data() + rt.stride(2), recv_counts(0, chunk_i) * nproc_),
      xn_dev.data() + y_start * xn_dev.stride(1),
      chunk_i,
      is_uniform,
      comm_);

    final_substitution(x, d, chunk, rtemp_rank, stream);
  }

  stream_pool_.fence();
}

template<typename ValueT>
void TridiagonalWang<ValueT>::reduce_and_backward_alltoall(
  solution_t                       x,
  coeff_t                          d,
  std::vector<std::pair<int, int>> ny_chunk_distribution,
  std::false_type /*mpi_can_access_default_space*/) const
{
  using nonstd::span;
  const auto                               n1 = this->n1_;
  MDView<ReduceX_t**, default_memory_pool> rtemp_rank(
    Kokkos::view_alloc("r_rank", Kokkos::WithoutInitializing), n1, this->n2_);

  deep_copy(stream_pool_.get_stream(1), xn_host, xn_dev);

  mpipp::barrier(comm_);

  AlltoallCountsDisplacements counts_displacements(
    nproc_, rank_, xn_host.extent_int(0), ny_chunk_distribution);
  auto const& [send_counts,
               recv_counts,
               send_displacements,
               recv_displacements,
               rt_displacements] = counts_displacements;

  const MDView<ReduceX_t*, memory_pool<host_space>> rt_storage(
    Kokkos::view_alloc(Kokkos::WithoutInitializing, "rt"),
    rt_displacements.back());
  auto rtemp_rank_recv = Kokkos::create_mirror_view(
    Kokkos::WithoutInitializing, memory_pool<host_space>(), rtemp_rank);

  stream_pool_.get_stream(1).fence(); // Wait for deep_copy to finish

  for (size_t chunk_i = 0; chunk_i < ny_chunk_distribution.size(); ++chunk_i) {
    auto const& chunk   = ny_chunk_distribution[chunk_i];
    auto const  y_start = chunk.first;
    auto const  y_end   = chunk.second;
    auto const  y_block = y_end - y_start;
    auto const  y_offset =
      y_start + send_displacements(rank_, chunk_i) / xn_host.extent_int(0);
    auto const is_uniform = (y_block % nproc_ == 0);

    auto stream = get_next_stream(stream_pool_);

    const MDView<ReduceX_t***, host_space> rt(
      Kokkos::view_wrap(rt_storage.data() + rt_displacements[chunk_i]),
      xn_host.extent(0),
      recv_counts(0, chunk_i) / xn_host.extent(0),
      nproc_ + 1);

    counts_displacements.exchange_forward(
      span(xn_host.data() + y_start * xn_host.stride(1),
           y_block * xn_host.stride(1)),
      rt.data() + rt.stride(2),
      chunk_i,
      is_uniform,
      comm_);

    Kokkos::Tools::pushRegion("Wang solve reduced system");
#pragma omp parallel for
    for (int j = 0; j < rt.extent_int(1); ++j) {
      for (int p = 1; p < nproc_; ++p) {
#pragma omp simd
        for (int i = 0; i < n1; ++i) {
          rt(i, j, p + 1) =
            ReduceX_t(ValueT(rt(i, j, p + 1))
                      - f_gather(i, j + y_offset, p) * ValueT(rt(i, j, p)));
        }
      }
    }

    counts_displacements.exchange_backward(
      span(rt.data(), recv_counts(0, chunk_i) * nproc_),
      rtemp_rank_recv.data() + rtemp_rank_recv.stride(1) * y_start,
      chunk_i,
      is_uniform,
      comm_);

    deep_copy(stream,
              subview(rtemp_rank, ALL, chunk),
              subview(rtemp_rank_recv, ALL, chunk));

#pragma omp parallel for
    for (int j = 0; j < rt.extent_int(1); ++j) {
      for (int p = nproc_ - 2; p >= 0; --p) {
#pragma omp simd
        for (int i = 0; i < n1; ++i) {
          rt(i, j, p + 1) = ReduceX_t(ValueT(rt(i, j, p + 1))
                                      - g_gather(i, j + y_offset, p + 1)
                                          * ValueT(rt(i, j, p + 2)));
        }
      }
    }
    Kokkos::Tools::popRegion();

    counts_displacements.exchange_backward(
      span(rt.data() + rt.stride(2), recv_counts(0, chunk_i) * nproc_),
      xn_host.data() + y_start * xn_host.stride(1),
      chunk_i,
      is_uniform,
      comm_);

    deep_copy(
      stream, subview(xn_dev, ALL, chunk), subview(xn_host, ALL, chunk));

    final_substitution(x, d, chunk, rtemp_rank, stream);
  }

  stream_pool_.fence();
}

template void TridiagonalWang<double>::reduce_and_backward_alltoall(
  solution_t                       x,
  coeff_t                          d,
  std::vector<std::pair<int, int>> ny_chunk_distribution,
  std::true_type) const;
template void TridiagonalWang<float>::reduce_and_backward_alltoall(
  solution_t                       x,
  coeff_t                          d,
  std::vector<std::pair<int, int>> ny_chunk_distribution,
  std::true_type) const;
template void TridiagonalWang<double>::reduce_and_backward_alltoall(
  solution_t                       x,
  coeff_t                          d,
  std::vector<std::pair<int, int>> ny_chunk_distribution,
  std::false_type) const;
template void TridiagonalWang<float>::reduce_and_backward_alltoall(
  solution_t                       x,
  coeff_t                          d,
  std::vector<std::pair<int, int>> ny_chunk_distribution,
  std::false_type) const;

} // namespace alps::solver
