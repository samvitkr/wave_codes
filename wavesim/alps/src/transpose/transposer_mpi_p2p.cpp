#include "transposer_mpi_p2p.h"

#include <common/async/event.h>
#include <common/base/logging.h>
#include <common/base/macros.h>
#include <common/device/device_traits.h>
#include <common/kokkos_abstraction/exec_policy.h>
#include <common/kokkos_abstraction/pool_space.h>

#include <Kokkos_Core.hpp>
#include <Kokkos_Profiling_ScopedRegion.hpp>
#include <mpipp/point2point.h>

#include <stdexcept>

namespace alps::transpose {

template<class DataType, class... Properties>
using MDView = Kokkos::View<DataType, Kokkos::LayoutLeft, Properties...>;

template<class T, typename ExecSpace>
TransposerMPIPoint2Point<T, ExecSpace>::TransposerMPIPoint2Point(
  mpipp::communicator comm,
  int                 n0,
  int                 n1,
  int /*max_nz_hint*/,
  TransposerOptions /*options*/)
  : base_t(n0, n1, comm.size())
  , comm_(comm)
  , pack_streams_{2}
  , event_pool_()
{
  event_pool_.resize(comm.size() + 2);
}

template<class T, typename ExecSpace>
TransposerMPIPoint2Point<T, ExecSpace>::~TransposerMPIPoint2Point() = default;

template class TransposerMPIPoint2Point<float, Kokkos::OpenMP>;
template class TransposerMPIPoint2Point<double, Kokkos::OpenMP>;

#if defined(KOKKOS_ENABLE_CUDA)
template class TransposerMPIPoint2Point<float, Kokkos::Cuda>;
template class TransposerMPIPoint2Point<double, Kokkos::Cuda>;
#elif defined(KOKKOS_ENABLE_HIP)
template class TransposerMPIPoint2Point<float, Kokkos::HIP>;
template class TransposerMPIPoint2Point<double, Kokkos::HIP>;
#endif

namespace {

template<class InType>
struct pack_functor
{
  using T         = typename InType::non_const_value_type;
  using ExecSpace = typename InType::execution_space;

  typename InType::const_type src;
  T* KOKKOS_RESTRICT          ptr_dst;
  int                         n0pp, n1pp, offset, howmany;

  pack_functor(T*            output,
               InType const& input,
               int           block_i,
               int           n0_per_proc,
               int           n1_per_proc,
               int           batch)
    : src(input)
    , ptr_dst{output}
    , n0pp{n0_per_proc}
    , n1pp{n1_per_proc}
    , offset{n0_per_proc * block_i}
    , howmany{batch}
  {}

  KOKKOS_FUNCTION void operator()(int i, int j, int k) const
  {
    ptr_dst[i + (j + k * n1pp) * n0pp] = src(i + offset, j, k);
  }

  void execute(ExecSpace const& stream) const
  {
    using tile_t        = Kokkos::Array<int, 3>;
    auto constexpr tile = is_host_execution_space_v<ExecSpace>
                          ? tile_t{128, 8, 4}
                          : tile_t{32, 4, 2};

    auto policy = LoopPolicy<3, ExecSpace>(
      stream, tile_t{0, 0, 0}, tile_t{n0pp, n1pp, howmany}, tile);
    Kokkos::parallel_for("pack", policy, *this);
  }
};

template<typename OutType, typename Op>
struct unpack_functor
{
  using T         = typename OutType::non_const_value_type;
  using ExecSpace = typename OutType::execution_space;

  OutType                  out;
  T const* KOKKOS_RESTRICT ptr_buffer;
  int                      n0pp, n1pp, howmany, offset;

  unpack_functor(OutType  output,
                 T const* input,
                 int      block_i,
                 int      n0_per_proc,
                 int      n1_per_proc,
                 int      batch,
                 Op /*op*/)
    : out(std::move(output))
    , ptr_buffer(input)
    , n0pp{n0_per_proc}
    , n1pp{n1_per_proc}
    , howmany{batch}
    , offset{n1_per_proc * block_i}
  {}

  KOKKOS_INLINE_FUNCTION void operator()(int i, int j, int k) const
  {
    Op::apply(out(i + offset, j, k), ptr_buffer[j + (i + k * n1pp) * n0pp]);
  }

  void execute(ExecSpace const& stream) const
  {
    using tile_t        = Kokkos::Array<int, 3>;
    auto constexpr tile = is_host_execution_space_v<ExecSpace>
                          ? tile_t{64, 8, 2}
                          : tile_t{16, 8, 2};
    auto policy         = LoopPolicy<3, ExecSpace>(
      stream, tile_t{0, 0, 0}, tile_t{n1pp, n0pp, howmany}, tile);

    Kokkos::parallel_for("unpack", policy, *this);
  }
};

std::pair<int, int> get_next_peers(int phase, int rank, int nproc)
{
  // mix the first half and the second half of phases
  auto s = phase / 2;
  if (phase % 2 != 0) {
    s += (nproc + 1) / 2;
  }

  // pairwise exchange
  return {(rank + s) % nproc, (rank - s + nproc) % nproc};
}
} // anonymous namespace

template<class T, class ExecSpace>
template<typename Op>
void TransposerMPIPoint2Point<T, ExecSpace>::execute_impl(
  const OutType&   out,
  const InType&    in,
  int              howmany,
  const ExecSpace& space) const
{
  using Kokkos::ALL;
  using dev_space  = memory_pool<typename ExecSpace::memory_space>;
  using host_space = memory_pool<Kokkos::SharedHostPinnedSpace>;
  using event_t    = typename async::event_pool<ExecSpace>::event_t;

  constexpr int  tag                = 64 << 8;
  constexpr auto mpi_can_access_dev = mpi_can_access_v<dev_space>;

  auto scope = Kokkos::Profiling::ScopedRegion("TransposeMPI p2p");

  // synchronize all activities with the passed in space
  auto blocking_event = event_pool_.acquire();
  enqueue(blocking_event, space);
  for (auto i = 1; i <= pack_streams_.actual_size(); ++i) {
    wait_for(blocking_event, pack_streams_.get_stream(i));
  }

  // reserve nproc send and recv requests respectively
  const auto           nproc = comm_.size();
  const auto           rank  = comm_.rank();
  mpipp::irequest_pool send_reqs;
  mpipp::irequest_pool all_reqs;
  send_reqs.reserve(nproc);
  all_reqs.reserve(nproc * 2);

  const auto blk_size    = n0np * n1np * howmany;
  const auto alloc_bytes = blk_size * nproc * sizeof(T);

  auto send_buf = MDView<T**, dev_space>(
    (T*)dev_space().allocate(alloc_bytes), blk_size, nproc);
  auto recv_buf = MDView<T**, dev_space>(
    (T*)dev_space().allocate(alloc_bytes), blk_size, nproc);
  auto send_buf_host = MDView<T**, host_space>();
  auto recv_buf_host = MDView<T**, host_space>();
  if constexpr (!mpi_can_access_dev) {
    send_buf_host = Kokkos::create_mirror(
      Kokkos::WithoutInitializing, host_space(), send_buf);
    recv_buf_host = Kokkos::create_mirror(
      Kokkos::WithoutInitializing, host_space(), recv_buf);
  }

  auto mpi_recv_block = [=, &recv_buf, &recv_buf_host](int m) {
    auto  stride = recv_buf.stride(1);
    auto* base   = mpi_can_access_dev ? recv_buf.data() : recv_buf_host.data();
    return nonstd::span(base + m * stride, blk_size);
  };

  std::vector<std::pair<int, event_t>> pack_events;
  // lambda function to test event completion and send out data
  auto send_completed = [&, blk_size, comm = this->comm_](
                          mpipp::irequest_pool& reqs, std::size_t sent) {
    for (; sent < pack_events.size(); ++sent) {
      auto& [dst, evt] = pack_events.at(sent);
      if (!is_complete(evt)) break;

      auto* ptr = mpi_can_access_dev ? send_buf.data() : send_buf_host.data();
      ptr += dst * send_buf.stride(1);
      reqs.push(mpipp::isend(nonstd::span(ptr, blk_size), dst, tag, comm));
    }
    return sent;
  };

  // Pack data one by one
  auto        self_event{event_pool_.acquire()};
  std::size_t send_count{0};
  int         remaining_recvs{0};
  for (int s = 0; s < nproc; ++s) {
    auto [dst, src]    = get_next_peers(s, rank, nproc);
    auto const  stream = pack_streams_.get_next_stream();
    auto* const ptr    = send_buf.data() + dst * send_buf.stride(1);
    pack_functor(ptr, in, dst, n0np, n1np, howmany).execute(stream);
    if (dst != rank && !mpi_can_access_dev) {
      // for non-self send, need to copy the packed data to host
      Kokkos::deep_copy(
        stream, subview(send_buf_host, ALL, dst), subview(send_buf, ALL, dst));
    }
    // record a event after each packing and possible d2h copy
    if (dst == rank) {
      enqueue(self_event, stream);
    } else {
      enqueue(pack_events.emplace_back(dst, event_pool_.acquire()).second,
              stream);
    }

    // post non-blocking receive
    if (src != rank) {
      all_reqs.push(mpipp::irecv(mpi_recv_block(src), src, tag, comm_));
      ++remaining_recvs;
    }

    // Try to send out previously packed data
    send_count = send_completed(send_reqs, send_count);
  }
  // Attach send requests to recv_reqs to simplify wait logic
  auto recv_reqs_size = (int)all_reqs.size();
  all_reqs.push(std::move(send_reqs));

  // unpack data from self
  auto stream0 = pack_streams_.get_next_stream();
  wait_for(self_event, stream0);
  auto* ptr0 = send_buf.data() + rank * send_buf.stride(1);
  unpack_functor(out, ptr0, rank, n0np, n1np, howmany, Op{}).execute(stream0);

  // unpack data from other processes
  for (auto finished_req = all_reqs.testsome();
       remaining_recvs > 0 || send_count < pack_events.size();
       finished_req = all_reqs.testsome()) {
    send_count = send_completed(all_reqs, send_count);

    for (auto i_finished : finished_req) {
      if (i_finished >= recv_reqs_size) continue; // skip send requests
      const int  m      = all_reqs.get_status(i_finished).source();
      const auto stream = pack_streams_.get_next_stream();

      if constexpr (!mpi_can_access_dev) {
        Kokkos::deep_copy(
          stream, subview(recv_buf, ALL, m), subview(recv_buf_host, ALL, m));
      }
      auto* ptr_m = recv_buf.data() + m * recv_buf.stride(1);
      unpack_functor(out, ptr_m, m, n0np, n1np, howmany, Op{}).execute(stream);
      --remaining_recvs;
    }
  }
  all_reqs.waitall();

  pack_streams_.fence(); // ensure all operations are done before deallocation
  dev_space().deallocate(recv_buf.data(), alloc_bytes);
  dev_space().deallocate(send_buf.data(), alloc_bytes);
}

template<class T, class ExecSpace>
void TransposerMPIPoint2Point<T, ExecSpace>::execute_impl(
  const OutType&   out,
  const InType&    in,
  int              howmany,
  TransposeOps     op,
  const ExecSpace& space) const
{
  std::visit(
    [&, this](auto&& op_) {
      using Op = std::decay_t<decltype(op_)>;
      this->execute_impl<Op>(out, in, howmany, space);
    },
    op);
}

} // namespace alps::transpose
