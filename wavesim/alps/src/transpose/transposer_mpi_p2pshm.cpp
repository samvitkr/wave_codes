#include "transposer_mpi_p2pshm.h"

#include <common/base/logging.h>
#include <common/device/device_traits.h>
#include <common/kokkos_abstraction/exec_policy.h>
#include <common/kokkos_abstraction/pool_space.h>

#include <Kokkos_Core.hpp>
#include <Kokkos_Profiling_ScopedRegion.hpp>
#include <fmt/format.h>
#include <mpipp/collectives.h>
#include <mpipp/comm.h>
#include <mpipp/point2point.h>
#include <mpipp/shared_window.h>
#include <nonstd/span.hpp>

#include <atomic>
#include <chrono>
#include <memory>
#include <new>
#include <stdexcept>
#include <string>
#include <thread>
#include <unistd.h>
#include <utility>

namespace alps::transpose {

// PIMPL Implementation class
template<class T, typename ExecSpace>
struct TransposerMPIPoint2PointSHM<T, ExecSpace>::Impl
{
#if __cpp_lib_hardware_interference_size >= 201603
  static constexpr std::size_t ALIGNMENT =
    std::hardware_destructive_interference_size;
#else
  static constexpr std::size_t ALIGNMENT = 64; // Fallback to 64 bytes on x86_64
#endif
  static_assert((ALIGNMENT & (ALIGNMENT - 1)) == 0,
                "Alignment must be a power of two");

  mpipp::communicator comm;
  mpipp::communicator shm_comm;
  std::vector<int>    peer_ranks_in_shm;

  mpipp::shared_window<void>          win_ctrl{MPI_WIN_NULL};
  std::vector<std::atomic<uint64_t>*> generation;
  uint64_t                            last_generation{0};

  mpipp::shared_window<void>   win{MPI_WIN_NULL};
  nonstd::span<T>              shm_buffer;
  std::vector<nonstd::span<T>> peer_buffers;

  nonstd::span<T> get_shm_buffer(std::size_t requested_size);

  static_assert(ALIGNMENT % sizeof(T) == 0,
                "Alignment must be a multiple of the data type size");

  static std::size_t aligned_bytes(std::size_t byte_size)
  {
    return (byte_size + ALIGNMENT - 1u) & -ALIGNMENT;
  }
};

template<class T, typename ExecSpace>
TransposerMPIPoint2PointSHM<T, ExecSpace>::TransposerMPIPoint2PointSHM(
  mpipp::communicator comm,
  int                 n0,
  int                 n1,
  int /*max_nz_hint*/,
  TransposerOptions /*options*/)
  : base_t(n0, n1, comm.size())
  , impl(std::make_unique<Impl>())
{
  impl->comm     = std::move(comm);
  impl->shm_comm = impl->comm.split_shared();
  std::vector<int> ranks_in_world(base_t::np, 0);
  for (int i = 0; i < base_t::np; ++i) {
    ranks_in_world.at(i) = i;
  }
  impl->peer_ranks_in_shm =
    mpipp::translate_ranks(impl->comm, ranks_in_world, impl->shm_comm);

  // Create a shared window for control variables (generations)
  // Allocate ALIGNMENT*2 bytes: worst-case std::align padding (ALIGNMENT-1)
  // plus sizeof(atomic<uint64_t>), ensuring the counter lands on its own cache
  // line to avoid false sharing with adjacent allocations.
  static_assert(std::atomic<uint64_t>::is_always_lock_free);
  static_assert(
    std::is_trivially_destructible_v<std::atomic<uint64_t>>,
    "atomic<uint64_t> must be trivially destructible for placement new in "
    "shared memory");
  impl->win_ctrl = {Impl::ALIGNMENT * 2, 1, impl->shm_comm};
  void* aligned  = [](mpipp::shared_window<void> const& win) -> void* {
    void*  base_ptr = win.base_address();
    size_t space    = win.size_bytes();
    void*  ptr      = std::align(
      Impl::ALIGNMENT, sizeof(std::atomic<uint64_t>), base_ptr, space);

    if (ptr == nullptr || space < sizeof(std::atomic<uint64_t>)) {
      throw std::runtime_error("Failed to align control shared memory buffer");
    }
    return ptr;
  }(impl->win_ctrl);
  // placement new to construct atomic in shared memory
  new (aligned) std::atomic<uint64_t>{0};

  // Gather all offsets
  std::ptrdiff_t offset =
    static_cast<std::byte*>(aligned)
    - static_cast<std::byte*>(impl->win_ctrl.base_address());
  std::vector<std::ptrdiff_t> all_offsets(impl->shm_comm.size(), 0);
  mpipp::allgather(offset, all_offsets.data(), impl->shm_comm);

  impl->generation.reserve(impl->shm_comm.size());
  for (int i = 0; i < impl->shm_comm.size(); ++i) {
    auto       peer_mem  = (impl->win_ctrl).template shared_query<std::byte>(i);
    std::byte* peer_base = peer_mem.data();
    auto*      peer_ctrl =
      reinterpret_cast<std::atomic<uint64_t>*>(peer_base + all_offsets.at(i));
    impl->generation.push_back(peer_ctrl);
  }
}

template<class T, typename ExecSpace>
TransposerMPIPoint2PointSHM<T, ExecSpace>::~TransposerMPIPoint2PointSHM() =
  default;

template class TransposerMPIPoint2PointSHM<float, Kokkos::OpenMP>;
template class TransposerMPIPoint2PointSHM<double, Kokkos::OpenMP>;

namespace {

template<typename InType>
struct pack_functor
{
  using T         = typename InType::non_const_value_type;
  using ExecSpace = typename InType::execution_space;
  using policy_t  = Kokkos::TeamPolicy<ExecSpace, Kokkos::IndexType<int>>;
  using member_t  = typename policy_t::member_type;

  typename InType::const_type in;
  T*                          out{nullptr};
  int                         n0pp, n1pp, offset, howmany;

  pack_functor(T*            output,
               InType const& input,
               int           i_block,
               int           n0_per_proc,
               int           n1_per_proc,
               int           batch)
    : in(input)
    , out(output)
    , n0pp{n0_per_proc}
    , n1pp{n1_per_proc}
    , offset{i_block * n0_per_proc}
    , howmany{batch}
  {}

  KOKKOS_INLINE_FUNCTION void operator()(const member_t& team) const
  {
    const auto m = static_cast<int>(team.league_rank());
    const auto k = m / n1pp;
    const auto j = m % n1pp;
    Kokkos::parallel_for(Kokkos::TeamVectorRange(team, n0pp), [&](int& i) {
      out[i + m * n0pp] = in(i + offset, j, k);
    });
  }

  void execute(policy_t const& policy) const
  {
    Kokkos::parallel_for("pack", policy, *this);
#pragma omp parallel
    {
      std::atomic_thread_fence(std::memory_order_release);
    }
  }
};

template<typename OutType, typename Op>
struct unpack_functor
{
  using T         = typename OutType::non_const_value_type;
  using ExecSpace = typename OutType::execution_space;
  using policy_t  = Kokkos::TeamPolicy<ExecSpace, Kokkos::IndexType<int>>;
  using member_t  = typename policy_t::member_type;

  T const* buffer;
  Kokkos::View<T***, Kokkos::LayoutLeft, typename ExecSpace::memory_space> out;
  int n0pp, n1pp, howmany, offset;

  unpack_functor(OutType const& output,
                 T const*       input,
                 int            i_block,
                 int            n0_per_proc,
                 int            n1_per_proc,
                 int            batch,
                 Op /*op*/)
    : buffer(input)
    , out(output)
    , n0pp{n0_per_proc}
    , n1pp{n1_per_proc}
    , howmany{batch}
    , offset{i_block * n1_per_proc}
  {}

  KOKKOS_INLINE_FUNCTION void operator()(const member_t& team) const
  {
    const auto  k     = static_cast<int>(team.league_rank());
    const auto* ptr_k = buffer + k * n0pp * n1pp;
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, n0pp), [&](int& j) {
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(team, n1pp), [&](int& i) {
        Op::apply(out(i + offset, j, k), ptr_k[j + i * n0pp]);
      });
    });
  }
};

template<typename T, typename MemSpace>
struct BufferDeleter
{
  std::size_t bytes;

  void operator()(T* ptr) const { MemSpace().deallocate(ptr, bytes); }
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

template<class T, typename ExecSpace>
nonstd::span<T> TransposerMPIPoint2PointSHM<T, ExecSpace>::Impl::get_shm_buffer(
  std::size_t requested_size)
{
  if (requested_size <= shm_buffer.size()) {
    return shm_buffer.subspan(0, requested_size);
  }

  if (win.is_valid()) {
    win.reset();
    shm_buffer = {};
    peer_buffers.clear();
  }

  auto const requested_bytes = requested_size * sizeof(T);
  // Total bytes needed for manual alignment
  auto const alloc_bytes = ALIGNMENT + requested_bytes;

  MPI_Info info{MPI_INFO_NULL};
  MPI_Info_create(&info);
  MPI_Info_set(info, "alloc_shared_noncontig", "true");
  win = mpipp::shared_window<void>(alloc_bytes, 1, shm_comm, info);
  MPI_Info_free(&info);

  // Use std::align to find the aligned buffer position after the header.
  shm_buffer = [=](void* const       base,
                   std::size_t const size_bytes) -> nonstd::span<T> {
    void*       base_ptr = base;
    size_t      space    = size_bytes;
    auto* const aligned =
      std::align(ALIGNMENT, requested_bytes, base_ptr, space);

    // Runtime check: ensure alignment succeeded and buffer fits
    if (aligned == nullptr || space < requested_bytes) {
      throw std::runtime_error("Failed to align shared memory buffer");
    }

    // Save the aligned buffer
    return {reinterpret_cast<T*>(aligned), requested_size};
  }(win.base_address(), win.size_bytes());

  // Ensure the buffer is touched
  auto pagesize = sysconf(_SC_PAGESIZE);
  for (size_t offset = 0; offset < shm_buffer.size_bytes();
       offset += pagesize) {
    volatile char* ptr = (volatile char*)shm_buffer.data() + offset;
    *ptr               = 0;
  }
  std::atomic_thread_fence(std::memory_order_release);

  std::ptrdiff_t offset = reinterpret_cast<std::byte*>(shm_buffer.data())
                        - reinterpret_cast<std::byte*>(win.base_address());
  std::vector<std::ptrdiff_t> all_offsets(shm_comm.size(), 0);
  mpipp::allgather(offset, all_offsets.data(), shm_comm);

  // Populate the peer buffer pointers by reading each peer's header
  peer_buffers.reserve(shm_comm.size());
  for (int i = 0; i < shm_comm.size(); ++i) {
    auto  peer_mem     = win.shared_query<std::byte>(i);
    auto* peer_base    = peer_mem.data();
    auto* aligned_peer = reinterpret_cast<T*>(peer_base + all_offsets.at(i));
    peer_buffers.emplace_back(aligned_peer, requested_size);
  }

  return shm_buffer;
}

template<class T, class ExecSpace>
template<typename Op>
void TransposerMPIPoint2PointSHM<T, ExecSpace>::execute_impl(
  const OutType&   out,
  const InType&    in,
  int              howmany,
  const ExecSpace& space) const
{
  using Kokkos::ALL;
  using policy_t     = Kokkos::TeamPolicy<ExecSpace, Kokkos::IndexType<int>>;
  using dev_buffer_t = memory_pool<typename ExecSpace::memory_space>;
  using deleter_t    = BufferDeleter<T, dev_buffer_t>;

  constexpr int tag = 64 << 8;

  auto scope = Kokkos::Profiling::ScopedRegion("TransposeMPI p2p");

  // reserve nproc send and recv requests respectively
  auto const           nproc       = base_t::np;
  auto const           rank        = impl->comm.rank();
  auto const           n_shm_procs = impl->shm_comm.size();
  mpipp::irequest_pool reqs;
  reqs.reserve(nproc * 2);

  size_t const blk_size = n0np * n1np * howmany;
  // Round up the block size, ensuring that each block is aligned to the cache
  // line size to avoid false sharing
  auto const blk_stride = Impl::aligned_bytes(blk_size * sizeof(T)) / sizeof(T);
  auto const recv_buf_size = blk_stride * nproc;
  auto const send_buf_size = blk_stride * (nproc - n_shm_procs);
  auto const ptr_send_buf  = std::unique_ptr<T, deleter_t>(
    (T*)dev_buffer_t().allocate(send_buf_size * sizeof(T)),
    deleter_t{send_buf_size * sizeof(T)});
  nonstd::span<T> send_buf(ptr_send_buf.get(), send_buf_size);
  nonstd::span<T> recv_buf = impl->get_shm_buffer(recv_buf_size);

  // Signal to peers that our shared buffer is ready for this generation's data
  // and our shared buffer is safe to overwrite
  auto const my_shm_rank = impl->shm_comm.rank();
  impl->generation.at(my_shm_rank)->fetch_add(1, std::memory_order_release);

  // Pack data one by one
  int              remaining_recvs{0};
  std::vector<int> shm_flags(n_shm_procs, 0);

  auto pack_policy = policy_t(space, n1np * howmany, Kokkos::AUTO());
  for (int i = 0, i_non_shm = 0; i < nproc; ++i) {
    auto [dst, src] = get_next_peers(i, rank, nproc);
    auto dst_in_shm = impl->peer_ranks_in_shm.at(dst);
    auto src_in_shm = impl->peer_ranks_in_shm.at(src);

    if (src_in_shm == MPI_UNDEFINED) {
      // Post receive if the source is not in shared memory
      reqs.push(mpipp::irecv(
        recv_buf.subspan(src * blk_stride, blk_size), src, tag, impl->comm));
    } else {
      reqs.push(mpipp::irecv(shm_flags.at(src_in_shm), src, tag, impl->comm));
    }
    ++remaining_recvs;

    nonstd::span<T> pack_blk{};
    if (dst_in_shm != MPI_UNDEFINED) {
      // Directly pack to the peer buffer if the destination is in shared memory
      // Wait until the peer has consumed the previous generation's data before
      // packing into its shared memory buffer
      using Clock                 = std::chrono::steady_clock;
      constexpr auto timeout      = std::chrono::seconds(30);
      auto const     t0           = Clock::now();
      uint64_t const expected_gen = impl->last_generation + 1;
      while (impl->generation.at(dst_in_shm)->load(std::memory_order_acquire)
             != expected_gen) {
        std::this_thread::yield();
        if (Clock::now() - t0 > timeout) {
          auto msg = fmt::format("Timeout waiting for peer shm_rank={} to "
                                 "advance generation (stuck at {})",
                                 dst_in_shm,
                                 impl->last_generation);
          throw std::runtime_error(msg);
        }
      }
      pack_blk =
        impl->peer_buffers.at(dst_in_shm).subspan(rank * blk_stride, blk_size);
    } else {
      pack_blk = send_buf.subspan(i_non_shm * blk_stride, blk_size);
      ++i_non_shm;
    }
    pack_functor(pack_blk.data(), in, dst, n0np, n1np, howmany)
      .execute(pack_policy);

    if (dst_in_shm == MPI_UNDEFINED) {
      // Post send for local data if the destination is not in shared memory
      reqs.push(mpipp::isend(pack_blk, dst, tag, impl->comm));
    } else {
      // Use send to notify the destination process
      reqs.push(mpipp::isend(1, dst, tag, impl->comm));
    }
  }

  // unpack data from other processes
  auto const unpack_policy = policy_t(space, howmany, Kokkos::AUTO());
  for (auto finished_req = reqs.waitsome(); finished_req.size() > 0;
       finished_req      = reqs.waitsome()) {
    for (auto i_finished : finished_req) {
      if (i_finished % 2 != 0) continue; // skip send completions

      auto m   = reqs.get_status(i_finished).source();
      auto blk = recv_buf.subspan(m * blk_stride, blk_size);
      Kokkos::parallel_for(
        "unpack",
        unpack_policy,
        unpack_functor(out, blk.data(), m, n0np, n1np, howmany, Op{}));
      --remaining_recvs;
    }
  }
  if (remaining_recvs != 0) {
    throw std::runtime_error("MPI receive incomplete");
  }

  // Advance last_generation to match the generation we just completed
  ++(impl->last_generation);
}

template<class T, class ExecSpace>
void TransposerMPIPoint2PointSHM<T, ExecSpace>::execute_impl(
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
