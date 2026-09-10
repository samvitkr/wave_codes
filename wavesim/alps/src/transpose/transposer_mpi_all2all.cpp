#include "transposer_mpi_all2all.h"

#include <common/base/logging.h>
#include <common/base/macros.h>
#include <common/device/device_traits.h>
#include <common/kokkos_abstraction/pool_space.h>

#include <Kokkos_Core.hpp>
#include <Kokkos_Profiling_ScopedRegion.hpp>
#include <mpipp/collectives.h>

#include <memory>
#include <stdexcept>

namespace alps::transpose {

template<class DataType, class... Properties>
using MDView = Kokkos::View<DataType, Kokkos::LayoutLeft, Properties...>;

template<class T, typename ExecSpace>
TransposerMPIAll2All<T, ExecSpace>::TransposerMPIAll2All(
  mpipp::communicator comm,
  int                 n0,
  int                 n1,
  int /*max_nz_hint*/,
  TransposerOptions /*options*/)
  : base_t(n0, n1, comm.size())
  , comm_(comm)
{}

template<class T, typename ExecSpace>
TransposerMPIAll2All<T, ExecSpace>::~TransposerMPIAll2All() = default;

template class TransposerMPIAll2All<float, Kokkos::OpenMP>;
template class TransposerMPIAll2All<double, Kokkos::OpenMP>;

#if defined(KOKKOS_ENABLE_CUDA)
template class TransposerMPIAll2All<float, Kokkos::Cuda>;
template class TransposerMPIAll2All<double, Kokkos::Cuda>;
#elif defined(KOKKOS_ENABLE_HIP)
template class TransposerMPIAll2All<float, Kokkos::HIP>;
template class TransposerMPIAll2All<double, Kokkos::HIP>;
#endif

namespace {

template<typename InType, typename OutType>
struct pack_functor
{
  using T         = typename OutType::non_const_value_type;
  using ExecSpace = typename OutType::execution_space;
  using policy_t  = Kokkos::TeamPolicy<ExecSpace, Kokkos::IndexType<int>>;
  using member_t  = typename policy_t::member_type;

  typename InType::const_type                    in;
  MDView<T***, typename ExecSpace::memory_space> out;
  int                                            n0pp, n1pp, np;

  pack_functor(OutType const& output,
               InType const&  input,
               int            n_proc,
               int            n0_per_proc,
               int            n1_per_proc,
               int            batch)
    : in(input)
    , out(output.data(), n0_per_proc, n1_per_proc * batch, n_proc)
    , n0pp{n0_per_proc}
    , n1pp{n1_per_proc}
    , np{n_proc}
  {}

  KOKKOS_FUNCTION void operator()(const member_t& team) const
  {
    using Kokkos::parallel_for;

    const auto m = static_cast<int>(team.league_rank());
    const auto k = m / n1pp;
    const auto j = m % n1pp;
    parallel_for(
      Kokkos::TeamThreadRange(team, np), KOKKOS_TR_LAMBDA(int& p) {
        const auto offset = p * n0pp;
        parallel_for(
          Kokkos::ThreadVectorRange(team, n0pp),
          KOKKOS_TR_LAMBDA(int& i) { out(i, m, p) = in(i + offset, j, k); });
      });
  }

  void execute(ExecSpace const& stream) const
  {
    constexpr auto vector_len = []() constexpr -> int {
      if constexpr (is_cuda_execution_space_v<ExecSpace>) return 16;
      if constexpr (is_hip_execution_space_v<ExecSpace>) return 32;
      return 1;
    }();
    auto policy = policy_t(stream, out.extent(1), Kokkos::AUTO(), vector_len);
    Kokkos::parallel_for("pack", policy, *this);
  }
};

template<typename InType, typename OutType, typename Op>
struct unpack_functor
{
  using T         = typename InType::non_const_value_type;
  using ExecSpace = typename InType::execution_space;
  using policy_t  = Kokkos::TeamPolicy<ExecSpace, Kokkos::IndexType<int>>;
  using member_t  = typename policy_t::member_type;

  MDView<T const***, typename ExecSpace::memory_space> buffer;
  MDView<T***, typename ExecSpace::memory_space>       out;
  int                                                  n0pp, n1pp, howmany, np;

  unpack_functor(OutType const& output,
                 InType const&  input,
                 int            n_proc,
                 int            n0_per_proc,
                 int            n1_per_proc,
                 int            batch,
                 Op /*op*/)
    : buffer(input.data(), n0_per_proc, n1_per_proc, batch * n_proc)
    , out(output)
    , n0pp{n0_per_proc}
    , n1pp{n1_per_proc}
    , howmany{batch}
    , np{n_proc}
  {}

  KOKKOS_INLINE_FUNCTION void operator()(const member_t& team) const
  {
    using Kokkos::parallel_for;

    const auto  m            = static_cast<int>(team.league_rank());
    const auto  k            = m % howmany;
    const auto  p            = m / howmany;
    const auto* ptr_buffer_m = buffer.data() + m * buffer.stride(2);
    const auto  offset       = p * n1pp;
    parallel_for(
      Kokkos::TeamThreadRange(team, n0pp), KOKKOS_TR_LAMBDA(int& j) {
        parallel_for(
          Kokkos::ThreadVectorRange(team, n1pp), KOKKOS_TR_LAMBDA(int& i) {
            Op::apply(out(i + offset, j, k), ptr_buffer_m[j + i * n0pp]);
          });
      });
  }

  void execute(ExecSpace const& stream) const
  {
    constexpr auto vector_len = []() constexpr -> int {
      if constexpr (is_cuda_execution_space_v<ExecSpace>) return 8;
      if constexpr (is_hip_execution_space_v<ExecSpace>) return 32;
      return 1;
    }();
    auto policy = policy_t(stream, howmany * np, Kokkos::AUTO(), vector_len);
    Kokkos::parallel_for("unpack", policy, *this);
  }
};

template<typename T, typename MemSpace>
struct BufferDeleter
{
  std::size_t bytes;

  void operator()(T* ptr) const { MemSpace().deallocate(ptr, bytes); }
};
} // anonymous namespace

template<class T, class ExecSpace>
void TransposerMPIAll2All<T, ExecSpace>::execute_impl(
  const OutType&   out,
  const InType&    in,
  int              howmany,
  TransposeOps     op,
  const ExecSpace& space) const
{
  auto scope = Kokkos::Profiling::ScopedRegion("TransposeMPI all2all");

  using dev_buffer_t = memory_pool<typename ExecSpace::memory_space>;
  using deleter_t    = BufferDeleter<T, dev_buffer_t>;

  const auto nproc       = comm_.size();
  const auto blk_size    = n0np * n1np * howmany;
  const auto buffer_size = blk_size * nproc;
  const auto alloc_bytes = buffer_size * sizeof(T);

  // allocate send and recv buffer
  std::unique_ptr<T, deleter_t> send_buf_ptr(
    (T*)dev_buffer_t().allocate(alloc_bytes), deleter_t{alloc_bytes});
  std::unique_ptr<T, deleter_t> recv_buf_ptr(
    (T*)dev_buffer_t().allocate(alloc_bytes), deleter_t{alloc_bytes});
  MDView<T*, dev_buffer_t> send_buf(send_buf_ptr.get(), buffer_size);
  MDView<T*, dev_buffer_t> recv_buf(recv_buf_ptr.get(), buffer_size);

  // pack data
  pack_functor(send_buf, in, nproc, n0np, n1np, howmany).execute(space);
  space.fence();

  if (::alps::mpi_can_access_v<dev_buffer_t>) {
    mpipp::alltoall(nonstd::span(send_buf.data(), send_buf.span()),
                    blk_size,
                    recv_buf.data(),
                    comm_);
  } else {
    using host_buf_space = memory_pool<Kokkos::SharedHostPinnedSpace>;
    auto send_buf_host   = Kokkos::create_mirror(
      Kokkos::WithoutInitializing, host_buf_space(), send_buf);
    auto recv_buf_host = Kokkos::create_mirror(
      Kokkos::WithoutInitializing, host_buf_space(), recv_buf);
    Kokkos::deep_copy(space, send_buf_host, send_buf);
    space.fence();
    mpipp::alltoall(nonstd::span(send_buf_host.data(), send_buf_host.span()),
                    blk_size,
                    recv_buf_host.data(),
                    comm_);
    Kokkos::deep_copy(space, recv_buf, recv_buf_host);
    space.fence(); // ensure recv_buf_host is not released too early
  }

  // unpack data
  std::visit(
    [=, n0pp = this->n0np, n1pp = this->n1np, &space](auto&& op_) {
      using Op = std::decay_t<decltype(op_)>;
      unpack_functor(out, recv_buf, nproc, n0pp, n1pp, howmany, Op{})
        .execute(space);
    },
    op);
  space.fence();
}
} // namespace alps::transpose
