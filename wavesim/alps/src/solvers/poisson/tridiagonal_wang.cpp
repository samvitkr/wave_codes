#include "tridiagonal_wang.h"

#include <common/async/streams.h>
#include <common/base/logging.h>
#include <common/base/macros.h>
#include <common/container/view_utils.h>
#include <common/device/device_traits.h>
#include <common/kokkos_abstraction/exec_policy.h>
#include <common/kokkos_abstraction/pool_space.h>
#include <common/utils/timeout_mpi_barrier.h>
#include <decomp/ghost_cell_exchange.h>

#include <Kokkos_Core.hpp>
#include <mpipp/collectives.h>
#include <mpipp/comm.h>
#include <mpipp/point2point.h>

#include <stdexcept>
#include <vector>

namespace alps {
namespace solver {

template<class DataType, class... Properties>
using MDView   = Kokkos::View<DataType, Kokkos::LayoutLeft, Properties...>;
using member_t = GridPolicy<>::member_type;
using Kokkos::ALL;
using Kokkos::parallel_for;
using Kokkos::subview;
using Kokkos::TeamThreadRange;

namespace {
auto constexpr default_tile_size = []() -> Kokkos::Array<range_index_t, 2> {
  using Device = Kokkos::DefaultExecutionSpace;
  if constexpr (is_cuda_execution_space_v<Device>) return {32, 4};
  if constexpr (is_hip_execution_space_v<Device>) return {64, 8};
  return {0, 0};
}();
} // namespace

template<typename ValueT>
TridiagonalWang<ValueT>::TridiagonalWang(int batch_count_1,
                                         int batch_count_2,
                                         int n_eqns,
                                         const mpipp::communicator& comm)
  : base_t(batch_count_1, batch_count_2, n_eqns, comm)
  , up_id{(rank_ + 1) < nproc_ ? rank_ + 1 : MPI_PROC_NULL}
  , down_id{rank_ > 0 ? rank_ - 1 : MPI_PROC_NULL}
  , stream_pool_{6}
{
  if (n_eqns <= 2) {
    throw std::invalid_argument("The Wang algorithm must have at least 2 nodes "
                                "per processor.");
  }
}

template<typename ValueT>
void TridiagonalWang<ValueT>::setup_impl(coeff_t const& d,
                                         coeff_t const& dl,
                                         coeff_t const& du)
{
  const auto n1 = this->n1_;
  const auto n2 = this->n2_;
  const auto n3 = this->nz_;

  if (d.extent_int(0) < n1 || d.extent_int(1) < n2 || d.extent_int(2) != n3) {
    throw std::invalid_argument("Poisson solver coefficient d size mismatch.");
  }
  if (dl.extent_int(0) < n1 || dl.extent_int(1) < n2
      || dl.extent_int(2) != n3) {
    throw std::invalid_argument("Poisson solver coefficient dl size mismatch.");
  }
  if (du.extent_int(0) < n1 || du.extent_int(1) < n2
      || du.extent_int(2) != n3) {
    throw std::invalid_argument("Poisson solver coefficient du size mismatch.");
  }

  // A few policies that are used frequently below
  const auto policy_2d   = LoopPolicy<2>({0, 0}, {n1, n2});
  const auto policy_n3_1 = LoopPolicy<3>({0, 0, 0}, {n1, n2, n3 - 1});
  const auto policy_n2   = GridPolicy<>(n2, Kokkos::AUTO);

  // eliminate c, and generate f
  // 1. pass top layer of du (upper diagonal) to the upper process (stored in
  // du_bottom)
  const auto                            du_top = subview(du, ALL, ALL, n3 - 1);
  MDView<ValueT**, default_memory_pool> du_bottom("du_bottom", du_top.layout());
  exchange(du_top, up_id, du_bottom, down_id, comm_, HaloMsgTagPrefix::Z);

  // 2. elimination
  f_ = MDView<ValueT***, default_memory_pool>(
    Kokkos::view_alloc("matF", Kokkos::WithoutInitializing),
    d.extent(0),
    d.extent(1),
    n3);
  const auto& f = f_;
  if (down_id != MPI_PROC_NULL) {
    parallel_for(
      policy_2d, KOKKOS_LAMBDA(int i, int j) { f(i, j, 0) = dl(i, j, 0); });
  }
  parallel_for(
    policy_n2, KOKKOS_LAMBDA(member_t team) {
      int j = team.league_rank();
      for (int k = 1; k < n3; ++k) {
        parallel_for(
          TeamThreadRange(team, n1), KOKKOS_TR_LAMBDA(int& i) {
            dl(i, j, k) /= d(i, j, k - 1);
            d(i, j, k) -= dl(i, j, k) * du(i, j, k - 1);
            f(i, j, k) = -f(i, j, k - 1) * dl(i, j, k);
          });
      }
    });
  policy_n2.space().fence();

  // all gather f
  f_gather = MDView<ValueT***, Kokkos::SharedSpace>(
    Kokkos::view_alloc("matF_all", Kokkos::WithoutInitializing),
    f.extent(0),
    f.extent(1),
    nproc_);
  {
    const auto f_top = subview(f, ALL, ALL, n3 - 1);
    const auto f_top_host =
      Kokkos::create_mirror_view(Kokkos::WithoutInitializing, f_top);
    deep_copy(f_top_host, f_top);
    mpipp::allgather(nonstd::span(f_top_host.data(), f_top_host.span()),
                     f_gather.data(),
                     comm_);
  }

  // eliminate b, and generate g
  // 1. pass top layer of d to the upper process (stored in d_bottom)
  const auto                            d_top = subview(d, ALL, ALL, n3 - 1);
  MDView<ValueT**, default_memory_pool> d_bottom("d_bottom", d_top.layout());
  exchange(d_top, up_id, d_bottom, down_id, comm_, HaloMsgTagPrefix::Z);

  // 2. eliminate b
  MDView<ValueT***, default_memory_pool> g(
    Kokkos::view_alloc("", Kokkos::WithoutInitializing), f_.layout());
  parallel_for(
    policy_2d,
    KOKKOS_LAMBDA(int i, int j) { g(i, j, n3 - 2) = du(i, j, n3 - 2); });
  parallel_for(
    policy_n2, KOKKOS_LAMBDA(member_t team) {
      int j = team.league_rank();
      for (int k = n3 - 2; k >= 1; --k) {
        parallel_for(
          TeamThreadRange(team, n1), KOKKOS_TR_LAMBDA(int& i) {
            du(i, j, k - 1) /= d(i, j, k);
            g(i, j, k - 1) = -du(i, j, k - 1) * g(i, j, k);
            f(i, j, k - 1) -= du(i, j, k - 1) * f(i, j, k);
          });
      }
    });
  MDView<ValueT**, default_memory_pool> g_bottom("g_bottom", d_top.layout());
  parallel_for(
    policy_2d, KOKKOS_LAMBDA(int i, int j) {
      du_bottom(i, j) /= d(i, j, 0);
      g_bottom(i, j) = -du_bottom(i, j) * g(i, j, 0);
      d_bottom(i, j) -= du_bottom(i, j) * f(i, j, 0);
    });
  policy_2d.space().fence();

  // pass d_bottom down to the upper layer of d on the lower process
  exchange(d_bottom, down_id, d_top, up_id, comm_, HaloMsgTagPrefix::Z);

  MDView<ValueT***, memory_pool<host_space>> d_gather(
    Kokkos::view_alloc("d_gather", Kokkos::WithoutInitializing),
    f_gather.layout());
  {
    auto const d_top_host =
      Kokkos::create_mirror_view(Kokkos::WithoutInitializing, d_top);
    Kokkos::deep_copy(d_top_host, d_top);
    mpipp::allgather(
      nonstd::span(d_top_host.data(), d.stride(2)), d_gather.data(), comm_);
  }

  exchange(du_bottom, down_id, du_top, up_id, comm_, HaloMsgTagPrefix::Z);

  g_gather = MDView<ValueT***, Kokkos::SharedSpace>(
    Kokkos::view_alloc("matG_all", Kokkos::WithoutInitializing),
    f_gather.layout());
  {
    const auto g_bottom_host =
      Kokkos::create_mirror_view(Kokkos::WithoutInitializing, g_bottom);
    deep_copy(g_bottom_host, g_bottom);
    mpipp::allgather(
      nonstd::span(g_bottom_host.data(), d.stride(2)), g_gather.data(), comm_);
  }

  using HostES = Kokkos::HostSpace::execution_space;
  for (int k = 0; k < nproc_ - 1; ++k) {
    parallel_for(LoopPolicy<2, HostES>({0, 0}, {n1, n2}),
                 [this, k, d_gather](int i, int j) {
                   f_gather(i, j, k + 1) /= d_gather(i, j, k);
                   d_gather(i, j, k + 1) -=
                     f_gather(i, j, k + 1) * g_gather(i, j, k + 1);
                 });
  }

  parallel_for(LoopPolicy<3, HostES>({0, 0, 1}, {n1, n2, nproc_}),
               [this, d_gather](int i, int j, int k) {
                 g_gather(i, j, k) /= d_gather(i, j, k);
               });
  HostES().fence();

  if (rank_ != 0) {
    const auto d_gather_rank_host = subview(d_gather, ALL, ALL, rank_ - 1);
    const auto d_gather_rank =
      Kokkos::create_mirror_view(default_memory_pool(), d_gather_rank_host);
    deep_copy(d_gather_rank, d_gather_rank_host);
    parallel_for(
      policy_n3_1, KOKKOS_LAMBDA(int i, int j, int k) {
        f(i, j, k) /= d_gather_rank(i, j);
      });
    policy_n3_1.space().fence();
  }

  auto d_gather_rank_host = subview(d_gather, ALL, ALL, rank_);
  auto d_gather_rank =
    Kokkos::create_mirror_view(default_memory_pool(), d_gather_rank_host);
  deep_copy(d_gather_rank, d_gather_rank_host);
  parallel_for(
    policy_n3_1,
    KOKKOS_LAMBDA(int i, int j, int k) { g(i, j, k) /= d_gather_rank(i, j); });
  parallel_for(
    policy_2d,
    KOKKOS_LAMBDA(int i, int j) { g_bottom(i, j) /= d_gather_rank(i, j); });
  policy_2d.space().fence();

  t_ = MDView<ValueT***, default_memory_pool>(
    Kokkos::view_alloc("mat_t", Kokkos::WithoutInitializing),
    d.extent(0),
    d.extent(1),
    n3 - 1);
  auto& t = t_;
  if (rank_ != 0) {
    parallel_for(
      policy_n3_1, KOKKOS_LAMBDA(int i, int j, int k) {
        t(i, j, k) = g(i, j, k) - f(i, j, k) * g_bottom(i, j);
      });
  } else {
    parallel_for(
      policy_n3_1,
      KOKKOS_LAMBDA(int i, int j, int k) { t(i, j, k) = g(i, j, k); });
  }
  policy_n3_1.space().fence();

  Kokkos::deep_copy(d_top, d_gather_rank);

  xn_host = MDView<ReduceX_t**, PoolSpace<host_space>>(
    Kokkos::view_alloc("xn_host", Kokkos::WithoutInitializing),
    n1,
    d_top.extent(1));

  tune_algorithms(d, dl, du);

  Kokkos::fence();
}

template<typename ValueT>
TridiagonalWang<ValueT>::~TridiagonalWang() = default;

template<typename ValueT>
void TridiagonalWang<ValueT>::solve_impl(solution_t const& x,
                                         coeff_t const&    d,
                                         coeff_t const&    dl,
                                         coeff_t const&    du)
{
  const auto n1 = this->n1_;
  const auto n2 = this->n2_;
  const auto n3 = this->nz_;
  const auto x0 = subview(x, ALL, ALL, 0);
  const auto xn = subview(x, ALL, ALL, n3 - 1);

  MDView<ValueT**, default_memory_pool> r(
    Kokkos::view_alloc("", Kokkos::WithoutInitializing), x0.layout());

  MDView<ValueT**, host_space> r_host =
    (!alps::mpi_can_access_v<typename decltype(r)::memory_space>)
      ? MDView<ValueT**, memory_pool<host_space>>("r_host_mirror", r.layout())
      : MDView<ValueT**, memory_pool<host_space>>();
  auto req =
    alps::mpi_can_access_v<typename decltype(r)::memory_space>
      ? mpipp::irecv(nonstd::span(r.data(), x.stride(2)), up_id, 3, comm_)
      : mpipp::irecv(nonstd::span(r_host.data(), x.stride(2)), up_id, 3, comm_);

  if (!check_same_layout_and_offset(x, d, dl, du)) {
    throw std::runtime_error("Mismatched extents in x, d, dl, and du");
  }

  // Perform the local solve
  auto policy2d = LoopPolicy<2>(
    stream_pool_.get_stream(1), {0, 0}, {n1, n2}, default_tile_size);
  parallel_for(
    "Wang local1", policy2d, KOKKOS_LAMBDA(int i, int j) {
      auto offset = i + j * (int)x.stride(1);
      auto stride = x.stride(2);
      for (int k = 1; k < n3; ++k) {
        auto k_stride = k * stride;
        x.data()[offset + k_stride] -=
          dl.data()[offset + k_stride] * x.data()[offset + k_stride - stride];
      }
      for (int k = n3 - 3; k >= 0; --k) {
        auto k_stride = k * stride;
        x.data()[offset + k_stride] -=
          du.data()[offset + k_stride] * x.data()[offset + k_stride + stride];
      }
    });

  if (alps::mpi_can_access_v<typename decltype(x0)::memory_space>) {
    policy2d.space().fence();
    mpipp::send(nonstd::span(x0.data(), x.stride(2)), down_id, 3, comm_);
  } else {
    const auto x0_host = Kokkos::create_mirror_view(
      Kokkos::WithoutInitializing, memory_pool<host_space>(), x0);
    policy2d.space().fence();
    if (down_id != MPI_PROC_NULL) {
      deep_copy(x0_host, x0);
    }
    mpipp::send(nonstd::span(x0_host.data(), x.stride(2)), down_id, 3, comm_);
  }

  // Prepare the reduced system
  if constexpr (std::is_same_v<
                  pool_base_space_t<default_memory_pool>,
                  typename std::decay_t<decltype(xn_host)>::memory_space>) {
    // when the underlying type of the default_memory_pool is the same as the
    // memory space of xn_host, we can directly use xn_host. This is the case
    // when the default_memory_pool is pool_space<Kokkos::HostSpace>. This
    // branch is needed because Kokkos::create_mirror_view does not recognize
    // default_memory_pool and Kokkos::HostSpace as the same memory space and
    // allocate a new view.
    xn_dev = decltype(xn_dev)(xn_host.data(), xn_host.layout());
  } else {
    xn_dev = Kokkos::create_mirror_view(
      Kokkos::WithoutInitializing, default_memory_pool(), xn_host);
  }
  const auto& xn_d = xn_dev;
  req.wait();
  if (rank_ != nproc_ - 1) {
    if constexpr (!alps::mpi_can_access_v<typename decltype(r)::memory_space>) {
      deep_copy(policy2d.space(), r, r_host);
    }
    const auto du_top = subview(du, ALL, ALL, n3 - 1);
    parallel_for(
      "Wang prepare reduced", policy2d, KOKKOS_LAMBDA(int i, int j) {
        xn_d(i, j) = ReduceX_t(xn(i, j) - du_top(i, j) * r(i, j));
      });
  } else {
    parallel_for(
      "Wang prepare reduced", policy2d, KOKKOS_LAMBDA(int i, int j) {
        xn_d(i, j) = (ReduceX_t)xn(i, j);
      });
  }

  policy2d.space().fence();
  if (reduction_method == ReductionMethod::MPI_Alltoall) {
    reduce_and_backward_alltoall(x,
                                 d,
                                 ny_chunk_distribution_,
                                 alps::mpi_can_access<default_memory_pool>());
  }

  xn_dev = {}; // deallocate xn_dev
}

template<typename ValueT>
void TridiagonalWang<ValueT>::final_substitution(
  const solution_t&                    x,
  const coeff_t&                       d,
  std::pair<int, int>                  y_range,
  const MDView<ReduceX_t const**>&     rtemp_rank,
  const Kokkos::DefaultExecutionSpace& space) const
{
  const auto  xn   = subview(x, ALL, ALL, this->nz_ - 1);
  const auto& xn_d = xn_dev;
  const auto  policy_n3_1 =
    LoopPolicy<3>(space,
                  {0, y_range.first, 0},
                  {this->n1_, y_range.second, this->nz_ - 1},
                  {default_tile_size[0],
                   default_tile_size[0] == 0 ? 0 : 1,
                   default_tile_size[1]});
  const auto policy2d = LoopPolicy<2>(
    space, {0, y_range.first}, {this->n1_, y_range.second}, default_tile_size);

  if (rank_ != 0) {
    const auto& f = this->f_;
    const auto& t = this->t_;
    parallel_for(
      "Wang backward", policy_n3_1, KOKKOS_LAMBDA(int i, int j, int k) {
        auto v = x(i, j, k) - f(i, j, k) * ValueT(rtemp_rank(i, j))
               - ValueT(xn_d(i, j)) * t(i, j, k);
        x(i, j, k) = v / d(i, j, k);
      });
  } else {
    const auto& t = this->t_;
    parallel_for(
      "Wang backward bottom", policy_n3_1, KOKKOS_LAMBDA(int i, int j, int k) {
        auto v     = x(i, j, k) - ValueT(xn_d(i, j)) * t(i, j, k);
        x(i, j, k) = v / d(i, j, k);
      });
  }

  const auto d_top = subview(d, ALL, ALL, this->nz_ - 1);
  parallel_for(
    "Wang top", policy2d, KOKKOS_LAMBDA(int i, int j) {
      xn(i, j) = ValueT(xn_d(i, j)) / d_top(i, j);
    });
}

template class TridiagonalWang<double>;
template class TridiagonalWang<float>;

} // namespace solver
} // namespace alps
