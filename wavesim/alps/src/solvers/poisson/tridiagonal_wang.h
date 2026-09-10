#pragma once

#include "tridiagonal_solver.h"
#include <common/async/streams.h>

#include <Kokkos_Core.hpp>
#include <mpipp/comm.h>

#include <utility>
#include <vector>

namespace alps {
namespace solver {

template<typename ValueT>
class TridiagonalWang final : public TridiagonalSolver<ValueT>
{
 private:
  using base_t = TridiagonalSolver<ValueT>;
  using base_t::comm_;
  using base_t::nproc_;
  using base_t::rank_;
  using typename base_t::coeff_t;
  using typename base_t::solution_t;

  // Data type used for x in reduced system, a smaller data type can marginally
  // reduce solution time but may cause stability issues
  // using ReduceX_t = ValueT;
  using ReduceX_t = Kokkos::Experimental::bhalf_t;

 public:
  TridiagonalWang(int                        batch_count_1,
                  int                        batch_count_2,
                  int                        n_eqns,
                  mpipp::communicator const& comm);

  ~TridiagonalWang() override;

  void
  setup_impl(coeff_t const& d, coeff_t const& dl, coeff_t const& du) override;
  void solve_impl(solution_t const& x,
                  coeff_t const&    d,
                  coeff_t const&    dl,
                  coeff_t const&    du) override;

  void final_substitution(
    const solution_t&                                          x,
    const coeff_t&                                             d,
    std::pair<int, int>                                        y_range,
    const Kokkos::View<ReduceX_t const**, Kokkos::LayoutLeft>& rtemp_rank,
    const Kokkos::DefaultExecutionSpace&                       space) const;

  enum class ReductionMethod
  {
    MPI_Alltoall
  };

  ReductionMethod reduction_method{ReductionMethod::MPI_Alltoall};

  void reduce_and_backward_alltoall(
    solution_t                       x,
    coeff_t                          d,
    std::vector<std::pair<int, int>> ny_chunk_distribution,
    std::true_type /*mpi_can_access_default_space*/) const;
  void reduce_and_backward_alltoall(
    solution_t                       x,
    coeff_t                          d,
    std::vector<std::pair<int, int>> ny_chunk_distribution,
    std::false_type /*mpi_can_access_default_space*/) const;

  void tune_algorithms(const coeff_t& d, const coeff_t& dl, const coeff_t& du);

  using host_space = Kokkos::SharedHostPinnedSpace;
  Kokkos::View<ValueT***, Kokkos::LayoutLeft, Kokkos::SharedSpace> f_gather;
  Kokkos::View<ValueT***, Kokkos::LayoutLeft, Kokkos::SharedSpace> g_gather;
  Kokkos::View<ValueT***, Kokkos::LayoutLeft>                      f_, t_;
  int up_id, down_id;

  Kokkos::View<ReduceX_t**, Kokkos::LayoutLeft, host_space> xn_host;
  Kokkos::View<ReduceX_t**, Kokkos::LayoutLeft>             xn_dev;

  static constexpr int minimum_ny_chunk_size =
    16; // Adjustable parameter for achieving good latency hiding
  std::vector<std::pair<int, int>> ny_chunk_distribution_;

  StreamPool<Kokkos::DefaultExecutionSpace> stream_pool_;
};

} // namespace solver
} // namespace alps
