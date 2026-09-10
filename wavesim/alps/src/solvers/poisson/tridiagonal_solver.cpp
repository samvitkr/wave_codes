#include "tridiagonal_solver.h"

#include "tridiagonal_mpi_lu.h"
#include "tridiagonal_single.h"
#include "tridiagonal_single_cusparse.h"
#include "tridiagonal_wang.h"

#include <common/base/logging.h>
#include <mpipp/comm.h>

#include <memory>

namespace alps {
namespace solver {

template<typename ValueT>
TridiagonalSolver<ValueT>::TridiagonalSolver(int batch_count_1,
                                             int batch_count_2,
                                             int n_eqns,
                                             const mpipp::communicator& comm)
  : n1_{batch_count_1}
  , n2_{batch_count_2}
  , batch_{batch_count_1 * batch_count_2}
  , nz_{n_eqns}
  , comm_{comm}
  , rank_{comm.rank()}
  , nproc_{comm.size()}
  , logger_{alps::get_logger("tdsolver")}
{}

template<typename ValueT>
TridiagonalSolver<ValueT>::~TridiagonalSolver() = default;

template<typename ValueT>
void TridiagonalSolver<ValueT>::setup(coeff_t const& d,
                                      coeff_t const& dl,
                                      coeff_t const& du)
{
  logger_->trace("setup tridiagonal system");

  Kokkos::Profiling::pushRegion("TridiagonalSolver::setup");
  this->setup_impl(d, dl, du);
  Kokkos::Profiling::popRegion();
}

template<typename ValueT>
void TridiagonalSolver<ValueT>::solve(solution_t const& x,
                                      coeff_t const&    d,
                                      coeff_t const&    dl,
                                      coeff_t const&    du)
{
  logger_->trace("solve tridiagonal system");

  Kokkos::Profiling::pushRegion("TridiagonalSolver::solve");
  this->solve_impl(x, d, dl, du);
  Kokkos::Profiling::popRegion();
}

template<typename ValueT>
std::unique_ptr<TridiagonalSolver<ValueT>>
create_tridiagonal_solver(int                        batch_count_1,
                          int                        batch_count_2,
                          int                        n_eqns,
                          const mpipp::communicator& comm,
                          TridiagonalAlgorithm       method)
{
  switch (method) {
    case TridiagonalAlgorithm::MpiLU:
      return std::make_unique<TridiagonalMpiLU<ValueT>>(
        batch_count_1, batch_count_2, n_eqns, comm);
    case TridiagonalAlgorithm::SeqLU:
      return std::make_unique<TridiagonalSeqLU<ValueT>>(
        batch_count_1, batch_count_2, n_eqns, comm);
    case TridiagonalAlgorithm::CuSparse:
      return std::make_unique<TridiagonalCuSparse<ValueT>>(
        batch_count_1, batch_count_2, n_eqns, comm);
    case TridiagonalAlgorithm::Wang:
      return std::make_unique<TridiagonalWang<ValueT>>(
        batch_count_1, batch_count_2, n_eqns, comm);
    case TridiagonalAlgorithm::Default:
      if (comm.size() == 1) {
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
        return std::make_unique<TridiagonalCuSparse<ValueT>>(
          batch_count_1, batch_count_2, n_eqns, comm);
#else
        return std::make_unique<TridiagonalSeqLU<ValueT>>(
          batch_count_1, batch_count_2, n_eqns, comm);
#endif
      }
  }
  return std::make_unique<TridiagonalWang<ValueT>>(
    batch_count_1, batch_count_2, n_eqns, comm);
}

template class TridiagonalSolver<double>;
template class TridiagonalSolver<float>;

template std::unique_ptr<TridiagonalSolver<double>>
create_tridiagonal_solver(int,
                          int,
                          int,
                          const mpipp::communicator& comm,
                          TridiagonalAlgorithm);
template std::unique_ptr<TridiagonalSolver<float>>
create_tridiagonal_solver(int,
                          int,
                          int,
                          const mpipp::communicator& comm,
                          TridiagonalAlgorithm);

} // namespace solver
} // namespace alps
