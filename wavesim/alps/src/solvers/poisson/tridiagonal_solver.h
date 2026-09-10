#pragma once

#include <common/base/logging_fwd.h>

#include <Kokkos_Core_fwd.hpp>
#include <mpipp/comm.h>

#include <memory>

// forward declaration
namespace Kokkos {
struct LayoutLeft;
}

namespace alps {
namespace solver {

template<typename ValueT>
class TridiagonalSolver
{
 public:
  using coeff_t    = Kokkos::View<ValueT***, Kokkos::LayoutLeft>;
  using solution_t = Kokkos::View<ValueT***, Kokkos::LayoutLeft>;

  TridiagonalSolver(int                        batch_count_1,
                    int                        batch_count_2,
                    int                        n_eqns,
                    mpipp::communicator const& comm);

  virtual ~TridiagonalSolver();

  /// One time setup of the solver
  void setup(coeff_t const& d, coeff_t const& dl, coeff_t const& du);

  /// Solve for the solution
  void solve(solution_t const& x,
             coeff_t const&    d,
             coeff_t const&    dl,
             coeff_t const&    du);

  /// Indicating whether the solver destroys the coefficient matrix during solve
  bool is_coeff_overwritten{false};

  int                 n1_;
  int                 n2_;
  int                 batch_;
  int                 nz_;
  mpipp::communicator comm_;
  int                 rank_;
  int                 nproc_;

  virtual void
  setup_impl(coeff_t const& d, coeff_t const& dl, coeff_t const& du) = 0;

  virtual void solve_impl(solution_t const& x,
                          coeff_t const&    d,
                          coeff_t const&    dl,
                          coeff_t const&    du) = 0;

 protected:
  Logger logger_;
};

enum class TridiagonalAlgorithm
{
  MpiLU    = 0,
  SeqLU    = 1,
  CuSparse = 2,
  Wang     = 3,
  Default  = -1
};

template<typename ValueT>
std::unique_ptr<TridiagonalSolver<ValueT>> create_tridiagonal_solver(
  int                        batch_count_1,
  int                        batch_count_2,
  int                        n_eqns,
  mpipp::communicator const& comm,
  TridiagonalAlgorithm       method = TridiagonalAlgorithm::Default);

} // namespace solver
} // namespace alps
