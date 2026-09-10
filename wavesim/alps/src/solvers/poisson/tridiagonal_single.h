#pragma once

#include "tridiagonal_solver.h"

#include <mpipp/comm.h>

namespace alps {
namespace solver {
template<typename ValueT>
class TridiagonalSeqLU final : public TridiagonalSolver<ValueT>
{
 private:
  using base_t = TridiagonalSolver<ValueT>;
  using typename base_t::coeff_t;
  using typename base_t::solution_t;

 public:
  TridiagonalSeqLU(int                        batch_count_1,
                   int                        batch_count_2,
                   int                        n_eqns,
                   mpipp::communicator const& comm);

  void
  setup_impl(coeff_t const& d, coeff_t const& dl, coeff_t const& du) override;
  void solve_impl(solution_t const& x,
                  coeff_t const&    d,
                  coeff_t const&    dl,
                  coeff_t const&    du) override;
};

} // namespace solver
} // namespace alps
