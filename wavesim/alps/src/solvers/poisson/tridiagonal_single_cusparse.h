#pragma once

#include "tridiagonal_solver.h"

#include <mpipp/comm.h>

#include <memory>

namespace alps {
namespace solver {

namespace detail {
template<typename T>
class GTSVImpl;
} // namespace detail

template<typename ValueT>
class TridiagonalCuSparse final : public TridiagonalSolver<ValueT>
{
 private:
  using base_t = TridiagonalSolver<ValueT>;
  using typename base_t::coeff_t;
  using typename base_t::solution_t;

 public:
  TridiagonalCuSparse(int                        batch_count_1,
                      int                        batch_count_2,
                      int                        n_eqns,
                      mpipp::communicator const& comm);

  ~TridiagonalCuSparse() override;

  void
  setup_impl(coeff_t const& d, coeff_t const& dl, coeff_t const& du) override;
  void solve_impl(solution_t const& x,
                  coeff_t const&    d,
                  coeff_t const&    dl,
                  coeff_t const&    du) override;

 private:
  using impl_t = detail::GTSVImpl<ValueT>;
  std::unique_ptr<impl_t> impl_;
};

} // namespace solver
} // namespace alps
