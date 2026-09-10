#pragma once

#include <common/container/view_types.h>
#include <common/real_type.h>
#include <mpipp/comm.h>
#include <solvers/ns/curvilinear_common/diffusion_cn_zeta_eqn.h>

#include <Kokkos_Core_fwd.hpp>

namespace alps::solver {

template<typename CoeffFunctor>
void solve_diffusion_cn_zeta_jacobi(
  MDView<Real***> const&               dst,
  MDView<Real const***> const&         rhs,
  CoeffFunctor const&                  coeff,
  int                                  nx,
  int                                  ny,
  int                                  n_eqns,
  DiffusionCNZetaSolverOptions const&  options,
  mpipp::communicator const&           comm,
  Kokkos::DefaultExecutionSpace const& stream);

} // namespace alps::solver
