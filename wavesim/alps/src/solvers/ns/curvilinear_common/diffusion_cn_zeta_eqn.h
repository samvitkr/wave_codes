#pragma once

#include <common/container/view_types.h>
#include <common/real_type.h>
#include <solvers/mesh/mesh_fwd.h>
#include <solvers/ns/channel/diffusion_cn.h>

#include <Kokkos_Core_fwd.hpp>

namespace alps::solver {

enum class DiffusionCNZetaAlgorithm
{
  TDMA,
  PCR,
  Jacobi,
};

struct DiffusionCNZetaSolverOptions
{
  DiffusionCNZetaAlgorithm algorithm{DiffusionCNZetaAlgorithm::TDMA};
  int                      jacobi_max_iters{100};
  double                   jacobi_abs_tol{1e-8};
  double                   jacobi_rel_tol{1e-12};
};

template<class VarLoc>
class DiffusionCNZetaEqn;

template<class VarLoc>
class DiffusionCNZetaEqn
{
 public:
  // Solves the 1D implicit CN line system along zeta for each (i,j):
  //
  //   [I - alpha * D_zeta( g33 * D_zeta )] u^{n+1} = rhs,
  //
  // where g33 = g^{33}(zeta, exr, eyr, invJ) is the contravariant metric term.
  // The u field is actually J^{-1}u, when used in the flow solver.
  //
  // Note: VarLoc = NodePt requires Dirichlet BCs at both boundaries; the
  // constructor throws for any other BC type.
  DiffusionCNZetaEqn(double                       alph,
                     DiffusionCNBCType            top_bc_type,
                     DiffusionCNBCType            bottom_bc_type,
                     DiffusionCNZetaSolverOptions options = {});

  void initialize() const {}

  void solve(MDView<Real***> const&               dst,
             MDView<Real const***> const&         rhs,
             BottomWaveMesh const&                mesh,
             Kokkos::DefaultExecutionSpace const& stream) const;

  double                       alpha;
  DiffusionCNBCType            top_bc;
  DiffusionCNBCType            bottom_bc;
  DiffusionCNZetaSolverOptions options;
};
} // namespace alps::solver
