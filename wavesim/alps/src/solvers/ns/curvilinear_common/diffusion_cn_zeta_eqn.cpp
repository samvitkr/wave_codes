#include "diffusion_cn_zeta_eqn.h"

#include <mpipp/collectives.h>
#include <solvers/ns/curvilinear_common/diffusion_cn_zeta/jacobi.h>
#include <solvers/ns/curvilinear_common/diffusion_cn_zeta/pcr.h>
#include <solvers/ns/curvilinear_common/diffusion_cn_zeta/tdma.h>
#include <solvers/ns/curvilinear_common/diffusion_cn_zeta_coeff.h>
#include <spectral/spectral.h>

#include <stdexcept>
#include <type_traits>

namespace alps::solver {

template<class VarLoc>
void DiffusionCNZetaEqn<VarLoc>::solve(
  MDView<Real***> const&               dst,
  MDView<Real const***> const&         rhs,
  BottomWaveMesh const&                mesh,
  Kokkos::DefaultExecutionSpace const& stream) const
{
  auto const nx     = mesh.extents()[0];
  auto const ny     = mesh.extents()[1];
  auto const n_eqns = [&] {
    auto const is_top = mesh.grid.comm().is_last(2);
    if constexpr (std::is_same_v<VarLoc, NodePt>) {
      return is_top ? mesh.extent(2) - 1 : mesh.extent(2);
    }
    return mesh.extent(2);
  }();

  if (rhs.extent_int(0) < nx || rhs.extent_int(1) < ny
      || rhs.extent_int(2) < n_eqns || dst.extent_int(0) < nx
      || dst.extent_int(1) < ny || dst.extent_int(2) < n_eqns) {
    throw std::invalid_argument("CN rhs/solution size mismatch");
  }

  // distributed tridiagonal solvers require at least 2 points, otherwise fall
  // back to TDMA solver with point-to-point communication
  int reduced_ok_local = n_eqns >= 2 ? 1 : 0;
  int reduced_ok       = 0;
  mpipp::allreduce(reduced_ok_local,
                   reduced_ok,
                   mpipp::min<int>(),
                   mesh.grid.comm().axis_comm[2]);

  auto const algorithm =
    reduced_ok == 0 ? DiffusionCNZetaAlgorithm::TDMA : options.algorithm;
  DiffusionCNZetaCoeffProvider<VarLoc, BottomWaveMesh> coeff(
    mesh, bottom_bc, top_bc, (Real)alpha);
  auto const& comm = mesh.grid.comm().axis_comm[2];
  switch (algorithm) {
    case DiffusionCNZetaAlgorithm::TDMA:
      solve_diffusion_cn_zeta_tdma(
        dst, rhs, coeff, nx, ny, n_eqns, comm, stream);
      break;
    case DiffusionCNZetaAlgorithm::PCR:
      solve_diffusion_cn_zeta_pcr(
        dst, rhs, coeff, nx, ny, n_eqns, comm, stream);
      break;
    case DiffusionCNZetaAlgorithm::Jacobi:
      solve_diffusion_cn_zeta_jacobi(
        dst, rhs, coeff, nx, ny, n_eqns, options, comm, stream);
      break;
  }
}

template<>
DiffusionCNZetaEqn<CenterPt>::DiffusionCNZetaEqn(
  double                       alph,
  DiffusionCNBCType            top_bc_type,
  DiffusionCNBCType            bottom_bc_type,
  DiffusionCNZetaSolverOptions options_)
  : alpha{alph}
  , top_bc{top_bc_type}
  , bottom_bc{bottom_bc_type}
  , options{options_}
{}

template<>
DiffusionCNZetaEqn<NodePt>::DiffusionCNZetaEqn(
  double                       alph,
  DiffusionCNBCType            top_bc_type,
  DiffusionCNBCType            bottom_bc_type,
  DiffusionCNZetaSolverOptions options_)
  : alpha{alph}
  , top_bc{top_bc_type}
  , bottom_bc{bottom_bc_type}
  , options{options_}
{
  if (top_bc_type != DiffusionCNBCType::Dirichlet
      || bottom_bc_type != DiffusionCNBCType::Dirichlet) {
    throw std::invalid_argument("diffusion CN in node requires Dirichlet BCs");
  }
}

template class DiffusionCNZetaEqn<CenterPt>;
template class DiffusionCNZetaEqn<NodePt>;

} // namespace alps::solver
