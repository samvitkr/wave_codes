#pragma once

#include <common/container/view_types.h>
#include <common/real_type.h>
#include <solvers/mesh/mesh_fwd.h>

#include <Kokkos_Core_fwd.hpp>

namespace alps::solver {

template<typename ValueT>
class TridiagonalSolver;

enum DiffusionCNBCType
{
  Dirichlet,
  Neumann
};

template<typename VarLoc>
class DiffusionCNEqn;

template<>
class DiffusionCNEqn<CenterPt>
{
 public:
  Mesh const&     mesh_;
  MDView<Real***> coeff_d, coeff_l, coeff_u;

  std::unique_ptr<TridiagonalSolver<Real>> solver;

  DiffusionCNEqn(Mesh const&       mesh,
                 double            alph,
                 DiffusionCNBCType top_bc_type,
                 DiffusionCNBCType bottom_bc_type);

  void initialize();

  void solve(MDView<Real***> const&               rhs,
             Kokkos::DefaultExecutionSpace const& stream) const;

  ~DiffusionCNEqn();

  double            alpha;
  DiffusionCNBCType top_bc;
  DiffusionCNBCType bottom_bc;
};

template<>
class DiffusionCNEqn<NodePt>
{
 public:
  Mesh const&     mesh_;
  MDView<Real***> coeff_d, coeff_l, coeff_u;

  std::unique_ptr<TridiagonalSolver<Real>> solver;

  DiffusionCNEqn(Mesh const& mesh, double alph);

  void initialize();

  void solve(MDView<Real***> const&               rhs,
             Kokkos::DefaultExecutionSpace const& stream) const;

  ~DiffusionCNEqn();

  double alpha;

  static constexpr DiffusionCNBCType top_bc{DiffusionCNBCType::Dirichlet};
  static constexpr DiffusionCNBCType bottom_bc{DiffusionCNBCType::Dirichlet};
};

/** @brief Set the coefficients of the CN operator for a cell-centered variable
 *
 * Set the tridiagonal matrix coefficients for the Crank-Nicolson diffusion
 * operator for a cell-centered variable $\phi$
 * \[ (I - alpha * Laplacian) \phi \]
 * in the form of:
 * \[ l_i * phi_{i-1} + d_i * phi_{i} + u_i * phi_{i+1}. \]
 */
void set_cn_operator_coefficients(MDView<Real***> const& coeff_d,
                                  MDView<Real***> const& coeff_l,
                                  MDView<Real***> const& coeff_u,
                                  double                 alpha,
                                  DiffusionCNBCType      top_bc,
                                  DiffusionCNBCType      bottom_bc,
                                  Mesh const&            mesh,
                                  CenterPt,
                                  Kokkos::DefaultExecutionSpace const& stream);

/** @brief Set the coefficients of the CN operator for a node variable
 *
 * Set the tridiagonal matrix coefficients for the Crank-Nicolson diffusion
 * operator for a node variable $\phi$
 * \[ (I - alpha * Laplacian) \phi \]
 * in the form of:
 * \[ l_i * phi_{i-1} + d_i * phi_{i} + u_i * phi_{i+1}. \]
 */
void set_cn_operator_coefficients(MDView<Real***> const& coeff_d,
                                  MDView<Real***> const& coeff_l,
                                  MDView<Real***> const& coeff_u,
                                  double                 alpha,
                                  DiffusionCNBCType      top_bc,
                                  DiffusionCNBCType      bottom_bc,
                                  Mesh const&            mesh,
                                  NodePt,
                                  Kokkos::DefaultExecutionSpace const& stream);

} // namespace alps::solver
