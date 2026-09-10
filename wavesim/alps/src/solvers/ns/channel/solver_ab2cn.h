#pragma once

#include <common/base/logging_fwd.h>
#include <solvers/field/flow_field.h>
#include <solvers/ns/channel_options.h>
#include <solvers/source_terms/body_force_types.h>
#include <solvers/source_terms/scalar_source_types.h>
#include <solvers/turbulence_model/models.h>

#include <filesystem>
#include <memory>
#include <string>

namespace alps {
namespace solver {

class PressureEqn;

template<typename VarLoc>
class DiffusionCNEqn;

class ChannelFlowSolverAB2CN
{
 public:
  using SolverType = ChannelFlowSolverAB2CN;

  struct URhs
  {
    Vector3Field<Real***> Ru;
    Vector3Field<Real***> Ru_viscous;
  };

  FlowField const&     flow_field;
  ChannelSolverOptions options;

  double      dt;
  mutable int step{-1};

  /// true if the solver is at the first step of AB2
  bool steps_initialized{false};
  /// saved values of Ru for AB2 at the previous step
  Vector3Field<Real***> Ru_saved;

  std::unique_ptr<DiffusionCNEqn<CenterPt>> ueqn;
  std::unique_ptr<DiffusionCNEqn<NodePt>>   weqn;
  std::unique_ptr<PressureEqn>              peqn;

  std::vector<std::unique_ptr<BodyForce>> body_forces;

  ChannelFlowSolverAB2CN(const FlowField&     flow_field,
                         ChannelSolverOptions config);

  void initialize();

  URhs calc_explicit_rhs() const;

  void calc_uhat(URhs rhs);

  void solve_ueqn() const;

  void project() const;

  void correct() const;

  /// Validate the solution and examine the maximum divergence and the indices
  void validate() const;

  double get_time() const { return flow_field.time; }

  void set_time(double t) const { flow_field.time = t; }

  /// Save the solver state for restarting simulations
  void save(std::filesystem::path grid_file,
            std::filesystem::path data_file,
            std::filesystem::path aux_file) const;

  /// Load from files: `grid_file`, `data_file` and an auxiliary file `aux_file`
  void load(std::filesystem::path grid_file,
            std::filesystem::path data_file,
            std::filesystem::path aux_file);

  void apply_bc(const Kokkos::DefaultExecutionSpace& space,
                WhichBoundary boundary = WhichBoundary::Both) const;

  void impose_u_eqn_top_bc(const Kokkos::DefaultExecutionSpace& space) const;
  void impose_u_eqn_bottom_bc(const Kokkos::DefaultExecutionSpace& space) const;

  std::string info() const;

  ~ChannelFlowSolverAB2CN();

 private:
  Logger logger;
};

} // namespace solver
} // namespace alps
