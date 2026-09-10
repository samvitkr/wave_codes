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
#include <variant>

namespace alps {
namespace solver {

class PressureEqn;

class ChannelFlowSolverAB2
{
 public:
  using SolverType = ChannelFlowSolverAB2;

  FlowField const&     flow_field;
  ChannelSolverOptions options;

  double      dt;
  mutable int step{-1};

  /// true if the solver is at the first step of AB2
  bool steps_initialized{false};
  /// saved values of Ru for AB2 at the previous step
  Vector3Field<Real***> Ru_saved;

  std::unique_ptr<PressureEqn> peqn;

  // supported SGS models
  using SGSModelVariants = std::variant<std::monostate,
                                        ConstantSmagorinsky,
                                        DynamicSmagorinsky,
                                        AnisotropicMinimumDissipation>;
  SGSModelVariants turbulence_model;

  std::vector<std::unique_ptr<BodyForce>> body_forces;

  std::unique_ptr<WallLayerModel> wall_model_bottom;

  // supported scalar SGS models
  using ScalarSGSModelVariants =
    std::variant<std::monostate,
                 ConstantTurbulentSc,
                 DynamicSmagorinskyScalar,
                 AnisotropicMinimumDissipationScalar>;
  std::vector<ScalarSGSModelVariants> scalar_sgs_models;

  std::vector<std::vector<std::unique_ptr<ScalarSource>>> scalar_sources;

  std::vector<MDView<Real***>> Rc_saved;

  ChannelFlowSolverAB2(const FlowField&     flow_field,
                       ChannelSolverOptions config);

  void initialize() const;

  Vector3Field<Real***> calc_explicit_rhs() const;

  void calc_uhat(Vector3Field<Real***> Ru);

  void project() const;

  void correct() const;

  /// Validate the solution and examine the maximum divergence and the indices
  void validate() const;

  void advance_scalars_ab2();

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

  void
  set_boundary_stress_flux(Tensor33Field<Real***> const& fluxes,
                           WhichBoundary boundary = WhichBoundary::Both) const;

  void apply_scalar_bc(std::size_t                          scalar_index,
                       const Kokkos::DefaultExecutionSpace& space,
                       WhichBoundary boundary = WhichBoundary::Both) const;

  void
  set_boundary_scalar_flux(Vector3Field<Real***> const& fluxes,
                           std::size_t                  scalar_index,
                           WhichBoundary boundary = WhichBoundary::Both) const;

  std::string info() const;

  ~ChannelFlowSolverAB2();

 private:
  SGSModelVariants create_sgs_model(SGSModelOptionVariants const& opt) const;

  ScalarSGSModelVariants
  create_scalar_sgs_model(ScalarSGSModelOptionVariants const& opt);

  Logger logger;
};

} // namespace solver
} // namespace alps
