//
// Created by xuananqing on 7/4/23.
//

#pragma once

#include <common/base/logging_fwd.h>
#include <solvers/field/free_surface_flow_field.h>
#include <solvers/ns/free_surface/options.h>
#include <solvers/source_terms/body_force_types.h>
#include <solvers/source_terms/scalar_source_types.h>
#include <solvers/turbulence_model/models.h>

#include <filesystem>
#include <memory>
#include <variant>

namespace alps::solver {

template<typename MeshType>
class PressureCurvilinearEqn;

class FreeSurfaceSolver
{
 public:
  using SolverType = FreeSurfaceSolver;

  const FreeSurfaceFlowField& flow_field;
  FreeSurfaceSolverOptions    options;

  double      dt;
  mutable int step{-1};

  /// true if the solver is at the first step of AB2
  bool steps_initialized{false};
  /// saved values of Ru for AB2 at the previous step
  Vector3Field<Real***> Ru_saved;

  std::unique_ptr<PressureCurvilinearEqn<TopWaveMesh>> peqn;

  // supported SGS models
  using SGSModelVariants = std::variant<std::monostate,
                                        ConstantSmagorinsky,
                                        DynamicSmagorinsky,
                                        AnisotropicMinimumDissipation>;
  SGSModelVariants turbulence_model;

  std::vector<std::unique_ptr<BodyForce>> body_forces;

  // supported scalar SGS models
  using ScalarSGSModelVariants =
    std::variant<std::monostate,
                 ConstantTurbulentSc,
                 DynamicSmagorinskyScalar,
                 AnisotropicMinimumDissipationScalar>;
  std::vector<ScalarSGSModelVariants> scalar_sgs_models;

  std::vector<std::vector<std::unique_ptr<ScalarSource>>> scalar_sources;

  std::vector<MDView<Real***>> Rc_saved;

  FreeSurfaceSolver(const FreeSurfaceFlowField& flow_field,
                    FreeSurfaceSolverOptions    config);

  void initialize() const;

  Vector3Field<Real***> calc_explicit_rhs() const;

  void calc_uhat(Vector3Field<Real***> Ru, int rk_stage);

  void project(int rk_stage) const;

  void correct(int rk_stage) const;

  void validate() const;

  void advance_scalars_ab2();

  void scale_scalars_by_J() const;

  double get_time() const { return flow_field.time; }

  void set_time(double t) const { flow_field.time = t; }

  void save(std::filesystem::path grid_file,
            std::filesystem::path data_file,
            std::filesystem::path aux_file) const; /// save restart file

  void load(std::filesystem::path grid_file,
            std::filesystem::path data_file,
            std::filesystem::path aux_file); /// load restart file

  void calc_uhat_rk2(Vector3Field<Real***> Ru, int stage);

  void apply_bottom_bc(const Kokkos::DefaultExecutionSpace& space) const;

  void apply_scalar_bc(std::size_t                          scalar_index,
                       const Kokkos::DefaultExecutionSpace& space,
                       WhichBoundary boundary = WhichBoundary::Both) const;

  void
  set_boundary_scalar_flux(Vector3Field<Real***> const& fluxes,
                           std::size_t                  scalar_index,
                           WhichBoundary boundary = WhichBoundary::Both) const;

  ~FreeSurfaceSolver();

 private:
  SGSModelVariants create_sgs_model(SGSModelOptionVariants const& opt) const;

  ScalarSGSModelVariants
  create_scalar_sgs_model(ScalarSGSModelOptionVariants const& opt);

  mutable MDView<float***> p2;

  Logger logger;
};

std::ostream& operator<<(std::ostream& os, FreeSurfaceSolver const& solver);
} // namespace alps::solver
