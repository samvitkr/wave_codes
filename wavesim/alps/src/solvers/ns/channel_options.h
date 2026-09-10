#pragma once

#include <common/program_options/config_table.h>
#include <solvers/source_terms/boussinesq_options.h>
#include <solvers/source_terms/coriolis_options.h>
#include <solvers/source_terms/pressure_grad_options.h>
#include <solvers/source_terms/rayleigh_damp_options.h>
#include <solvers/turbulence_model/models_options.h>
#include <solvers/turbulence_model/wall_models.h>

#include <vector>

namespace alps::solver {

class ScalarScOrPr
{
 public:
  std::string label() const { return label_; }

  void setSc(double Sc);
  void setPr(double Pr);

  auto value() const { return value_; }

  ScalarScOrPr();

  ScalarScOrPr(double value, std::string label);

  static ScalarScOrPr parse_from(ConfigTable const& config);

 private:
  std::string label_{"unset"};
  double      value_{};
};

struct ChannelSolverOptions;

struct ScalarOptions
{
  std::size_t unique_id;
  std::string label;

  ScalarScOrPr ScPr;

  ScalarSGSModelOptionVariants sgs_options{std::monostate{}};

  ScalarRayleighDampOptions rayleigh_damp;

  static ScalarOptions parse_from(ConfigTable const& config,
                                  std::size_t        scalar_i);
};

struct ChannelSolverOptions
{
  double Re{};

  int    n_steps{};
  double max_time{std::numeric_limits<double>::max()};
  double dt{};

  PressureGradOptions pressure_grad;
  CoriolisOptions     coriolis;
  BoussinesqOptions   boussinesq;
  RayleighDampOptions rayleigh_damp;

  // Turbulence model options
  SGSModelOptionVariants turbulence_model{std::monostate{}};
  // Smagorinsky type specific parameters
  int C0_update_frequency{1};

  // Wall model options
  WallLayerModelOptionVariants bottom_wall_model{std::monostate{}};

  std::vector<ScalarOptions> scalars;

  static ChannelSolverOptions parse_from(ConfigTable const& config);
};

} // namespace alps::solver
