#include "channel_options.h"

#include <solvers/turbulence_model/models_options.h>

#include <string>

namespace alps::solver {

ChannelSolverOptions parse_sgs_model(ChannelSolverOptions const& orig_options,
                                     ConfigTable const&          config)
{
  if (!config.contains("LES") || !config.contains("LES.SGS")) {
    return orig_options;
  }

  ChannelSolverOptions options = orig_options;

  const auto sgs_table     = config.extract_table("LES.SGS");
  options.turbulence_model = parse_sgs_model_options(sgs_table);

  options.C0_update_frequency =
    sgs_table.get_value_or("Cs_update_frequency", 1);

  return options;
}

ChannelSolverOptions parse_wall_model(ChannelSolverOptions const& orig_options,
                                      ConfigTable const&          config)
{
  if (!config.contains("LES") || !config.contains("LES.Wallmodel")
      || !config.contains("LES.Wallmodel.Bottom")) {
    return orig_options;
  }

  ChannelSolverOptions options = orig_options;

  const auto wallmodel_table = config.extract_table("LES.Wallmodel.Bottom");
  options.bottom_wall_model  = WallLayerModel::parse_options(wallmodel_table);

  return options;
}

ChannelSolverOptions parse_body_forces(ChannelSolverOptions const& orig_options,
                                       ConfigTable const&          config)
{
  ChannelSolverOptions options = orig_options;

  options.pressure_grad = PressureGradOptions::parse_from(config);

  options.coriolis = CoriolisOptions::parse_from(config);

  options.boussinesq = BoussinesqOptions::parse_from(config);

  options.rayleigh_damp = RayleighDampOptions::parse_from(config);

  return options;
}

ChannelSolverOptions ChannelSolverOptions::parse_from(ConfigTable const& config)
{
  ChannelSolverOptions options{};
  options.Re      = config.get_value<double>("Re");
  options.n_steps = config.get_value<int>("NStep");
  options.max_time =
    config.get_value_or<double>("MaxTime", std::numeric_limits<double>::max());
  options.dt = config.get_value<double>("dt");

  options = parse_body_forces(options, config);
  options = parse_sgs_model(options, config);
  options = parse_wall_model(options, config);

  if (config.contains("scalars")) {
    auto const scalar_configs = config.extract_array_of_tables("scalars");
    if (scalar_configs.size() >= 10) {
      throw std::runtime_error("Too many scalars specified");
    }
    for (std::size_t scalar_i = 0; scalar_i < scalar_configs.size();
         ++scalar_i) {
      auto const new_scalar_options =
        ScalarOptions::parse_from(scalar_configs[scalar_i], scalar_i);
      options.scalars.push_back(new_scalar_options);
    }
  }

  return options;
}

void ScalarScOrPr::setSc(double const Sc)
{
  value_ = Sc;
  label_ = "Sc";
}
void ScalarScOrPr::setPr(double Pr)
{
  value_ = Pr;
  label_ = "Pr";
}

ScalarScOrPr::ScalarScOrPr() = default;

ScalarScOrPr::ScalarScOrPr(double value, std::string label)
  : label_{std::move(label)}
  , value_{value}
{}

ScalarScOrPr ScalarScOrPr::parse_from(ConfigTable const& config)
{

  if (config.contains("Sc") && config.contains("Pr")) {
    throw std::runtime_error(
      "Cannot specify both Schmidt and Prandtl numbers for scalar");
  }
  if (!config.contains("Sc") && !config.contains("Pr")) {
    throw std::runtime_error(
      "Must specify either Schmidt (Sc) or Prandtl (Pr) number for scalar");
  }

  ScalarScOrPr result = config.contains("Sc")
                        ? ScalarScOrPr{config.get_value<double>("Sc"), "Sc"}
                        : ScalarScOrPr{config.get_value<double>("Pr"), "Pr"};

  if (result.value() < 0) {
    throw std::runtime_error("Schmidt or Prandtl number must be positive");
  }

  return result;
}

ScalarOptions ScalarOptions::parse_from(ConfigTable const& config,
                                        std::size_t        scalar_i)
{
  ScalarOptions options{};

  options.unique_id = scalar_i;
  options.label =
    config.get_value_or<std::string>("label", "c" + std::to_string(scalar_i));

  options.ScPr = ScalarScOrPr::parse_from(config);

  // SGS model
  if (config.contains("LES") && config.contains("LES.SGS")) {
    const auto sgs_table = config.extract_table("LES.SGS");
    options.sgs_options  = parse_scalar_sgs_model_options(sgs_table);
  }

  // Source terms
  options.rayleigh_damp = ScalarRayleighDampOptions::parse_from(config);

  return options;
}
} // namespace alps::solver
