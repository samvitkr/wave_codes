#pragma once

#include <solvers/field/bc_types.h>
#include <solvers/field/prescribed_waves.h>

namespace alps::solver {

class FlowOverWaveField;

// clang-format off
/**
 * @brief Class describing a monochromatic wave boundary condition
 *
 * This BC is specified in the config file as:
 * ```toml
 * BC.bottom = { type = "monochromaticwavewall", wave_n = 1, wave_ak = 1.0, wave_c = 1.0 }
 * ```
 * `wave_n` is the number of wave periods in the x-direction, `wave_ak` is the
 * wave steepness, and `wave_c` is the phase speed.
 */
// clang-format on
struct MonochromaticWaveWall : public solver::NoSlipWallVarying
{

  using base_t = solver::NoSlipWallVarying;

  MonochromaticWave                params;
  solver::FlowOverWaveField const& field;

  MonochromaticWaveWall(solver::FlowOverWaveField const& field_,
                        MonochromaticWave                parameters);

  void load_from(const ConfigTable& config) override;

  std::string info() const override;

  std::unique_ptr<solver::VelocityBC> get_updated_bc(double time) override;
};

} // namespace alps::solver
