//
// Created by xuanx004 on 5/26/23.
//

#pragma once

#include "wall_models.h"

#include <common/container/matrix_field.h>
#include <common/container/view_types.h>
#include <common/real_type.h>

#include <string_view>

namespace alps::solver {

class FlowField;

/// @brief Class describing a log law wall model
/** see @ref LogLawWallModelOptions for the parameters */
class LogLawWallModel : public WallLayerModel
{
 public:
  LogLawWallModel(LogLawWallModelOptions options);

  std::string info() const override;

  static constexpr std::string_view name = "log law wall model";

  LogLawWallModelOptions options_;
};

/// @brief Apply the wall shear stress to the momentum fluxes
/**
 * @param fluxes Momentum fluxes
 * @param flow Flow field
 * @param wall_model
 */
void apply_bottom_wall_shear_flux(Tensor33Field<Real***> const& fluxes,
                                  FlowField const&              flow,
                                  LogLawWallModel const&        wall_model);

/// @brief Calculate the magnitude of the wall shear stress using the given wall
/// model
/**
 * @param shear_flux_mag Result of the wall shear stress magnitude
 * @param u_relative Off-wall velocity relative to the wall velocity
 * @param delta_z Distance from the wall
 * @param wall_model
 */
void calculate_wall_shear_flux(MDView<Real**> const&           shear_flux_mag,
                               MDView<Real const** [2]> const& u_relative,
                               Real                            delta_z,
                               LogLawWallModel const&          wall_model);

} // namespace alps::solver
