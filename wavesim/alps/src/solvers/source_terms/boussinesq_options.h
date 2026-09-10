//
// Created by xuanx004 on 7/22/24.
//

#pragma once

#include <common/program_options/config_table.h>

#include <limits>
#include <string>
#include <string_view>

namespace alps::solver {

struct BoussinesqOptions
{
  bool enabled{false};

  // Scalar field id or label for computing the density
  int         scalar_id{-1};
  std::string scalar_label;

  double Ri{
    0}; // Richardson coefficient from scalar variation to buoyancy forcing
  double ref_scalar{
    std::numeric_limits<double>::max()}; // reference scalar value

  static constexpr std::string_view config_key{"BodyForce.Boussinesq"};

  [[nodiscard]] static BoussinesqOptions parse_from(ConfigTable const& config);
};

} // namespace alps::solver
