//
// Created by xuanx004 on 7/16/24.
//

#pragma once

#include <common/program_options/config_table.h>

#include <string_view>

namespace alps::solver {

struct RayleighDampOptionsCommon
{
  bool   enabled{false};
  double factor; /// damping strength is multiplied by this factor
  double z_zero; /// zeta coordinate where damping starts
  double z_one;  /// zeta coordinate where damping reaches its maximum

  [[nodiscard]] static RayleighDampOptionsCommon
  parse_from(ConfigTable const& config);
};

struct RayleighDampOptions : public RayleighDampOptionsCommon
{
  double u_b{0}; // background state u = b + a * z
  double u_a{0}; // background state u = b + a * z
  double v_b{0}; // background state v = b + a * z
  double v_a{0}; // background state v = b + a * z

  static constexpr std::string_view config_key{"BodyForce.RayleighDamping"};

  [[nodiscard]] static RayleighDampOptions
  parse_from(ConfigTable const& config);

  RayleighDampOptions() = default;

 private:
  RayleighDampOptions(RayleighDampOptionsCommon const& common)
    : RayleighDampOptionsCommon{common}
  {}
};

struct ScalarRayleighDampOptions : public RayleighDampOptionsCommon
{
  double b{0}; // background state y = b + a * z
  double a{0}; // background state y = b + a * z

  static constexpr std::string_view config_key{"Source.RayleighDamping"};

  [[nodiscard]] static ScalarRayleighDampOptions
  parse_from(ConfigTable const& config);

  ScalarRayleighDampOptions() = default;

 private:
  ScalarRayleighDampOptions(RayleighDampOptionsCommon const& common)
    : RayleighDampOptionsCommon{common}
  {}
};

} // namespace alps::solver
