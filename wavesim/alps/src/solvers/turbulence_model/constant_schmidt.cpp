//
// Created by xuanx004 on 7/14/24.
//

#include "constant_schmidt.h"

#include <common/kokkos_abstraction/exec_policy.h>
#include <common/runtime/async_utils.h>

#include <fmt/format.h>

namespace alps::solver {

ConstantTurbulentScOptions
ConstantTurbulentScOptions::parse_from(ConfigTable const& config)
{
  ConstantTurbulentScOptions options{};

  constexpr std::string_view sct_field = "Sct";
  constexpr std::string_view prt_field = "Prt";
  if (!config.contains(sct_field) && !config.contains(prt_field)) {
    throw std::runtime_error(
      "ConstantTurbulentScOptions: missing 'Sct' or 'Prt' in config");
  }
  if (config.contains(sct_field) && config.contains(prt_field)) {
    throw std::runtime_error(
      "ConstantTurbulentScOptions: both 'Sct' and 'Prt' in config");
  }

  if (config.contains(sct_field)) {
    options.Sct = config.get_value_or(sct_field, (double)options.Sct);
  } else {
    options.Sct = config.get_value_or(prt_field, (double)options.Sct);
  }

  if (options.Sct <= 0) {
    throw std::runtime_error(
      "ConstantTurbulentScOptions: 'Sct' or 'Prt' must be > 0");
  }

  return options;
}

ConstantTurbulentSc::ConstantTurbulentSc(ConstantTurbulentScOptions options)
  : options_{std::move(options)}
{}

std::string ConstantTurbulentSc::show() const
{
  return fmt::format("{} (Sc_t(Pr_t) = {})", name, options_.Sct);
}

void calculate_eddy_diffusivity(HaloView<Real***> const&   nuD,
                                ConstantTurbulentSc const& sgs_model,
                                FlowField const&           flow)
{
  auto const& nut  = flow.nu_t;
  auto const& mesh = flow.mesh;

  auto const policy =
    LoopPolicy<3>(get_next_stream(),
                  {0, 0, 0},
                  {mesh.extent(0), mesh.extent(1), mesh.extent(2)});
  auto const Sct = (Real)sgs_model.options_.Sct;
  Kokkos::parallel_for(
    "nuD const", policy, KOKKOS_LAMBDA(int i, int j, int k) {
      nuD(i, j, k) = nut(i, j, k) * Sct;
    });

  policy.space().fence();
}

} // namespace alps::solver
