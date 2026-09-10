#include "smooth.h"

#include <common/utils/to_lower_case.h>
#include <spectral/spectral.h>

#include <fmt/format.h>

namespace alps::solver::hos {

std::unique_ptr<Smoother> make_smoother(const alps::ConfigTable& config)
{
  auto const type = to_lower_case(config.get_value<std::string>("type"));

  if (type == "none") {
    return std::make_unique<Smoother>();
  }
  if (type == "lowpass") {
    auto factor = config.get_value_or<Real>("factor", 0.5);
    return std::make_unique<LowPassFilter>(factor);
  }
  throw std::runtime_error("Unknown smoother: " + type);
}

Smoother::~Smoother() {}

void Smoother::apply(decltype(HOSField::value)            solution,
                     Grid const&                          grid,
                     Kokkos::DefaultExecutionSpace const& space) const
{
  Kokkos::Profiling::pushRegion("apply " + name());

  apply_impl(solution, grid, space);

  Kokkos::Profiling::popRegion();
}

void Smoother::apply_impl(decltype(HOSField::value)            solution,
                          Grid const&                          grid,
                          Kokkos::DefaultExecutionSpace const& space) const
{
  // do nothing
  (void)solution; // suppress unused parameter warning
  (void)grid;     // suppress unused parameter warning
  (void)space;    // suppress unused parameter warning
}

std::string Smoother::name() const
{
  return "empty smoother";
}

std::string Smoother::info() const
{
  return "empty smoother";
}

LowPassFilter::LowPassFilter(Real factor)
  : factor_{factor}
{}

void LowPassFilter::apply_impl(decltype(HOSField::value)            solution,
                               Grid const&                          grid,
                               Kokkos::DefaultExecutionSpace const& space) const
{
  auto kc_x = static_cast<int>((int)(solution.extent_int(0) / 2) * factor_);
  auto kc_y = static_cast<int>((int)(solution.extent_int(1) / 2) * factor_);
  spectral::cutoff_xy(solution, kc_x, kc_y, grid, space);

  space.fence();
}

std::string LowPassFilter::name() const
{
  return "LowPassFilter";
}

std::string LowPassFilter::info() const
{
  return fmt::format("LowPassFilter (factor = {})", factor_);
}

} // namespace alps::solver::hos
