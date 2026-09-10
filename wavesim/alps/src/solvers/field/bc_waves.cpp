#include "bc_waves.h"

#include <common/kokkos_abstraction/exec_policy.h>
#include <solvers/field/flow_over_wave_field.h>

#include <fmt/format.h>

namespace alps::solver {

MonochromaticWaveWall::MonochromaticWaveWall(FlowOverWaveField const& field_,
                                             MonochromaticWave parameters)
  : base_t(field_.mesh.extent(0), field_.mesh.extent(1), true)
  , params{parameters}
  , field{field_}
{}

void MonochromaticWaveWall::load_from(const ConfigTable& config)
{
  params = MonochromaticWave::parse_from(config);
}

std::string MonochromaticWaveWall::info() const
{
  return fmt::format("Monochromatic Wave (nwaves={}, ak={}, c⁺={})",
                     params.n_waves,
                     params.ak,
                     params.phase_c);
}

std::unique_ptr<solver::VelocityBC>
MonochromaticWaveWall::get_updated_bc(double time)
{
  auto new_bc = std::make_unique<MonochromaticWaveWall>(field, this->params);

  const auto& bottom_u = new_bc->u_;
  const auto& bottom_v = new_bc->v_;
  const auto& bottom_w = new_bc->w_;
  const auto& eta      = new_bc->eta_;
  const auto& eta_t    = new_bc->eta_t_;

  const auto c_p = params.phase_c;
  const auto kx0 = field.mesh.pex;
  const auto amp = double(params.ak / params.n_waves / (double)kx0);
  const auto k   = (double)kx0 * params.n_waves;
  const auto Lx  = Kokkos::numbers::pi_v<double> * 2 / (double)kx0;
  const auto nx  = eta.extent_int(0);
  const auto ny  = eta.extent_int(1);

  const auto policy = alps::LoopPolicy<2>({0, 0}, {nx, ny});
  Kokkos::parallel_for(
    policy, KOKKOS_LAMBDA(int i, int j) {
      const auto x     = i * (Lx / nx);
      const auto phase = k * (x - c_p * time);
      eta(i, j)        = static_cast<Real>(amp * Kokkos::cos(phase));
      bottom_u(i, j)   = static_cast<Real>(amp * c_p * k * Kokkos::cos(phase));
      bottom_w(i, j)   = static_cast<Real>(amp * c_p * k * Kokkos::sin(phase));
      bottom_v(i, j)   = 0;
      eta_t(i, j)      = bottom_w(i, j);
    });

  policy.space().fence(); // fence to ensure all updates before return

  return new_bc;
}

} // namespace alps::solver
