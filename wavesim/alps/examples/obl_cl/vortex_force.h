//
// Created by xuanx004 on 12/27/23.
//

#pragma once

#include <common/program_options/config_table.h>
#include <solvers/ns/channel/solver.h>
#include <solvers/source_terms/body_force_types.h>

#include <functional>

struct StokesDriftMonochromaticWave
{
  KOKKOS_FORCEINLINE_FUNCTION alps::Real Us(alps::Real z) const noexcept
  {
    return static_cast<alps::Real>(Us0 * Kokkos::exp(2 * k_wave * (double)z));
  }

  KOKKOS_FORCEINLINE_FUNCTION alps::Real Vs(alps::Real /*z*/) const noexcept
  {
    return 0;
  }

  KOKKOS_FORCEINLINE_FUNCTION alps::Real Ws(alps::Real /*z*/) const noexcept
  {
    return 0;
  }

  double k_wave;
  double Us0;
};

class VortexForce : public alps::solver::BodyForce
{
 public:
  VortexForce(const alps::solver::ChannelFlowSolverAB2& flow_solver,
              StokesDriftMonochromaticWave              Us_model);

  std::string info() const override;

  void add_forces(const alps::Vector3Field<alps::Real***>& u,
                  const Kokkos::DefaultExecutionSpace& space) const override;

  alps::solver::ChannelFlowSolverAB2 const& solver;
  StokesDriftMonochromaticWave              stokes_drift;
};

StokesDriftMonochromaticWave
create_stokes_drift_from_config(alps::ConfigTable const& config);
