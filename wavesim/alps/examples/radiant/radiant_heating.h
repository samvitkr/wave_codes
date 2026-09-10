#pragma once

#include <common/container/view_types.h>
#include <common/program_options/config_table.h>
#include <solvers/source_terms/scalar_source_types.h>

namespace alps::solver {
class Mesh;
class ChannelFlowSolverAB2;
} // namespace alps::solver

/**
 * @class RadiantHeating
 * @brief Represents a radiant heating source term from Ohlmann & Siegel (JPO,
 * 2000).
 */
class RadiantHeating : public alps::solver::ScalarSource
{
 private:
  using Real = alps::Real;

 public:
  static RadiantHeating
  parse_from_config(alps::ConfigTable const&            config,
                    alps::solver::ChannelFlowSolverAB2& solver);

  /**
   * @brief Constructs a RadiantHeating object.
   *
   * This constructor initializes a RadiantHeating object with the given
   * parameters.
   *
   * @param flow_solver The ChannelFlowSolverAB2 object used for solving the
   * flow.
   * @param R0 The net incoming incident radiation.
   * @param Ai A vector of Real values representing the Ai values.
   * @param Ki A vector of Real values representing the Ki values.
   */
  RadiantHeating(alps::solver::ChannelFlowSolverAB2& flow_solver,
                 Real                                R0_,
                 std::vector<Real>                   Ai_,
                 std::vector<Real>                   Ki_);

  std::string info() const override;

  void add_source(Kokkos::View<Real***, Kokkos::LayoutLeft> const& fc,
                  Kokkos::DefaultExecutionSpace const& space) const override;

 private:
  void update_En(Kokkos::DefaultExecutionSpace const& space) const;

  alps::solver::Mesh const* mesh;
  Real                      R0;
  std::vector<Real>         Ai;
  std::vector<Real>         Ki;

  alps::HaloView<Real*> En;
};
