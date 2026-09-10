#include "radiant_heating.h"

#include <common/container/view_types.h>
#include <common/kokkos_abstraction/exec_policy.h>
#include <common/kokkos_abstraction/pool_space.h>
#include <decomp/mdcomm.h>
#include <solvers/ns/channel/solver.h>

#include <fmt/ranges.h>

RadiantHeating::RadiantHeating(alps::solver::ChannelFlowSolverAB2& flow_solver,
                               Real                                R0_,
                               std::vector<Real>                   Ai_,
                               std::vector<Real>                   Ki_)
  : mesh{&flow_solver.flow_field.mesh}
  , R0{R0_}
  , Ai{std::move(Ai_)}
  , Ki{std::move(Ki_)}
  , En(alps::HaloView<Real*, alps::default_memory_pool>(
      "En",
      {mesh->zw.begin(0), mesh->zw.end(0) - 1}))
{
  if (Ai.size() != Ki.size()) {
    throw std::runtime_error("Ai and Ki must have the same size");
  }
  if (Ai.size() == 0) {
    throw std::runtime_error("Ai and Ki must have at least one element");
  }
}

RadiantHeating
RadiantHeating::parse_from_config(alps::ConfigTable const&            config,
                                  alps::solver::ChannelFlowSolverAB2& solver)
{
  auto const R0 = config.get_value<double>("R0");
  auto const Ai = config.get_value<std::vector<double>>("Ai");
  auto const Ki = config.get_value<std::vector<double>>("Ki");

  if (Ai.size() != Ki.size()) {
    throw std::runtime_error("Ai and Ki must have the same size");
  }
  if (Ai.size() == 0) {
    throw std::runtime_error("Ai and Ki must have at least one element");
  }

  std::vector<Real> Ai_real(Ai.begin(), Ai.end());
  std::vector<Real> Ki_real(Ki.begin(), Ki.end());
  return {solver, static_cast<Real>(R0), Ai_real, Ki_real};
}

std::string RadiantHeating::info() const
{
  return fmt::format("RadiantHeating: R0 = {}, Ai = ({}), Ki = ({})",
                     R0,
                     fmt::join(Ai, ","),
                     fmt::join(Ki, ","));
}

struct AddRadiantHeating
{
  using Real = alps::Real;

  alps::MDView<Real***>       fc;
  alps::HaloView<Real const*> En;
  alps::HaloView<Real const*> zw;

  Real hbar;
  int  nx;
  int  ny;
  int  z_begin;
  int  z_end;

  AddRadiantHeating(alps::MDView<Real***> const&       fc_,
                    alps::HaloView<Real const*> const& En_,
                    alps::solver::Mesh const&          mesh)
    : fc(fc_)
    , En(En_)
    , zw(mesh.zw)
    , hbar(mesh.hbar)
    , nx(mesh.extent(0))
    , ny(mesh.extent(1))
    , z_begin(mesh.comm().is_first(2) ? 1 : 0)
    , z_end(mesh.comm().is_last(2) ? mesh.extent(2) - 1 : mesh.extent(2))
  {}

  KOKKOS_FUNCTION void
  operator()(alps::GridPolicy<>::member_type const& team) const
  {
    auto const k  = team.league_rank() / ny + z_begin;
    auto const j  = team.league_rank() % ny;
    Real       dE = (En(k) - En(k - 1)) / ((zw(k) - zw(k - 1)) * hbar);
    Kokkos::parallel_for(Kokkos::TeamVectorRange(team, nx),
                         [&](int i) { fc(i, j, k) += dE; });
  }

  void run(Kokkos::DefaultExecutionSpace const& space) const
  {
    alps::GridPolicy<> policy(space, ny * (z_end - z_begin), Kokkos::AUTO());
    Kokkos::parallel_for("add radiant source", policy, *this);
  }
};

void RadiantHeating::update_En(Kokkos::DefaultExecutionSpace const& space) const
{
  auto const En_host =
    Kokkos::create_mirror_view(Kokkos::WithoutInitializing,
                               alps::PoolSpace<Kokkos::SharedHostPinnedSpace>(),
                               En);
  for (auto k = En.begin(0); k < En.end(0); ++k) {
    Real sum{0};
    for (size_t ii = 0; ii < Ai.size(); ++ii) {
      sum += Ai[ii] * std::exp(Ki[ii] * (mesh->zw_h(k) - 1) * mesh->hbar);
    }
    En_host(k) = sum * R0;
  }
  Kokkos::deep_copy(space, En.view(), En_host.view());
}

void RadiantHeating::add_source(
  Kokkos::View<Real***, Kokkos::LayoutLeft> const& fc,
  Kokkos::DefaultExecutionSpace const&             space) const
{
  // Update En at each invokcation because it should be relatively cheap
  // also at the construction of the RadiantHeating object, the mesh may not be
  // initialized, so this guarantees that En is correctly computed
  update_En(space);
  AddRadiantHeating functor(fc, En, *mesh);

  functor.run(space);
  space.fence();
}
