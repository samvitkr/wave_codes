#include "solver_ab2.h"

#include "bc.h"
#include "pressure.h"
#include <common/base/logging.h>
#include <common/base/macros.h>
#include <common/container/view_utils.h>
#include <common/device/devices.h>
#include <common/kokkos_abstraction/pool_space.h>
#include <common/runtime/async_utils.h>
#include <decomp/block_partition.h>
#include <solvers/field/flow_field.h>
#include <solvers/operators/div.h>
#include <solvers/source_terms/boussinesq.h>
#include <solvers/source_terms/coriolis.h>
#include <solvers/source_terms/pressure_grad.h>
#include <solvers/source_terms/rayleigh_damp.h>
#include <solvers/turbulence_model/log_law_wall_model.h>
#include <solvers/turbulence_model/models.h>

#include <Kokkos_Core.hpp>
#include <fmt/format.h>
#include <mpipp/collectives.h>

#include <stdexcept>
#include <utility>
#include <variant>

#include <solvers/ns/validate_max_div.h>

namespace alps {
namespace solver {

ChannelFlowSolverAB2::ChannelFlowSolverAB2(const FlowField&     flow,
                                           ChannelSolverOptions config)
  : flow_field{flow}
  , options(std::move(config))
  , dt{options.dt}
  , peqn{std::make_unique<PressureEqn>(flow_field)}
  , logger{get_logger("NS")}
{
  Rc_saved.clear();
  Rc_saved.resize(flow_field.scalars.size());

  scalar_sources.resize(flow_field.scalars.size());

  turbulence_model = create_sgs_model(options.turbulence_model);
  if (options.C0_update_frequency < 1) options.C0_update_frequency = 1;

  // Allocate storage needed for the turbulence model
  if (!std::holds_alternative<std::monostate>(turbulence_model)) {
    flow_field.nu_t = HaloView<Real***, default_memory_pool>(
      Kokkos::view_alloc("nu_t", Kokkos::WithoutInitializing),
      flow.pp.layout(),
      {begin(flow.pp, 0), begin(flow.pp, 1), begin(flow.pp, 2)});
  }

  if (auto const* option =
        std::get_if<LogLawWallModelOptions>(&options.bottom_wall_model);
      option != nullptr) {
    wall_model_bottom = std::make_unique<LogLawWallModel>(*option);
  }

  scalar_sgs_models.resize(flow_field.scalars.size());
  for (std::size_t is = 0; is < flow_field.scalars.size(); ++is) {
    auto const& so_sgs       = options.scalars.at(is).sgs_options;
    scalar_sgs_models.at(is) = create_scalar_sgs_model(so_sgs);

    if (!std::holds_alternative<std::monostate>(scalar_sgs_models.at(is))) {
      flow_field.scalars.at(is).nuD = HaloView<Real***, default_memory_pool>(
        Kokkos::view_alloc(fmt::format("nuD{}", is),
                           Kokkos::WithoutInitializing),
        flow.pp.layout(),
        {begin(flow.pp, 0), begin(flow.pp, 1), begin(flow.pp, 2)});
    }
  }

  // Add body forces
  if (options.pressure_grad.enabled) {
    body_forces.emplace_back(
      std::make_unique<ConstantPressureGradient<SolverType>>(
        *this, options.pressure_grad));
  }
  if (options.coriolis.enabled) {
    body_forces.emplace_back(
      std::make_unique<CoriolisForce<SolverType>>(*this, options.coriolis));
  }
  if (options.boussinesq.enabled) {
    body_forces.emplace_back(
      std::make_unique<BoussinesqForce<SolverType>>(*this, options.boussinesq));
  }
  if (options.rayleigh_damp.enabled) {
    body_forces.emplace_back(
      std::make_unique<RayleighDamp<SolverType>>(*this, options.rayleigh_damp));
  }

  // Add scalar sources
  for (std::size_t is = 0; is < flow_field.scalars.size(); ++is) {
    auto const& so = options.scalars.at(is);
    if (so.rayleigh_damp.enabled) {
      scalar_sources.at(is).emplace_back(
        std::make_unique<RayleighDampScalar<SolverType>>(
          *this, is, so.rayleigh_damp));
    }
  }
}

ChannelFlowSolverAB2::SGSModelVariants
ChannelFlowSolverAB2::create_sgs_model(SGSModelOptionVariants const& opt) const
{
  auto try_match = [&](auto&& arg) -> SGSModelVariants {
    using OptT = std::decay_t<decltype(arg)>;
    // Explicitly check for all supported SGS models
    if constexpr (std::is_same_v<OptT, std::monostate>) {
      return std::monostate{};
    } else if constexpr (std::is_same_v<OptT, ConstantSmagorinskyOptions>) {
      return ConstantSmagorinsky{arg};
    } else if constexpr (std::is_same_v<OptT, DynamicSmagorinskyOptions>) {
      return DynamicSmagorinsky{arg};
    } else if constexpr (std::is_same_v<OptT,
                                        AnisotropicMinimumDissipationOptions>) {
      return AnisotropicMinimumDissipation{arg};
    } else {
      std::string_view model_name = OptT::ModelType::name;
      throw std::runtime_error("Unsupported SGS model type: "
                               + std::string(model_name));
    }
    ALPS_UNREACHABLE(std::monostate{});
  };
  return std::visit(try_match, opt);
}

ChannelFlowSolverAB2::ScalarSGSModelVariants
ChannelFlowSolverAB2::create_scalar_sgs_model(
  ScalarSGSModelOptionVariants const& opt)
{
  auto check_flow_sgs = [&] {
    if (std::holds_alternative<std::monostate>(turbulence_model)) {
      throw std::runtime_error(
        "SGS model for scalar requires a flow SGS model");
    }
  };
  auto try_match = [&](auto&& arg) -> ScalarSGSModelVariants {
    using OptT = std::decay_t<decltype(arg)>;
    // Explicitly check for all supported SGS models
    if constexpr (std::is_same_v<OptT, std::monostate>) {
      return std::monostate{};
    } else if constexpr (std::is_same_v<OptT, ConstantTurbulentScOptions>) {
      check_flow_sgs();
      return ConstantTurbulentSc{arg};
    } else if constexpr (std::is_same_v<OptT,
                                        DynamicSmagorinskyScalarOptions>) {
      auto* tu = std::get_if<DynamicSmagorinsky>(&turbulence_model);
      if (tu == nullptr) {
        throw std::runtime_error(
          "Scalar dynamic Smagorinsky SGS model should be paired a dynamic "
          "Smagorinsky SGS model for velocities");
      }
      return DynamicSmagorinskyScalar{*tu};
    } else if constexpr (std::is_same_v<
                           OptT,
                           AnisotropicMinimumDissipationScalarOptions>) {
      auto* tu = std::get_if<AnisotropicMinimumDissipation>(&turbulence_model);
      if (tu == nullptr) {
        throw std::runtime_error(
          "Scalar anisotropic minimum dissipation SGS model should be paired "
          "an anisotropic minimum dissipation SGS model for velocities");
      }
      return AnisotropicMinimumDissipationScalar{*tu};
    } else {
      std::string_view model_name = OptT::ModelType::name;
      throw std::runtime_error("Unsupported scalar SGS model type: "
                               + std::string(model_name));
    }
    ALPS_UNREACHABLE(std::monostate{});
  };
  return std::visit(try_match, opt);
}

void ChannelFlowSolverAB2::initialize() const
{
  peqn->initialize();
}

ChannelFlowSolverAB2::~ChannelFlowSolverAB2() = default;

void ChannelFlowSolverAB2::apply_bc(const Kokkos::DefaultExecutionSpace& space,
                                    const WhichBoundary boundary) const
{
  auto        stream1 = get_next_stream();
  auto const& flow    = this->flow_field;

  if (boundary == WhichBoundary::TopBC || boundary == WhichBoundary::Both) {
    auto event = get_device_event();
    enqueue(event, space);
    wait_for(event, stream1);

    bool  bc_set   = false;
    auto* base_ptr = flow.top_bc.get();
    if (auto const* bc = dynamic_cast<NoSlipWall*>(base_ptr); bc != nullptr) {
      apply_top_bc(flow, *bc, stream1);
      bc_set = true;
    }
    if (auto const* bc = dynamic_cast<GradientWall*>(base_ptr); bc != nullptr) {
      apply_top_bc(flow, *bc, stream1);
      bc_set = true;
    }
    if (auto const* bc = dynamic_cast<TangentialStressWall*>(base_ptr);
        bc != nullptr) {
      // specify gradient from stress
      using real_t = std::decay_t<decltype(flow.u.x)>::non_const_value_type;
      GradientWall bc_grad;
      bc_grad.grad_1 = real_t((double)bc->tau_1 * options.Re);
      bc_grad.grad_2 = real_t((double)bc->tau_2 * options.Re);
      apply_top_bc(flow, bc_grad, stream1);
      bc_set = true;
    }

    if (!bc_set) {
      logger->warn("No top boundary condition set.");
    }
  }

  if (boundary == WhichBoundary::BottomBC || boundary == WhichBoundary::Both) {
    bool  bc_set   = false;
    auto* base_ptr = flow.bottom_bc.get();
    if (auto const* bc = dynamic_cast<NoSlipWall*>(base_ptr); bc != nullptr) {
      apply_bottom_bc(flow, *bc, space);
      bc_set = true;
    }
    if (auto const* bc = dynamic_cast<GradientWall*>(base_ptr); bc != nullptr) {
      apply_bottom_bc(flow, *bc, space);
      bc_set = true;
    }
    if (auto const* bc = dynamic_cast<TangentialStressWall*>(base_ptr);
        bc != nullptr) {
      throw std::runtime_error(
        "TangentialStressWall not implemented for bottom boundary");
    }

    if (!bc_set) {
      logger->warn("No bottom boundary condition set");
    }
  }

  space.fence();
  stream1.fence();
}

void ChannelFlowSolverAB2::validate() const
{
  auto const [max_div, max_loc] = validate_max_div(flow_field);
  if (flow_field.mesh.comm().comm.rank() == 0) {
    logger->info("Max div: {:.3e} @ ({}, {}, {})",
                 max_div,
                 max_loc[0],
                 max_loc[1],
                 max_loc[2]);
  }
}

std::string ChannelFlowSolverAB2::info() const
{
  std::string os{"Channel flow solver (Cartesian grid)\n"};

  const auto& mesh = flow_field.mesh;
  const auto  kx0  = mesh.pex;
  const auto  ky0  = mesh.pey;
  const auto  Lz   = mesh.hbar;
  os += fmt::format("Nx x Ny x Nz = {} x {} x {}\n",
                    mesh.global_extent(0),
                    mesh.global_extent(1),
                    mesh.global_extent(2));
  os += fmt::format("Lx x Ly x Lz = {} x {} x {}\n",
                    2 * Kokkos::numbers::pi_v<Real> / kx0,
                    2 * Kokkos::numbers::pi_v<Real> / ky0,
                    Lz);
  os += fmt::format("Re = {}\n", options.Re);
  os += fmt::format("Bottom BC: {}; Top BC: {}",
                    flow_field.bottom_bc->info(),
                    flow_field.top_bc->info());
  for (const auto& body_force : body_forces) {
    os += "\n";
    os += body_force->info();
  }
  os += "\nTurbulence model: ";
  std::visit(
    [&](auto&& arg) {
      using ModelT = std::decay_t<decltype(arg)>;
      if constexpr (!std::is_same_v<ModelT, std::monostate>) {
        os += arg.show();
      } else {
        os += "none";
      }
    },
    turbulence_model);
  os += "\n";

  if (wall_model_bottom != nullptr) {
    os += fmt::format("Bottom wall model: {}", wall_model_bottom->info());
  } else {
    os += "Bottom wall model: None";
  }

  auto const& scalars        = flow_field.scalars;
  auto const& scalar_options = options.scalars;
  for (std::size_t is = 0; is < scalars.size(); ++is) {
    auto const& scalar  = scalars.at(is);
    auto const& so      = scalar_options.at(is);
    auto const& sources = scalar_sources.at(is);
    os += "\n";
    os += fmt::format("Scalar {}: \"{}\", {} = {}\n",
                      is,
                      scalar.label(),
                      so.ScPr.label(),
                      so.ScPr.value());
    os += fmt::format("          Bottom BC: {}; Top BC: {}\n",
                      scalar.bottom_bc->info(),
                      scalar.top_bc->info());
    for (auto const& source : sources) {
      os += fmt::format("          Source term: {}\n", source->info());
    }
    os += "          SGS model: ";
    std::visit(
      [&](auto&& arg) {
        using ModelT = std::decay_t<decltype(arg)>;
        if constexpr (!std::is_same_v<ModelT, std::monostate>) {
          os += arg.show();
        } else {
          os += "none";
        }
      },
      scalar_sgs_models.at(is));
  }

  return os;
}

} // namespace solver
} // namespace alps
