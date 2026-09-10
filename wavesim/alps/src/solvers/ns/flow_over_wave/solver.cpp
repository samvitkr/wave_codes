//
// Created by xuanx004 on 10/4/22.
//

#include "solver.h"

#include <common/base/logging.h>
#include <common/base/macros.h>
#include <common/device/devices.h>
#include <common/kokkos_abstraction/pool_space.h>
#include <common/runtime/async_utils.h>
#include <decomp/block_partition.h>
#include <solvers/field/bc_types.h>
#include <solvers/ns/curvilinear_common/diffusion_cn_zeta_eqn.h>
#include <solvers/ns/curvilinear_common/pressure.h>
#include <solvers/operators/div_curvilinear.h>
#include <solvers/source_terms/boussinesq.h>
#include <solvers/source_terms/coriolis.h>
#include <solvers/source_terms/pressure_grad.h>
#include <solvers/source_terms/rayleigh_damp.h>

#include <fmt/format.h>
#include <mpipp/collectives.h>

#include <solvers/ns/validate_max_div.h>

namespace alps::solver {

namespace {
DiffusionCNBCType velocity_bc_to_cn_bc_type(VelocityBC const* bc)
{
  if (dynamic_cast<NoSlipWall const*>(bc) != nullptr) {
    return DiffusionCNBCType::Dirichlet;
  }
  if (dynamic_cast<NoSlipWallVarying const*>(bc) != nullptr) {
    return DiffusionCNBCType::Dirichlet;
  }
  if (dynamic_cast<GradientWall const*>(bc) != nullptr) {
    return DiffusionCNBCType::Neumann;
  }
  if (dynamic_cast<TangentialStressWall const*>(bc) != nullptr) {
    return DiffusionCNBCType::Neumann;
  }

  throw std::runtime_error(
    "Unsupported velocity BC type for CN diffusion equation");
}
} // namespace

FlowOverWaveSolver::FlowOverWaveSolver(const FlowOverWaveField&  flow,
                                       FlowOverWaveSolverOptions config)
  : flow_field{flow}
  , options(std::move(config))
  , dt{options.dt}
  , peqn{std::make_unique<PressureCurvilinearEqn<BottomWaveMesh>>(flow.mesh)}
  , logger{get_logger("NS")}
{
  if (options.integrator == "ab2cn") {
    DiffusionCNZetaSolverOptions solver_options;
    if (options.DiffusionCNZetaAlgorithm == "tdma") {
      solver_options.algorithm = DiffusionCNZetaAlgorithm::TDMA;
    } else if (options.DiffusionCNZetaAlgorithm == "pcr") {
      solver_options.algorithm = DiffusionCNZetaAlgorithm::PCR;
    } else if (options.DiffusionCNZetaAlgorithm == "jacobi") {
      solver_options.algorithm = DiffusionCNZetaAlgorithm::Jacobi;
    } else {
      throw std::invalid_argument("Unknown DiffusionCNZetaAlgorithm: "
                                  + options.DiffusionCNZetaAlgorithm);
    }
    solver_options.jacobi_max_iters = options.DiffusionCNZetaJacobiMaxIters;
    solver_options.jacobi_abs_tol   = options.DiffusionCNZetaJacobiAbsTol;
    solver_options.jacobi_rel_tol   = options.DiffusionCNZetaJacobiRelTol;

    ueqn = std::make_unique<DiffusionCNZetaEqn<CenterPt>>(
      dt / options.Re / 2,
      velocity_bc_to_cn_bc_type(flow.top_bc.get()),
      velocity_bc_to_cn_bc_type(flow.bottom_bc.get()),
      solver_options);
    weqn =
      std::make_unique<DiffusionCNZetaEqn<NodePt>>(dt / options.Re / 2,
                                                   DiffusionCNBCType::Dirichlet,
                                                   DiffusionCNBCType::Dirichlet,
                                                   solver_options);
  }

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
      {flow.pp.begin(0), flow.pp.begin(1), flow.pp.begin(2)});
  }

  scalar_sgs_models.resize(flow_field.scalars.size());
  for (std::size_t is = 0; is < flow_field.scalars.size(); ++is) {
    auto const& scalar       = flow_field.scalars.at(is);
    auto const& so_sgs       = options.scalars.at(is).sgs_options;
    scalar_sgs_models.at(is) = create_scalar_sgs_model(so_sgs);

    if (!std::holds_alternative<std::monostate>(scalar_sgs_models.at(is))) {
      scalar.nuD = HaloView<Real***, default_memory_pool>(
        Kokkos::view_alloc(fmt::format("nuD{}", is),
                           Kokkos::WithoutInitializing),
        scalar.array.layout(),
        {scalar.array.begin(0), scalar.array.begin(1), scalar.array.begin(2)});
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

FlowOverWaveSolver::SGSModelVariants
FlowOverWaveSolver::create_sgs_model(SGSModelOptionVariants const& opt) const
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

FlowOverWaveSolver::ScalarSGSModelVariants
FlowOverWaveSolver::create_scalar_sgs_model(
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
          "anisotropic minimum dissipation SGS model for velocities");
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

void FlowOverWaveSolver::initialize() const
{
  auto stream = get_next_stream();
  flow_field.mesh.update_metric_coefficients(flow_field.eta, stream);
  stream.fence();
  if (options.integrator == "ab2cn") {
    ueqn->initialize();
    weqn->initialize();
  }
  peqn->initialize();
}

void FlowOverWaveSolver::project() const
{
  Kokkos::Profiling::pushRegion("Project U: div");
  auto div_u = div(flow_field.u, flow_field.mesh, CenterPt());
  Kokkos::Profiling::popRegion();

  Kokkos::Profiling::pushRegion("Project U: poisson");
  PressureCurvilinearEqnOptions const opts{Real(options.PressureEqnAbsTol),
                                           Real(options.PressureEqnRelTol),
                                           options.PressureEqnIterations};
  peqn->solve(flow_field.pp, div_u, (Real)dt, opts);

  ALPS_CHECK_LAST_DEVICE_ERROR();
  Kokkos::Profiling::popRegion();
}

void FlowOverWaveSolver::validate() const
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

std::ostream& operator<<(std::ostream& os, FlowOverWaveSolver const& solver)
{
  const auto& options = solver.options;
  const auto& mesh    = solver.flow_field.mesh;
  os << "Flow over wave solver (Curvilinear grid)\n";
  const auto kx0 = mesh.pex;
  const auto ky0 = mesh.pey;
  const auto Lz  = mesh.hbar;
  os << fmt::format("Nx x Ny x Nz = {}x{}x{}\n",
                    mesh.global_extent(0),
                    mesh.global_extent(1),
                    mesh.global_extent(2));
  os << fmt::format("Lx x Ly x Lz = {}x{}x{}\n",
                    2 * Kokkos::numbers::pi_v<Real> / kx0,
                    2 * Kokkos::numbers::pi_v<Real> / ky0,
                    Lz);
  os << fmt::format("Re = {}\n", options.Re);
  os << fmt::format("Integrator = {}\n", options.integrator);
  os << "Bottom BC: " << solver.flow_field.bottom_bc->info()
     << "; Top BC: " << solver.flow_field.top_bc->info() << "\n";
  for (const auto& body_force : solver.body_forces) {
    os << body_force->info() << "\n";
  }

  os << "Turbulence model: ";
  std::visit(
    [&](auto&& arg) {
      using ModelT = std::decay_t<decltype(arg)>;
      if constexpr (!std::is_same_v<ModelT, std::monostate>) {
        os << arg.show();
      } else {
        os << "none";
      }
    },
    solver.turbulence_model);
  os << "\n";

  os << fmt::format(
    "Pressure iteration: abs_tol = {}, rel_tol = {}, max_itrs = {}",
    options.PressureEqnAbsTol,
    options.PressureEqnRelTol,
    options.PressureEqnIterations);

  auto const& scalars        = solver.flow_field.scalars;
  auto const& scalar_options = options.scalars;
  for (std::size_t is = 0; is < scalars.size(); ++is) {
    auto const& scalar  = scalars.at(is);
    auto const& so      = scalar_options.at(is);
    auto const& sources = solver.scalar_sources.at(is);
    os << "\n";
    os << fmt::format("Scalar {}: \"{}\", {} = {}\n",
                      is,
                      scalar.label(),
                      so.ScPr.label(),
                      so.ScPr.value());
    os << "          Bottom BC: " << scalar.bottom_bc->info() << "; "
       << "Top BC: " << scalar.top_bc->info() << "\n";
    for (auto const& source : sources) {
      os << "          Source term: " << source->info() << "\n";
    }
    os << "          SGS model: ";
    std::visit(
      [&](auto&& arg) {
        using ModelT = std::decay_t<decltype(arg)>;
        if constexpr (!std::is_same_v<ModelT, std::monostate>) {
          os << arg.show();
        } else {
          os << "none";
        }
      },
      solver.scalar_sgs_models.at(is));
  }

  return os;
}

FlowOverWaveSolver::~FlowOverWaveSolver() = default;

} // namespace alps::solver
