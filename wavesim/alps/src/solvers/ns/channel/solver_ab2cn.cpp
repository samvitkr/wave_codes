#include "solver_ab2cn.h"

#include "pressure.h"
#include <common/base/logging.h>
#include <common/base/macros.h>
#include <common/container/view_utils.h>
#include <common/device/devices.h>
#include <common/kokkos_abstraction/pool_space.h>
#include <common/runtime/async_utils.h>
#include <decomp/block_partition.h>
#include <decomp/ghost_cell_exchange.h>
#include <io/hdf5.h>
#include <solvers/field/flow_field.h>
#include <solvers/ns/channel/bc.h>
#include <solvers/ns/channel/diffusion_cn.h>
#include <solvers/operators/div.h>
#include <solvers/poisson/tridiagonal_solver.h>
#include <solvers/source_terms/coriolis.h>
#include <solvers/source_terms/pressure_grad.h>
#include <solvers/source_terms/rayleigh_damp.h>
#include <spectral/spectral.h>

#include <Kokkos_Core.hpp>
#include <Kokkos_Profiling_ScopedRegion.hpp>
#include <fmt/format.h>
#include <mpipp/collectives.h>

#include <stdexcept>
#include <utility>

#include <solvers/ns/validate_max_div.h>

namespace alps::solver {

namespace {
auto VelocityBCtoCNBCType(VelocityBC const* bc) -> DiffusionCNBCType
{
  if (dynamic_cast<NoSlipWall const*>(bc) != nullptr
      || dynamic_cast<NoSlipWallVarying const*>(bc) != nullptr) {
    return DiffusionCNBCType::Dirichlet;
  }
  if (dynamic_cast<GradientWall const*>(bc) != nullptr
      || dynamic_cast<TangentialStressWall const*>(bc) != nullptr) {
    return DiffusionCNBCType::Neumann;
  }
  throw std::runtime_error("Unknown velocity BC type in DiffusionCNEqn");
}
} // anonymous namespace

ChannelFlowSolverAB2CN::ChannelFlowSolverAB2CN(const FlowField&     flow,
                                               ChannelSolverOptions config)
  : flow_field{flow}
  , options(std::move(config))
  , dt{options.dt}
  , ueqn{std::make_unique<DiffusionCNEqn<CenterPt>>(
      flow.mesh,
      dt / options.Re / 2,
      VelocityBCtoCNBCType(flow.top_bc.get()),
      VelocityBCtoCNBCType(flow.bottom_bc.get()))}
  , weqn{std::make_unique<DiffusionCNEqn<NodePt>>(flow.mesh,
                                                  dt / options.Re / 2)}
  , peqn{std::make_unique<PressureEqn>(flow_field)}
  , logger{get_logger("NS")}
{
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
  if (options.rayleigh_damp.enabled) {
    body_forces.emplace_back(
      std::make_unique<RayleighDamp<SolverType>>(*this, options.rayleigh_damp));
  }
}

void ChannelFlowSolverAB2CN::initialize()
{
  ueqn->initialize();
  weqn->initialize();
  peqn->initialize();
}

ChannelFlowSolverAB2CN::~ChannelFlowSolverAB2CN() = default;

void ChannelFlowSolverAB2CN::apply_bc(
  const Kokkos::DefaultExecutionSpace& space,
  const WhichBoundary                  boundary) const
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
      throw std::runtime_error(
        "TangentialStressWall not implemented for top boundary");
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

void ChannelFlowSolverAB2CN::validate() const
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

void ChannelFlowSolverAB2CN::save(std::filesystem::path const grid_file,
                                  std::filesystem::path const data_file,
                                  std::filesystem::path const aux_file) const
{
  Kokkos::fence();

  const auto& mesh = flow_field.mesh;

  mesh.save_to_file(grid_file);

  flow_field.save(data_file);

  {
    const auto& partition = mesh.partition();
    auto        aux_data  = io::hdf5::open_file_with_mpi(
      aux_file, HighFive::File::Overwrite, mesh.comm().raw_handle());

    if (this->Ru_saved.x.span() > 1) {
      io::hdf5::write3D_xyz(aux_data, partition, "hu", Ru_saved.x.view());
      io::hdf5::write3D_xyz(aux_data, partition, "hv", Ru_saved.y.view());
      io::hdf5::write3D_xyz(aux_data, partition, "hw", Ru_saved.z.view());
    }
  }

  if (mesh.comm().rank() == 0) {
    logger->info("Restart saved to " + grid_file.string() + ", "
                 + data_file.string() + ", " + aux_file.string());
  }

  Kokkos::fence();
}

void ChannelFlowSolverAB2CN::load(std::filesystem::path const grid_file,
                                  std::filesystem::path const data_file,
                                  std::filesystem::path const aux_file)
{
  if (!grid_file.has_filename() || !data_file.has_filename()
      || !aux_file.has_filename()) {
    throw std::invalid_argument("input files ('" + grid_file.string() + "' or '"
                                + data_file.string() + "' or '"
                                + aux_file.string() + "') are not valid paths");
  }

  Kokkos::fence();

  const auto& mesh = flow_field.mesh;
  if (mesh.comm().rank() == 0) {
    logger->info("Loading from " + grid_file.string() + ", "
                 + data_file.string() + ", " + aux_file.string());
  }

  mesh.read_from_file(grid_file);

  const auto& partition = mesh.partition();
  auto        nz        = local_end(flow_field.u.x, 2);

  if (std::filesystem::exists(data_file)) {
    auto h5file = io::hdf5::open_file_with_mpi(
      data_file.string(), HighFive::File::ReadOnly, mesh.comm().raw_handle());

    auto field = h5file.getDataSet("time");
    field.read(flow_field.time);

    io::hdf5::read3D_xyz(
      h5file,
      partition,
      "u",
      subview(flow_field.u.x, ALL, ALL, index_range(0, nz)).view());

    io::hdf5::read3D_xyz(
      h5file,
      partition,
      "v",
      subview(flow_field.u.y, ALL, ALL, index_range(0, nz)).view());

    io::hdf5::read3D_xyz(
      h5file,
      partition,
      "w",
      subview(flow_field.u.z, ALL, ALL, index_range(0, nz)).view());

    io::hdf5::read3D_xyz(
      h5file,
      partition,
      "pp",
      subview(flow_field.pp, ALL, ALL, index_range(0, nz)).view());
  } else {
    throw std::runtime_error("data file '" + data_file.string()
                             + "' does not exist");
  }

  update_halo_z(partition, flow_field.u.x, 1);
  update_halo_z(partition, flow_field.u.y, 1);
  update_halo_z(partition, flow_field.u.z, 1);
  update_halo_z(partition, flow_field.pp, 1);

  if (std::filesystem::exists(aux_file)) {
    auto aux_data = io::hdf5::open_file_with_mpi(
      aux_file, HighFive::File::ReadOnly, mesh.comm().raw_handle());

    if (aux_data.exist("/hu") && aux_data.exist("/hv")
        && aux_data.exist("/hw")) {
      this->steps_initialized = true;
      Ru_saved                = Vector3Field<Real***, default_memory_pool>(
        Kokkos::view_alloc("Ru", Kokkos::WithoutInitializing), mesh.extents());
      io::hdf5::read3D_xyz(aux_data, partition, "hu", Ru_saved.x.view());
      io::hdf5::read3D_xyz(aux_data, partition, "hv", Ru_saved.y.view());
      io::hdf5::read3D_xyz(aux_data, partition, "hw", Ru_saved.z.view());
    } else {
      if (mesh.comm().rank() == 0) {
        logger->warn("hu, hv or hw not found in restart data");
      }
    }
  }

  Kokkos::fence();
}

std::string ChannelFlowSolverAB2CN::info() const
{
  std::string os{"Channel flow solver (Cartesian grid, AB2-CN scheme)\n"};

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

  return os;
}

} // namespace alps::solver
