#include "flow_over_wave_field.h"

#include "flow_field.h"
#include <common/base/logging.h>
#include <common/container/view_utils.h>
#include <common/kokkos_abstraction/pool_space.h>
#include <common/utils/to_lower_case.h>
#include <decomp/block_partition.h>
#include <io/hdf5.h>
#include <solvers/field/bc_waves.h>
#include <solvers/mesh/curvilinear_mesh.h>

namespace alps {
namespace solver {
FlowOverWaveField::FlowOverWaveField(
  const alps::solver::BottomWaveMesh& local_mesh)
  : mesh{local_mesh}
  , u(Vector3Field<Real***, default_memory_pool>(
      Kokkos::view_alloc("u", Kokkos::WithoutInitializing),
      local_mesh.extents(),
      {0, 0, 1}))
  , pp(HaloView<Real***, default_memory_pool>(
      Kokkos::view_alloc("p", Kokkos::WithoutInitializing),
      u.x.layout(),
      {0, 0, -1}))
  , eta(MDView<Real**, default_memory_pool>(
      Kokkos::view_alloc("eta", Kokkos::WithoutInitializing),
      local_mesh.extent(0),
      local_mesh.extent(1)))
  , logger(get_logger("field"))
{
  // Zeroes and touches the variables
  Kokkos::deep_copy(u, 0);
  Kokkos::deep_copy(pp, 0);

  top_bc    = std::make_unique<GradientWall>();
  bottom_bc = std::make_unique<NoSlipWall>();
}

FlowField FlowOverWaveField::as_channel_flow_field() const
{
  return FlowField(mesh, u, pp, nu_t, scalars.shallow_copy(), time);
}

std::unique_ptr<VelocityBC>
FlowOverWaveField::parse_bc(const ConfigTable& config, std::string bc_key) const
{
  auto const bc_config = config.extract_table(bc_key);
  auto const bc_type   = bc_config.get_value<std::string>("type");
  // convert std::string bc_type to lower case
  auto const bc_type_lower = to_lower_case(bc_type);

  if (bc_type_lower == "monochromaticwavewall") {
    auto bc =
      std::make_unique<MonochromaticWaveWall>(*this, MonochromaticWave());
    bc->load_from(bc_config);
    return bc;
  }

  // call FlowField::parse_bc to handle the rest of the BCs
  return as_channel_flow_field().parse_bc(config, bc_key);
}

void FlowOverWaveField::parse_bcs_from(const ConfigTable& config,
                                       WhichBoundary      which)
{
  if (which == WhichBoundary::Both || which == WhichBoundary::BottomBC) {
    parse_bottom_bc_from(config);
  }

  if (which == WhichBoundary::Both || which == WhichBoundary::TopBC) {
    parse_top_bc_from(config);
  }
}

void FlowOverWaveField::parse_bottom_bc_from(const ConfigTable& config)
{
  auto constexpr bc_key = "BC.bottom";
  if (!config.contains(bc_key)) {
    throw std::runtime_error("No boundary condition specified for bottom");
  }
  auto new_bc = parse_bc(config, bc_key);
  set_bottom_bc(std::move(new_bc));
}

void FlowOverWaveField::parse_top_bc_from(const ConfigTable& config)
{
  auto constexpr bc_key = "BC.top";
  if (!config.contains(bc_key)) {
    throw std::runtime_error("No boundary condition specified for top");
  }
  auto new_bc = parse_bc(config, bc_key);
  set_top_bc(std::move(new_bc));
}

void FlowOverWaveField::set_bottom_bc(std::unique_ptr<VelocityBC> bc) const
{
  bottom_bc = std::move(bc);
}

void FlowOverWaveField::set_top_bc(std::unique_ptr<VelocityBC> bc) const
{
  top_bc = std::move(bc);
}

void FlowOverWaveField::save(HighFive::File& h5_file) const
{
  Kokkos::fence();

  this->as_channel_flow_field().save(h5_file);

  const auto            is_bot    = mesh.comm().is_first(2);
  const BlockPartition& partition = mesh.partition();

  // save eta and eta_t, only bottom processors write
  const auto block_shape_xy =
    is_bot ? std::vector{mesh.extent(0), mesh.extent(1)} : std::vector{0, 0};
  const auto offset_xy =
    is_bot ? std::vector{partition.offsets[0], partition.offsets[1]}
           : std::vector{0, 0};
  const std::vector total_shape_xy{partition.global_extents[0],
                                   partition.global_extents[1]};

  io::hdf5::write_blocks(
    h5_file, "eta", eta, total_shape_xy, block_shape_xy, offset_xy);
  io::hdf5::write_blocks(
    h5_file, "eta_t", mesh.et, total_shape_xy, block_shape_xy, offset_xy);

  if (mesh.comm().rank() == 0) {
    logger->info("eta, eta_t saved to {}", h5_file.getName());
  }

  Kokkos::fence();
}

void FlowOverWaveField::save(std::filesystem::path filename) const
{
  if (!filename.has_filename()) {
    throw std::invalid_argument("output filename " + filename.string()
                                + " is not valid");
  }

  if (mesh.comm().rank() == 0) {}

  HighFive::File h5file = io::hdf5::open_file_with_mpi(
    filename.string(), HighFive::File::Overwrite, mesh.comm().raw_handle());

  save(h5file);
}

void FlowOverWaveField::initialize_scalars(ConfigTable const& config)
{
  if (!config.contains("scalars")) return;

  auto scalar_configs = config.extract_array_of_tables("scalars");
  scalars.initialize(mesh, scalar_configs);
}

} // namespace solver
} // namespace alps
