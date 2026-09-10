#include "flow_field.h"

#include <common/base/logging.h>
#include <common/container/view_utils.h>
#include <common/kokkos_abstraction/pool_space.h>
#include <common/program_options/config_table.h>
#include <common/utils/to_lower_case.h>
#include <io/hdf5.h>
#include <solvers/field/bc_types.h>
#include <solvers/mesh/mesh.h>

#include <algorithm>

namespace alps {
namespace solver {
FlowField::FlowField(const Mesh& local_mesh)
  : mesh{local_mesh}
  , u(Vector3Field<Real***, default_memory_pool>(
      Kokkos::view_alloc("u", Kokkos::WithoutInitializing),
      local_mesh.extents(),
      {0, 0, 1}))
  , pp(HaloView<Real***, default_memory_pool>(
      Kokkos::view_alloc("p", Kokkos::WithoutInitializing),
      u.x.layout(),
      {0, 0, -1}))
  , logger(get_logger("field"))
{
  // Zeroes and touches the variables
  Kokkos::deep_copy(u, 0);
  Kokkos::deep_copy(pp, 0);

  top_bc    = std::make_unique<GradientWall>();
  bottom_bc = std::make_unique<NoSlipWall>();
}

FlowField::FlowField(Mesh const&                  mesh_,
                     Vector3Field<Real***> const& u_,
                     HaloView<Real***> const&     pp_,
                     HaloView<Real***> const&     nu_t_,
                     ScalarFields                 scalars_,
                     double                       time_)
  : mesh{mesh_}
  , u{u_}
  , pp{pp_}
  , nu_t{nu_t_}
  , scalars{std::move(scalars_)}
  , time{time_}
  , logger{get_logger("field")}
{}

std::unique_ptr<VelocityBC> FlowField::parse_bc(const ConfigTable& config,
                                                std::string        bc_key) const
{
  auto const bc_config = config.extract_table(bc_key);
  auto const bc_type   = bc_config.get_value<std::string>("type");
  // convert std::string bc_type to lower case
  auto const bc_type_lower = to_lower_case(bc_type);

  if (bc_type_lower == "gradientwall") {
    auto bc = std::make_unique<GradientWall>();
    bc->load_from(bc_config);
    return bc;
  }
  if (bc_type_lower == "noslipwall") {
    auto bc = std::make_unique<NoSlipWall>();
    bc->load_from(bc_config);
    return bc;
  }
  if (bc_type_lower == "tangentialstresswall") {
    auto bc = std::make_unique<TangentialStressWall>();
    bc->load_from(bc_config);
    return bc;
  }
  if (bc_type_lower == "custom" or bc_type_lower == "coded") {
    // do nothing
    // Custom BC is set externally by the user, usually by explicitly coding the
    // BC in the main program
    return nullptr;
  }

  throw std::runtime_error("Cannot handle boundary condition of type " + bc_type
                           + " at key " + bc_key);
  ALPS_UNREACHABLE(nullptr);
}

void FlowField::parse_bcs_from(const ConfigTable& config, WhichBoundary which)
{
  if (which == WhichBoundary::BottomBC || which == WhichBoundary::Both) {
    parse_bottom_bc_from(config);
  }

  if (which == WhichBoundary::TopBC || which == WhichBoundary::Both) {
    parse_top_bc_from(config);
  }
}

void FlowField::parse_bottom_bc_from(const ConfigTable& config)
{
  auto constexpr bc_key = "BC.bottom";
  if (!config.contains(bc_key)) {
    throw std::runtime_error("No boundary condition specified for bottom");
  }
  auto new_bc = parse_bc(config, bc_key);
  set_bottom_bc(std::move(new_bc));
}

void FlowField::parse_top_bc_from(const ConfigTable& config)
{
  auto constexpr bc_key = "BC.top";
  if (!config.contains(bc_key)) {
    throw std::runtime_error("No boundary condition specified for top");
  }
  auto new_bc = parse_bc(config, bc_key);
  set_top_bc(std::move(new_bc));
}

void FlowField::set_bottom_bc(std::unique_ptr<VelocityBC> bc) const
{
  bottom_bc = std::move(bc);
}

void FlowField::set_top_bc(std::unique_ptr<VelocityBC> bc) const
{
  top_bc = std::move(bc);
}

void FlowField::save(HighFive::File& h5_file) const
{
  Kokkos::fence();

  const auto& partition = mesh.partition();
  auto        nz        = local_end(u.x, 2);

  using time_t = std::decay_t<decltype(time)>;
  auto field =
    h5_file.createDataSet<time_t>("time", HighFive::DataSpace::From(time));
  if (mesh.comm().rank() == 0) {
    field.write(time);
  }

  mesh.save_to_file(h5_file);

  io::hdf5::write3D_xyz(
    h5_file, partition, "u", subview(u.x, ALL, ALL, index_range(0, nz)).view());

  io::hdf5::write3D_xyz(
    h5_file, partition, "v", subview(u.y, ALL, ALL, index_range(0, nz)).view());

  io::hdf5::write3D_xyz(
    h5_file, partition, "w", subview(u.z, ALL, ALL, index_range(0, nz)).view());

  io::hdf5::write3D_xyz(
    h5_file, partition, "pp", subview(pp, ALL, ALL, index_range(0, nz)).view());

  if (nu_t.span() > 1) {
    io::hdf5::write3D_xyz(h5_file,
                          partition,
                          "nu_t",
                          subview(nu_t, ALL, ALL, index_range(0, nz)).view());
  }

  scalars.save(h5_file);

  if (mesh.comm().rank() == 0) {
    logger->info("Solution saved to {}", h5_file.getName());
  }

  Kokkos::fence();
}

void FlowField::save(std::filesystem::path filename) const
{
  if (!filename.has_filename()) {
    throw std::invalid_argument("output filename " + filename.string()
                                + " is not valid");
  }

  HighFive::File h5file = io::hdf5::open_file_with_mpi(
    filename.string(), HighFive::File::Overwrite, mesh.comm().raw_handle());

  save(h5file);
}

void FlowField::initialize_scalars(ConfigTable const& config)
{
  if (!config.contains("scalars")) return;

  auto scalar_configs = config.extract_array_of_tables("scalars");
  scalars.initialize(mesh, scalar_configs);
}

} // namespace solver
} // namespace alps
