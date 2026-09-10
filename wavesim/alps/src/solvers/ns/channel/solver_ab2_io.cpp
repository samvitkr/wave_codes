#include "solver_ab2.h"

#include <common/base/logging.h>
#include <common/container/view_utils.h>
#include <decomp/ghost_cell_exchange.h>
#include <io/hdf5.h>
#include <solvers/field/flow_field.h>

#include <Kokkos_Core.hpp>

#include <stdexcept>

namespace alps::solver {

void ChannelFlowSolverAB2::save(std::filesystem::path const grid_file,
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

    for (std::size_t is = 0; is < flow_field.scalars.size(); ++is) {
      auto const dst_name = "/hc" + std::to_string(is);
      if (Rc_saved.at(is).span() > 1) {
        io::hdf5::write3D_xyz(aux_data, partition, dst_name, Rc_saved.at(is));
      }
    }
  }

  if (mesh.comm().rank() == 0) {
    logger->info("Restart saved to " + grid_file.string() + ", "
                 + data_file.string() + ", " + aux_file.string());
  }

  Kokkos::fence();
}

void ChannelFlowSolverAB2::load(std::filesystem::path const grid_file,
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

  std::vector<int> scalar_loaded(flow_field.scalars.size(), 0);
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

    // Load scalar fields
    for (std::size_t is = 0; is < flow_field.scalars.size(); ++is) {
      auto dst_name = "c" + std::to_string(is);
      if (!h5file.exist(dst_name)) continue;

      if (mesh.comm().rank() == 0) {
        logger->info("Loading scalar {}", dst_name);
      }
      flow_field.scalars.at(is).load(h5file, dst_name, logger);
      scalar_loaded.at(is) = 1;
    }
  } else {
    throw std::runtime_error("data file '" + data_file.string()
                             + "' does not exist");
  }

  update_halo_z(partition, flow_field.u.x, 1);
  update_halo_z(partition, flow_field.u.y, 1);
  update_halo_z(partition, flow_field.u.z, 1);
  update_halo_z(partition, flow_field.pp, 1);
  for (auto const& scalar : flow_field.scalars.storage()) {
    update_halo_z(partition, scalar.array, 1);
  }

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

    for (std::size_t is = 0; is < flow_field.scalars.size(); ++is) {
      auto const dst_name = "hc" + std::to_string(is);
      if (!aux_data.exist(dst_name)) continue;

      if (scalar_loaded.at(is) == 0) {
        if (mesh.comm().rank() == 0) {
          logger->error("Ignoring {} for scalar c{} because c{} is not "
                        "present in restart data",
                        dst_name,
                        is,
                        is);
        }
        continue;
      }

      if (mesh.comm().rank() == 0) {
        logger->info("Loading scalar rhs {}", dst_name);
      }
      Rc_saved.at(is) = MDView<Real***, default_memory_pool>(
        Kokkos::view_alloc("Rc" + std::to_string(is),
                           Kokkos::WithoutInitializing),
        mesh.extent(0),
        mesh.extent(1),
        mesh.extent(2));
      io::hdf5::read3D_xyz(aux_data, partition, dst_name, Rc_saved.at(is));
    }
  }

  Kokkos::fence();
}

} // namespace alps::solver
