//
// Created by xuanx004 on 10/4/22.
//

#include "solver.h"

#include <common/base/logging.h>
#include <decomp/ghost_cell_exchange.h>
#include <io/hdf5.h>

namespace alps::solver {

void FlowOverWaveSolver::save(std::filesystem::path grid_file,
                              std::filesystem::path data_file,
                              std::filesystem::path aux_file) const
{
  Kokkos::fence();

  const auto& mesh = flow_field.mesh;
  mesh.save_to_file(grid_file);

  flow_field.save(data_file);

  {
    auto const& partition = mesh.partition();
    auto        aux_data  = io::hdf5::open_file_with_mpi(
      aux_file, HighFive::File::Overwrite, mesh.comm().raw_handle());

    if (this->Ru_saved.x.span() > 1) {
      io::hdf5::write3D_xyz(aux_data, partition, "hu", Ru_saved.x.view());
      io::hdf5::write3D_xyz(aux_data, partition, "hv", Ru_saved.y.view());
      io::hdf5::write3D_xyz(aux_data, partition, "hw", Ru_saved.z.view());
    }

    for (std::size_t is = 0; is < flow_field.scalars.size(); ++is) {
      auto const dst_name = fmt::format("/hc{}", is);
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

void FlowOverWaveSolver::load(std::filesystem::path grid_file,
                              std::filesystem::path data_file,
                              std::filesystem::path aux_file)
{
  Kokkos::fence();

  auto const& flow      = flow_field;
  auto const& mesh      = flow.mesh;
  auto const& partition = mesh.partition();
  if (mesh.comm().rank() == 0) {
    logger->info("Loading from " + grid_file.string() + ", "
                 + data_file.string() + ", " + aux_file.string());
  }

  mesh.read_from_file(grid_file);

  const auto nz = mesh.extent(2);

  std::vector<int> scalar_loaded(flow_field.scalars.size(), 0);
  if (std::filesystem::exists(data_file)) {
    auto h5file = io::hdf5::open_file_with_mpi(
      data_file.string(), HighFive::File::ReadOnly, mesh.comm().raw_handle());

    auto field = h5file.getDataSet("time");
    field.read(flow_field.time);

    auto read_3d_data = [&](std::string const&     name,
                            MDView<Real***> const& variable) {
      io::hdf5::read3D_xyz(h5file, partition, name, variable);
    };
    read_3d_data("u", subview(flow.u.x, ALL, ALL, index_range(0, nz)).view());
    read_3d_data("v", subview(flow.u.y, ALL, ALL, index_range(0, nz)).view());
    read_3d_data("w", subview(flow.u.z, ALL, ALL, index_range(0, nz)).view());
    read_3d_data("pp", subview(flow.pp, ALL, ALL, index_range(0, nz)).view());

    if (h5file.exist("eta")) {
      const std::vector total_shape_xy{mesh.global_extent(0),
                                       mesh.global_extent(1)};
      const std::vector block_shape_xy{mesh.extent(0), mesh.extent(1)};
      const std::vector offset_xy{partition.offsets[0], partition.offsets[1]};
      io::hdf5::read_blocks(
        h5file, "eta", flow.eta, total_shape_xy, block_shape_xy, offset_xy);
      io::hdf5::read_blocks(
        h5file, "eta_t", mesh.et, total_shape_xy, block_shape_xy, offset_xy);
    } else {
      if (mesh.comm().rank() == 0) {
        logger->warn("'eta' and 'eta_t' are not found in the data file; they "
                     "are set to 0.");
      }
      Kokkos::deep_copy(flow.eta, 0);
      Kokkos::deep_copy(flow.mesh.et, 0);
    }

    // Load scalar fields
    for (std::size_t is = 0; is < flow_field.scalars.size(); ++is) {
      auto dst_name = fmt::format("c{}", is);
      if (!h5file.exist(dst_name)) continue;

      if (mesh.comm().rank() == 0) {
        logger->info("Loading scalar {}", dst_name);
      }
      flow_field.scalars.at(is).load(h5file, dst_name, logger);
      scalar_loaded.at(is) = 1;
    }
  } else {
    throw std::runtime_error("data file '" + data_file.string()
                             + "' does not exist.");
  }

  update_halo_z(partition, flow.u.x, 1);
  update_halo_z(partition, flow.u.y, 1);
  update_halo_z(partition, flow.u.z, 1);
  update_halo_z(partition, flow.pp, 1);
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
      auto const dst_name = fmt::format("hc{}", is);
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
        Kokkos::view_alloc(fmt::format("Rc{}", is),
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
