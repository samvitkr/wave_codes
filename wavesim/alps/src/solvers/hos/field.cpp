#include "field.h"

#include <common/kokkos_abstraction/pool_space.h>
#include <spectral/spectral.h>

namespace alps::solver::hos {

HOSField::HOSField(Grid const& spectral_grid)
  : value(MDView<Real** [2], default_memory_pool>(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "hos solution"),
      spectral_grid.extent(0),
      spectral_grid.extent(1)))
  , pa(MDView<Real**, default_memory_pool>(
      Kokkos::view_alloc(Kokkos::WithoutInitializing, "P_air"),
      spectral_grid.extent(0),
      spectral_grid.extent(1)))
  , grid{spectral_grid}
  , time{0}
  , kx0{static_cast<Real>(spectral_grid.pex)}
  , ky0{static_cast<Real>(spectral_grid.pey)}
{
  Kokkos::deep_copy(pa, 0);
}

void HOSField::save(HighFive::File& h5_file) const
{
  Kokkos::fence();

  auto field = h5_file.createDataSet<decltype(this->time)>(
    "time", HighFive::DataSpace::From(this->time));
  auto pex_dset = h5_file.createDataSet<std::decay_t<decltype(kx0)>>(
    "pex", HighFive::DataSpace::From(kx0));
  auto pey_dset = h5_file.createDataSet<std::decay_t<decltype(ky0)>>(
    "pey", HighFive::DataSpace::From(ky0));
  if (grid.comm().rank() == 0) {
    field.write(this->time);
    pex_dset.write(kx0);
    pey_dset.write(ky0);
  }

  const BlockPartition& partition = grid.partition();

  const std::vector total_shape_xy{partition.global_extents[0],
                                   partition.global_extents[1]};
  const std::vector block_shape_xy{partition.extents[0], partition.extents[1]};
  const std::vector offset_xy{partition.offsets[0], partition.offsets[1]};
  io::hdf5::write_blocks(
    h5_file, "eta", eta(), total_shape_xy, block_shape_xy, offset_xy);

  io::hdf5::write_blocks(
    h5_file, "vps", vps(), total_shape_xy, block_shape_xy, offset_xy);

  Kokkos::fence();
}

void HOSField::save(std::filesystem::path filename) const
{
  if (!filename.has_filename()) {
    throw std::invalid_argument("output filename " + filename.string()
                                + " is not valid");
  }

  HighFive::File h5file = io::hdf5::open_file_with_mpi(
    filename.string(), HighFive::File::Overwrite, grid.comm().raw_handle());

  save(h5file);
}

} // namespace alps::solver::hos
