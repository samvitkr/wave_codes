#include "mesh.h"

#include <common/container/view_utils.h>
#include <common/kokkos_abstraction/pool_space.h>
#include <io/hdf5.h>
#include <spectral/spectral.h>

#include <Kokkos_Core.hpp>
#include <fmt/format.h>
#include <mpipp/collectives.h>

namespace alps {
namespace solver {

Mesh::Mesh(const Grid& grid_, Real Lz, int n_ghost)
  : grid{grid_}
  , pex{static_cast<Real>(grid_.pex)}
  , pey{static_cast<Real>(grid_.pey)}
  , hbar{Lz}
  , zz(HaloView<Real*, default_memory_pool>(
      "zeta",
      std::pair(-n_ghost, grid.extent(2) + n_ghost - 1)))
  , dz(HaloView<Real*, default_memory_pool>(
      "dzeta",
      std::pair(-n_ghost, grid.extent(2) + n_ghost - 1)))
  , zw(HaloView<Real*, default_memory_pool>(
      "zeta@w",
      std::pair{-n_ghost, grid.extent(2) + n_ghost - 1}))
  , dzw(HaloView<Real*, default_memory_pool>(
      "dzeta@w",
      std::pair{-n_ghost, grid.extent(2) + n_ghost - 1}))
  , zz_h(HaloView<Real*, default_host_memory_pool>(
      "zeta",
      std::pair{-n_ghost, grid.extent(2) + n_ghost - 1}))
  , dz_h(HaloView<Real*, default_host_memory_pool>(
      "dzeta",
      std::pair{-n_ghost, grid.extent(2) + n_ghost - 1}))
  , zw_h(HaloView<Real*, default_host_memory_pool>(
      "zeta@w",
      std::pair{-n_ghost, grid.extent(2) + n_ghost - 1}))
  , dzw_h(HaloView<Real*, default_host_memory_pool>(
      "dzeta@w",
      std::pair{-n_ghost, grid.extent(2) + n_ghost - 1}))
{}

namespace {
void convert_node_mesh_to_centers(nonstd::span<Real>       z_center,
                                  nonstd::span<Real const> z_node)
{
  auto const nz = z_center.size();
  if (nz < 3) {
    throw std::invalid_argument("Size of the input mesh too small");
  }

  for (std::size_t i = 1; i < nz - 1; ++i) {
    z_center[i] = (z_node[i] + z_node[i - 1]) / 2;
  }
  z_center[0]      = 0;
  z_center[nz - 1] = 1;
}

void calculate_mesh_spacings(nonstd::span<Real> dz, nonstd::span<Real const> z)
{
  auto const nz = z.size();
  if (nz < 2) {
    throw std::invalid_argument("Size of the input mesh too small");
  }

  for (std::size_t i = 0; i < nz - 1; ++i) {
    dz[i] = z[i + 1] - z[i];
  }
  dz[nz - 1] = 0;
}
} // namespace

void Mesh::load(nonstd::span<Real const> z_node) const
{
  // verify the input grid
  if ((int)z_node.size() != grid.global_extent(2)) {
    throw std::invalid_argument(
      fmt::format("Size of the input mesh (Nz={}) is incompatible with the "
                  "assigned mesh (Nz={})",
                  z_node.size(),
                  grid.global_extent(2)));
  }
  auto       z_node_h = std::vector<Real>(z_node.begin(), z_node.end());
  auto const nz       = z_node.size();
  if (z_node_h.at(0) != 0
      || std::abs(z_node_h[nz - 2] - 1)
           > std::numeric_limits<Real>::epsilon()) {
    throw std::invalid_argument(
      fmt::format("Mesh must be normalized to 0..1 (Input: {}..{})",
                  z_node_h[0],
                  z_node_h[nz - 2]));
  }
  z_node_h[nz - 1] = z_node_h[nz - 2];

  std::vector<Real> z_center_h(z_node_h.size(), 0);
  std::vector<Real> dz_center(z_node_h.size(), 0);
  std::vector<Real> dz_node(z_node_h.size(), 0);

  convert_node_mesh_to_centers(z_center_h, z_node_h);
  calculate_mesh_spacings(dz_center, z_center_h);
  calculate_mesh_spacings(dz_node, z_node_h);

  auto k1 = zz.begin(0);
  auto k2 = zz.end(0);
  if (grid.comm().is_first(2)) k1 = 0;
  if (grid.comm().is_last(2)) k2 = local_end(zz, 0);

  for (auto k = k1; k < k2; ++k) {
    auto const offset = grid.offset(2);
    zz_h(k)           = z_center_h[k + offset];
    zw_h(k)           = z_node_h[k + offset];
    dz_h(k)           = dz_center[k + offset];
    dzw_h(k)          = dz_node[k + offset];
  }

  // set boundary values
  // spacings beyond the boundaries are set to zero
  if (grid.comm().is_first(2)) { // bottom boundary
    for (auto k = zz_h.begin(0); k < 0; ++k) {
      zz_h(k)  = 0;
      zw_h(k)  = 0;
      dz_h(k)  = 0;
      dzw_h(k) = 0;
    }
  }
  if (grid.comm().is_last(2)) { // top boundary
    for (auto k = local_end(zz_h, 0) + 1; k < end(zz_h, 0); ++k) {
      zz_h(k) = zz_h(local_end(zz_h, 0));
    }
    for (auto k = local_end(zw_h, 0); k < end(zw_h, 0); ++k) {
      zw_h(k)  = zw_h(local_end(zw_h, 0) - 1);
      dz_h(k)  = 0;
      dzw_h(k) = 0;
    }
    dzw_h(local_end(dzw_h, 0) - 1) = 0;
  }

  deep_copy(zz.view(), zz_h.view());
  deep_copy(zw.view(), zw_h.view());
  deep_copy(dz.view(), dz_h.view());
  deep_copy(dzw.view(), dzw_h.view());
}

void Mesh::save_to_file(HighFive::File& hdf_file) const
{
  auto nz = local_end(zz, 0);

  auto pex_dset = hdf_file.createDataSet<decltype(pex)>(
    "/pex", HighFive::DataSpace::From(pex));
  auto pey_dset = hdf_file.createDataSet<decltype(pey)>(
    "/pey", HighFive::DataSpace::From(pey));
  auto hbar_dset = hdf_file.createDataSet<decltype(hbar)>(
    "/hbar", HighFive::DataSpace::From(hbar));
  if (grid.comm().rank() == 0) {
    pex_dset.write(pex);
    pey_dset.write(pey);
    hbar_dset.write(hbar);
  } // Currently pex, pey and hbar are only saved for redundancy.

  auto const nz_distribution = partition().get_distribution(2);
  auto const nz_offset       = partition().get_offsets(2);

  auto zz_dset =
    hdf_file.createDataSet<Real>("zz", HighFive::DataSpace(global_extent(2)));
  auto zw_dset =
    hdf_file.createDataSet<Real>("zw", HighFive::DataSpace(global_extent(2)));
  auto dz_dset =
    hdf_file.createDataSet<Real>("dz", HighFive::DataSpace(global_extent(2)));
  auto dzw_dset =
    hdf_file.createDataSet<Real>("dzw", HighFive::DataSpace(global_extent(2)));
  for (auto [dset, data] : {std::tie(zz_dset, zz_h),
                            std::tie(zw_dset, zw_h),
                            std::tie(dz_dset, dz_h),
                            std::tie(dzw_dset, dzw_h)}) {
    if (grid.comm().is_first(1)) {
      std::vector<Real> gathered(global_extent(2));
      mpipp::gatherv(nonstd::span(&data(0), nz),
                     gathered.data(),
                     nonstd::span(nz_distribution),
                     nonstd::span(nz_offset),
                     0,
                     comm().axis_comm[2]);
      if (grid.comm().is_first(2)) {
        dset.write(gathered);
      }
    }
  }
}

void Mesh::read_and_interp(const HighFive::File& hdf_file) const
{
  auto const                             nz_global = grid.global_extent(2);
  MDView<Real*, Kokkos::HostSpace> const zw_all("zw all", nz_global);

  const std::vector total_shape_z{nz_global};
  const std::vector offset_xy{0};
  io::hdf5::read_blocks(
    hdf_file, "zw", zw_all, total_shape_z, total_shape_z, offset_xy);
  load(zw_all);
}

void Mesh::save_to_file(std::filesystem::path filename) const
{
  auto h5file = io::hdf5::open_file_with_mpi(
    filename.string(), HighFive::File::Overwrite, grid.comm().raw_handle());

  save_to_file(h5file);
}

void Mesh::read_from_file(std::filesystem::path filename) const
{

  auto h5file = io::hdf5::open_file_with_mpi(
    filename.string(), HighFive::File::ReadOnly, grid.comm().raw_handle());

  read_and_interp(h5file);
}

int Mesh::extent(int axis) const
{
  return grid.extent(axis);
}
Kokkos::Array<int, 3> Mesh::extents() const
{
  return grid.extents();
}
int Mesh::local_extent(int axis) const
{
  return grid.extent(axis);
}
Kokkos::Array<int, 3> Mesh::local_extents() const
{
  return grid.extents();
}
int Mesh::global_extent(int axis) const
{
  return grid.global_extent(axis);
}
Kokkos::Array<int, 3> Mesh::global_extents() const
{
  return {grid.global_extent(0), grid.global_extent(1), grid.global_extent(2)};
}

const MPIComm3D& Mesh::comm() const
{
  return grid.comm();
}

const BlockPartition& Mesh::partition() const
{
  return grid.partition();
}
} // namespace solver
} // namespace alps
