#pragma once

#include <io/io_fwd.h>

// Some versions of HDF5 headers defined OMPI_SKIP_MPICXX and the compiler may
// complain about redefinition, the following lines wrap HDF5 includes inside a
// guard
#if defined(OMPI_SKIP_MPICXX)
#define OLD_OMPI_SKIP_MPICXX OMPI_SKIP_MPICXX
#undef OMPI_SKIP_MPICXX
#endif
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>
#if defined(OLD_OMPI_SKIP_MPICXX)
#if !defined(OMPI_SKIP_MPICXX)
#define OMPI_SKIP_MPICXX OLD_OMPI_SKIP_MPICXX
#endif
#undef OLD_OMPI_SKIP_MPICXX
#endif

#include <common/base/macros.h>
#include <decomp/block_partition.h>

#include <Kokkos_Core.hpp>

#include <vector>

namespace alps {
namespace io::hdf5 {

namespace detail {
std::vector<size_t> left_to_right_layout(std::vector<int> const& arr);

void check_global_shape(HighFive::DataSet const&   dset,
                        std::vector<size_t> const& expected);
} // namespace detail

HighFive::File open_file_with_mpi(const std::string&         filename,
                                  HighFive::File::AccessMode openFlags,
                                  MPI_Comm                   comm,
                                  MPI_Info info = MPI_INFO_NULL);

/// Write a distributed Kokkos View in parallel to HDF5 file
/**
 *  @note If the View is inaccessible from the host, it will be copied to the
 * host before writing
 *
 *  @param file The HDF5 file to write to
 *  @param dataset_name The name of the dataset to write to
 *  @param data The data to write (must be a Kokkos View with LayoutLeft)
 *  @param total_shape The shape of the data (in LayoutLeft order)
 *  @param block_shape The shape of the local blocks (in LayoutLeft order)
 *  @param block_offset The offset of the blocks (in LayoutLeft order)
 */
template<typename DT, typename... ViewProperties>
auto write_blocks(HighFive::File& file,
                  std::string     dataset_name,
                  Kokkos::View<DT, Kokkos::LayoutLeft, ViewProperties...> data,
                  std::vector<int> total_shape,
                  std::vector<int> block_shape,
                  std::vector<int> block_offset)
{
  auto constexpr n_dims = decltype(data)::rank;
  if (total_shape.size() != n_dims || block_shape.size() != n_dims
      || block_offset.size() != n_dims) {
    throw std::runtime_error("total_shape, block_shape and block_offset must "
                             "have the same number of dimensions as the data");
  }
  // Change the shape and offsets to C (LayoutRight) order
  auto total  = detail::left_to_right_layout(total_shape);
  auto count  = detail::left_to_right_layout(block_shape);
  auto offset = detail::left_to_right_layout(block_offset);

  using T   = typename decltype(data)::non_const_value_type;
  auto dset = file.createDataSet<T>(dataset_name, HighFive::DataSpace(total));

  auto xfer_prop = HighFive::DataTransferProps{};
  xfer_prop.add(HighFive::UseCollectiveIO{});
  auto data_host = Kokkos::create_mirror_view(data);
  Kokkos::deep_copy(data_host, data);
  dset.select(offset, count).write_raw(data_host.data(), xfer_prop);
  return dset;
}

/// Read blocks from a HDF5 dataset in parallel into Kokkos View
/**
 *  @note If the View is inaccessible from the host, the data will be read into
 * host memory and then copied to the device
 *
 *  @param file The HDF5 file to read from
 *  @param dataset_name The name of the dataset to read from
 *  @param data The view to be read into (must be a Kokkos View with LayoutLeft)
 *  @param total_shape The shape of the data (in LayoutLeft order)
 *  @param block_shape The shape of the local blocks (in LayoutLeft order)
 *  @param block_offset The offset of the blocks (in LayoutLeft order)
 */
template<typename DT, typename... ViewProperties>
auto read_blocks(HighFive::File const& file,
                 std::string           dataset_name,
                 Kokkos::View<DT, Kokkos::LayoutLeft, ViewProperties...> data,
                 std::vector<int> total_shape,
                 std::vector<int> block_shape,
                 std::vector<int> block_offset)
{
  auto constexpr n_dims = decltype(data)::rank;
  if (total_shape.size() != n_dims || block_shape.size() != n_dims
      || block_offset.size() != n_dims) {
    throw std::runtime_error("total_shape, block_shape and block_offset must "
                             "have the same number of dimensions as the data");
  }
  // Change the shape and offsets to C (LayoutRight) order
  auto total  = detail::left_to_right_layout(total_shape);
  auto count  = detail::left_to_right_layout(block_shape);
  auto offset = detail::left_to_right_layout(block_offset);

  const auto dset = file.getDataSet(dataset_name);
  detail::check_global_shape(dset, total);

  auto xfer_prop = HighFive::DataTransferProps{};
  xfer_prop.add(HighFive::UseCollectiveIO{});
  auto data_host = Kokkos::create_mirror_view(data);
  dset.select(offset, count).read_raw(data_host.data(), xfer_prop);
  Kokkos::deep_copy(data, data_host);
  return dset;
}

template<typename DT, typename... ViewProperties>
auto write3D_xyz(HighFive::File&       file,
                 const BlockPartition& decomp,
                 std::string           dataset_name,
                 Kokkos::View<DT, Kokkos::LayoutLeft, ViewProperties...> data)
{
  auto constexpr n_dims = decltype(data)::rank;
  static_assert(n_dims == 3, "write3D_xyz only supports 3D Views");

  const std::vector<int> shape{
    decomp.extents[0], decomp.extents[1], decomp.extents[2]};
  const std::vector<int> offset{
    decomp.offsets[0], decomp.offsets[1], decomp.offsets[2]};
  const std::vector<int> global{decomp.global_extents[0],
                                decomp.global_extents[1],
                                decomp.global_extents[2]};
  return write_blocks(file, dataset_name, data, global, shape, offset);
}

template<typename DT, typename... ViewProperties>
auto write1D_z(HighFive::File&       file,
               const BlockPartition& decomp,
               std::string           dataset_name,
               Kokkos::View<DT, Kokkos::LayoutLeft, ViewProperties...> data,
               bool                                                    output)
{
  auto constexpr n_dims = decltype(data)::rank;
  static_assert(n_dims == 1, "write1D_z only supports 1D Views");

  const std::vector<int> shape{output ? decomp.extents[2] : 0};
  const std::vector<int> offset{output ? decomp.offsets[2] : 0};
  const std::vector<int> global{decomp.global_extents[2]};

  return write_blocks(file, dataset_name, data, global, shape, offset);
}

template<typename DT, typename... ViewProperties>
auto read3D_xyz(HighFive::File const& file,
                const BlockPartition& decomp,
                std::string           dataset_name,
                Kokkos::View<DT, Kokkos::LayoutLeft, ViewProperties...> data)
{
  auto constexpr n_dims = decltype(data)::rank;
  static_assert(n_dims == 3, "read3D_xyz only supports 3D Views");

  const std::vector<int> shape{
    decomp.extents[0], decomp.extents[1], decomp.extents[2]};
  const std::vector<int> offset{
    decomp.offsets[0], decomp.offsets[1], decomp.offsets[2]};
  const std::vector<int> global{decomp.global_extents[0],
                                decomp.global_extents[1],
                                decomp.global_extents[2]};
  return read_blocks(file, dataset_name, data, global, shape, offset);
}

template<typename DT, typename... ViewProperties>
auto read1D_z(HighFive::File const& file,
              const BlockPartition& decomp,
              std::string           dataset_name,
              Kokkos::View<DT, Kokkos::LayoutLeft, ViewProperties...> data)
{
  auto constexpr n_dims = decltype(data)::rank;
  static_assert(n_dims == 1, "read1D_z only supports 1D Views");

  const std::vector<int> shape{decomp.extents[2]};
  const std::vector<int> offset{decomp.offsets[2]};
  const std::vector<int> global{decomp.global_extents[2]};

  return read_blocks(file, dataset_name, data, global, shape, offset);
}

HighFive::DataSet write_string(HighFive::File& file,
                               std::string     dataset_name,
                               std::string     content);

} // namespace io::hdf5
} // namespace alps
