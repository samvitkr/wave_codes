#include "hdf5.h"

#include <common/base/logging.h>

namespace alps::io::hdf5 {
HighFive::File open_file_with_mpi(const std::string&         filename,
                                  HighFive::File::AccessMode openFlags,
                                  MPI_Comm                   comm,
                                  MPI_Info                   info)
{
  auto fapl = HighFive::FileAccessProps();
  fapl.add(HighFive::MPIOFileAccess(comm, info));
  return HighFive::File(filename, openFlags, fapl);
}

HighFive::DataSet write_string(HighFive::File& file,
                               std::string     dataset_name,
                               std::string     content)
{
  auto dataset = file.createDataSet<std::string>(
    dataset_name, HighFive::DataSpace::From(content));
  dataset.write(content);
  return dataset;
}

namespace detail {
std::vector<size_t> left_to_right_layout(std::vector<int> const& arr)
{
  std::vector<size_t> result(arr.size());
  for (size_t i = 0; i < arr.size(); ++i) {
    if (arr[i] < 0) {
      throw std::runtime_error(
        "Negative dimension sizes or offsets are not allowed");
    }
    result[arr.size() - 1 - i] = static_cast<size_t>(arr[i]);
  }
  return result;
}

void check_global_shape(HighFive::DataSet const&   dset,
                        std::vector<size_t> const& expected)
{
  const auto dset_dims = dset.getDimensions();
  const auto dset_path = dset.getPath();
  if (dset_dims.size() != expected.size()) {
    throw std::runtime_error(
      "Dataset " + dset_path
      + " has different number of dimensions than expected");
  }
  int all_equal = 1;
  for (size_t i = 0; i < expected.size(); ++i) {
    if ((std::int64_t)dset_dims[i] < (std::int64_t)expected[i]) {
      all_equal = -1;
      break;
    }
    if ((std::int64_t)dset_dims[i] != (std::int64_t)expected[i]) {
      all_equal = 0;
    }
  }
  if (all_equal < 0) {
    throw std::runtime_error("HDF5 dataset " + dset_path
                             + " has smaller shape than expected");
  }
  if (all_equal == 0) {
    auto logger = get_logger("hdf5");
    logger->warn("One or more dimensions of the HDF5 dataset " + dset_path
                 + " is larger than expected");
  }
}
} // namespace detail

} // namespace alps::io::hdf5
