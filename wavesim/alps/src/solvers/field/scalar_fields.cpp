//
// Created by xuanx004 on 7/12/24.
//

#include "scalar_fields.h"

#include <common/base/logging.h>
#include <common/container/view_utils.h>
#include <common/kokkos_abstraction/pool_space.h>
#include <common/utils/to_lower_case.h>
#include <io/hdf5.h>
#include <solvers/mesh/mesh.h>

namespace alps::solver {

void ScalarField::load(HighFive::File const& file,
                       std::string const&    path,
                       Logger const&         logger) const
{
  // do nothing if the dataset does not exist
  if (!file.exist(path)) return;

  auto const ds = io::hdf5::read3D_xyz(
    file,
    mesh->partition(),
    path,
    subview(array, Kokkos::ALL, Kokkos::ALL, std::pair(0, mesh->extent(2)))
      .view());

  auto const attr_ds_name = path + "_label";
  if (!file.exist(attr_ds_name)) return;
  auto const attr = file.getDataSet(attr_ds_name);
  if (mesh->comm().comm.rank() == 0) {
    // the fixed-length string stored in HDF5 is null padded
    // need to re-construct a string to remove the padding
    auto const label_padded = attr.read<std::string>();
    auto const label        = std::string(label_padded.c_str());
    if (!label.empty() && label != array.label()) {
      logger->warn(
        "Label mismatch: expected '{}', got '{}'", array.label(), label);
    }
  }
}

void ScalarField::save(HighFive::File& file, std::string const& path) const
{
  auto const& partition = mesh->partition();

  auto ds = io::hdf5::write3D_xyz(
    file, partition, path, create_inner_view(array).view());

  // the scalar label is stored as a fixed-length string in HDF5
  auto const attr_ds_name = path + "_label";
  auto const scalar_dataspace =
    HighFive::DataSpace(HighFive::DataSpace::dataspace_scalar);
  auto const fixed_str_nullpad =
    HighFive::FixedLengthStringType(128, HighFive::StringPadding::NullPadded);
  auto attr =
    file.createDataSet(attr_ds_name, scalar_dataspace, fixed_str_nullpad);
  if (partition.comm.rank() == 0) {
    attr.write(label());
  }

  if (nuD.span() > 1) {
    io::hdf5::write3D_xyz(
      file, partition, path + "_nuD", create_inner_view(nuD).view());
  }
}

void ScalarFields::push_back(ScalarField field)
{
  fields_.push_back(std::move(field));
}

void ScalarFields::clear()
{
  fields_.clear();
}

ScalarField* ScalarFields::find(std::string const& name) noexcept
{
  for (auto& field : fields_) {
    if (field.label() == name) {
      return &field;
    }
  }
  return nullptr;
}

ScalarField const* ScalarFields::find(std::string const& name) const noexcept
{
  for (auto const& field : fields_) {
    if (field.label() == name) {
      return &field;
    }
  }
  return nullptr;
}

std::size_t ScalarFields::size() const noexcept
{
  return fields_.size();
}

std::unique_ptr<ScalarBC> parse_scalar_bc(const ConfigTable& config,
                                          std::string        bc_key)
{
  auto const bc_config = config.extract_table(bc_key);
  auto const bc_type   = bc_config.get_value<std::string>("type");
  // convert std::string bc_type to lower case
  auto const bc_type_lower = to_lower_case(bc_type);

  if (bc_type_lower == "constantdirichlet") {
    auto bc = std::make_unique<ConstantDirichletBC>();
    bc->load_from(bc_config);
    return bc;
  }
  if (bc_type_lower == "constantgradient") {
    auto bc = std::make_unique<ConstantGradientBC>();
    bc->load_from(bc_config);
    return bc;
  }
  if (bc_type_lower == "constantflux") {
    auto bc = std::make_unique<ConstantFluxBC>();
    bc->load_from(bc_config);
    return bc;
  }
  if (bc_type_lower == "custom" or bc_type_lower == "coded") {
    return nullptr;
  }

  throw std::runtime_error("Cannot handle boundary condition of type " + bc_type
                           + " at key " + bc_key);
}

void ScalarFields::initialize(Mesh const&                     mesh,
                              std::vector<ConfigTable> const& scalar_configs)
{
  if (!fields_.empty()) {
    throw std::runtime_error("ScalarFields already initialized");
  }

  auto const n_scalars = scalar_configs.size();
  // make sure labels are unique
  std::vector<std::string> labels;
  for (std::size_t scalar_i = 0; scalar_i < n_scalars; ++scalar_i) {
    auto label = scalar_configs[scalar_i].get_value_or<std::string>(
      "label", "c" + std::to_string(scalar_i));
    for (auto const& existing_label : labels) {
      if (label == existing_label) {
        throw std::runtime_error("Duplicate scalar label: " + label);
      }
    }
    labels.emplace_back(label);
  }

  for (std::size_t scalar_i = 0; scalar_i < n_scalars; ++scalar_i) {
    ScalarField field;
    field.unique_id = scalar_i;
    field.array =
      HaloView<Real***, default_memory_pool>(labels.at(scalar_i),
                                             std::pair(0, mesh.extent(0) - 1),
                                             std::pair(0, mesh.extent(1) - 1),
                                             std::pair(-1, mesh.extent(2)));
    field.mesh = &mesh;
    fields_.push_back(std::move(field));
  }

  // process boundary conditions
  for (std::size_t scalar_i = 0; scalar_i < n_scalars; ++scalar_i) {
    auto const& config = scalar_configs[scalar_i];
    if (auto constexpr bc_key = "BC.bottom"; config.contains(bc_key)) {
      fields_[scalar_i].bottom_bc = parse_scalar_bc(config, bc_key);
    }
    if (auto constexpr bc_key = "BC.top"; config.contains(bc_key)) {
      fields_[scalar_i].top_bc = parse_scalar_bc(config, bc_key);
    }
  }
}

void ScalarFields::save(HighFive::File& file) const
{
  for (std::size_t is = 0; is < size(); ++is) {
    auto const dst_name = "c" + std::to_string(is);
    at(is).save(file, dst_name);
  }
}

ScalarFields ScalarFields::shallow_copy() const
{
  ScalarFields copy;
  for (auto const& field : fields_) {
    ScalarField copy_field;
    copy_field.unique_id = field.unique_id;
    copy_field.array     = field.array;
    copy_field.nuD       = field.nuD;
    copy_field.mesh      = field.mesh;
    copy.push_back(std::move(copy_field));
  }
  return copy;
}

} // namespace alps::solver
