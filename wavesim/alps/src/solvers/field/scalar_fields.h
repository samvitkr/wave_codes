//
// Created by xuanx004 on 7/12/24.
//

#pragma once

#include <common/base/logging_fwd.h>
#include <common/container/view_types.h>
#include <common/program_options/config_table.h>
#include <common/real_type.h>
#include <solvers/field/scalar_bc_types.h>

#include <vector>

// Forward declarations
namespace HighFive {
class File;
} // namespace HighFive

namespace alps::solver {

class Mesh;
class FlowOverWaveField;
class FreeSurfaceFlowField;

class ScalarField
{
 public:
  std::size_t               unique_id;
  HaloView<Real***>         array;
  std::unique_ptr<ScalarBC> top_bc;
  std::unique_ptr<ScalarBC> bottom_bc;

  mutable HaloView<Real***> nuD; // mutable to allow lazy allocation

  Mesh const* mesh{nullptr};

  std::string label() const { return array.label(); }

  void load(HighFive::File const& file,
            std::string const&    path,
            Logger const&         logger) const;

  void save(HighFive::File& file, std::string const& path) const;
};

class ScalarFields
{
 public:
  ScalarField const& at(std::size_t i) const { return fields_.at(i); }

  ScalarField& at(std::size_t i) { return fields_.at(i); }

  void push_back(ScalarField field);

  void clear();

  ScalarField* find(std::string const& name) noexcept;

  ScalarField const* find(std::string const& name) const noexcept;

  std::size_t size() const noexcept;

  auto& storage() { return fields_; }

  auto const& storage() const { return fields_; }

  void initialize(Mesh const&                     mesh,
                  std::vector<ConfigTable> const& scalar_configs);

  void save(HighFive::File& file) const;

  auto begin() { return fields_.begin(); }

  auto end() { return fields_.end(); }

  auto begin() const { return fields_.begin(); }

  auto end() const { return fields_.end(); }

 private:
  ScalarFields shallow_copy() const;

  std::vector<ScalarField> fields_;

  friend class FlowOverWaveField;
  friend class FreeSurfaceFlowField;
};

} // namespace alps::solver
