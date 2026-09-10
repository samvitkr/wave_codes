#pragma once

#include "flow_field.h"
#include "scalar_fields.h"
#include <common/base/logging_fwd.h>
#include <common/container/vector_field.h>
#include <common/container/view_types.h>
#include <io/io_fwd.h>
#include <solvers/field/bc_types.h>
#include <solvers/field/traits.h>
#include <solvers/mesh/curvilinear_mesh.h>

#include <filesystem>
#include <memory>

namespace alps {
namespace solver {

class FlowOverWaveField
{
 public:
  explicit FlowOverWaveField(const BottomWaveMesh& local_mesh);

  // cast to FlowField
  [[nodiscard]] FlowField as_channel_flow_field() const;

  /// Set the boundary conditions by parsing the config file
  void parse_bcs_from(const ConfigTable& config,
                      WhichBoundary      which = WhichBoundary::Both);

  /// Set the top boundary condition by parsing the config file
  void parse_top_bc_from(const ConfigTable& config);

  /// Set the bottom boundary condition by parsing the config file
  void parse_bottom_bc_from(const ConfigTable& config);

  /// Parse a boundary condition from an item in the config file
  std::unique_ptr<VelocityBC> parse_bc(const ConfigTable& config,
                                       std::string        bc_key) const;

  /// Set the bottom boundary to a given BC class
  /** @note This function is const to allow the BC to be changed in the solver
   */
  void set_bottom_bc(std::unique_ptr<VelocityBC> bc) const;

  /// Set the top boundary to a given BC class
  /** @note This function is const to allow the BC to be changed in the solver
   */
  void set_top_bc(std::unique_ptr<VelocityBC> bc) const;

  void initialize_scalars(ConfigTable const& config);

  void save(HighFive::File& h5_file) const;

  void save(std::filesystem::path filename) const;

  const BottomWaveMesh& mesh;
  Vector3Field<Real***> u;
  HaloView<Real***>     pp;
  MDView<Real**>        eta;

  mutable HaloView<Real***> nu_t;

  mutable std::unique_ptr<VelocityBC> top_bc;
  mutable std::unique_ptr<VelocityBC> bottom_bc;

  ScalarFields scalars;

  mutable double time{};

 private:
  Logger logger;
};

template<>
struct is_curvilinear_field<FlowOverWaveField> : std::true_type
{};

} // namespace solver
} // namespace alps
