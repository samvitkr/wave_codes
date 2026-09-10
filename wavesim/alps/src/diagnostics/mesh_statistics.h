//
// Created by xuananqing on 4/1/23.
//

#pragma once

#include "field_statistics.h"
#include <common/real_type.h>
#include <solvers/mesh/mesh_fwd.h>

#include <Kokkos_Core.hpp>

namespace alps::diagnostics {

/// A class for computing and writing out mesh coordinates conforming a field
/// statistics
class MeshStatistics : public FieldStatistics
{
 public:
  using base_t = FieldStatistics;
  using mesh_t = alps::solver::Mesh;

  explicit MeshStatistics(mesh_t const& mesh_);

  void calculate() override;

  std::vector<int> get_shape(FieldStatisticsType type) const override;

  std::vector<std::string>
  get_variable_labels(FieldStatisticsType type) const override;

  void write_to_tecplot(std::ostream&       out,
                        FieldStatisticsType type) const override;

 protected:
  ///  Write coordinates after averaging the mesh in the x- and y-directions
  void write_to_tecplot_xy(std::ostream& out) const;

  Kokkos::View<Real*, Kokkos::LayoutLeft, Kokkos::HostSpace> gathered_zz;

 private:
  ///  Write coordinates after averaging the mesh in the y-direction
  void write_to_tecplot_y(std::ostream& out) const;

  mesh_t const& mesh;
};

/// A class for computing and writing out the coordinates of a curvilinear mesh
/// needed by a field statistics
template<typename MT>
class CurvilinearMeshStatistics : public MeshStatistics
{
 public:
  using base_t = MeshStatistics;
  using mesh_t = MT;

  explicit CurvilinearMeshStatistics(mesh_t const& mesh_);

  void calculate() override;

  void write_to_tecplot(std::ostream&       out,
                        FieldStatisticsType type) const override;

  void cleanup() override;

  ~CurvilinearMeshStatistics() override;

 protected:
  ///  Write coordinates after averaging the mesh in the y-direction
  void write_to_tecplot_y(std::ostream& out) const;

 private:
  mesh_t const& mesh;

  Kokkos::View<Real**, Kokkos::LayoutLeft, Kokkos::HostSpace> invJ_mean_y;
};
} // namespace alps::diagnostics
