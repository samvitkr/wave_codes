//
// Created by xuananqing on 4/1/23.
//

#include "mesh_statistics.h"

#include <common/kokkos_abstraction/pool_space.h>
#include <decomp/block_partition.h>
#include <linear_algebra/mean_variance.h>
#include <solvers/mesh/curvilinear_mesh.h>
#include <solvers/mesh/mesh.h>

#include <mpipp/collectives.h>
#include <spdlog/fmt/ostr.h>

namespace alps::diagnostics {

MeshStatistics::MeshStatistics(mesh_t const& mesh_)
  : base_t({})
  , mesh{mesh_}
{}

void MeshStatistics::calculate()
{
  if (gathered_zz.extent_int(0) == mesh.global_extent(2)) return;

  // Collect and store the coordinates in the z direction
  auto const& comm            = mesh.comm();
  auto const  nz              = mesh.global_extent(2);
  auto const  nz_distribution = mesh.partition().get_distribution(2);

  gathered_zz = decltype(gathered_zz)("gathered zz", nz);

  mpipp::allgatherv(nonstd::span(&mesh.zz_h(0), mesh.local_extent(2)),
                    &gathered_zz(0),
                    nonstd::span(nz_distribution),
                    comm.axis_comm[2]);
}

std::vector<int> MeshStatistics::get_shape(FieldStatisticsType type) const
{
  // on processors without the gathered statistics, return empty
  auto const& comm = mesh.comm();
  if (!comm.is_first(1) || !comm.is_first(2)) return {};

  if (type == FieldStatisticsType::XY) {
    return {mesh.global_extent(2)};
  }
  if (type == FieldStatisticsType::Y) {
    return {mesh.global_extent(0), mesh.global_extent(2)};
  }
  return {};
}

std::vector<std::string>
MeshStatistics::get_variable_labels(FieldStatisticsType type) const
{
  if (type == FieldStatisticsType::XY) {
    return {"z"};
  }
  if (type == FieldStatisticsType::Y) {
    return {"x", "z"};
  }
  return {};
}

void MeshStatistics::write_to_tecplot(std::ostream&       out,
                                      FieldStatisticsType type) const
{
  if (type == FieldStatisticsType::XY) {
    write_to_tecplot_xy(out);
  }
  if (type == FieldStatisticsType::Y) {
    write_to_tecplot_y(out);
  }
}

void MeshStatistics::write_to_tecplot_y(std::ostream& out) const
{
  auto const& comm = mesh.comm();
  if (!comm.is_first(1) || !comm.is_first(2)) {
    mpipp::barrier(comm);
    return;
  }

  auto const nx = mesh.global_extent(0);
  auto const nz = mesh.global_extent(2);
  auto const Lx = 2 * Kokkos::numbers::pi_v<Real> / mesh.pex;
  for (int k = 0; k < nz; ++k) {
    for (int i = 0; i < nx; ++i) {
      auto const x = i * (Lx / nx);
      fmt::print(out, "{:.4e} ", x);
    }
    fmt::print(out, "\n");
  }
  for (int k = 0; k < nz; ++k) {
    for (int i = 0; i < nx; ++i) {
      auto const z = gathered_zz(k) * mesh.hbar;
      fmt::print(out, "{:.4e} ", z);
    }
    fmt::print(out, "\n");
  }

  mpipp::barrier(comm);
}

void MeshStatistics::write_to_tecplot_xy(std::ostream& out) const
{
  auto const& comm = mesh.comm();
  if (!comm.is_first(1) || !comm.is_first(2)) {
    mpipp::barrier(comm);
    return;
  }

  auto const nz = mesh.global_extent(2);
  for (int k = 0; k < nz; ++k) {
    auto const z = gathered_zz(k) * mesh.hbar;
    fmt::print(out, "{:.4e} ", z);
  }
  fmt::print(out, "\n");

  mpipp::barrier(comm);
}

template<typename MT>
CurvilinearMeshStatistics<MT>::CurvilinearMeshStatistics(mesh_t const& mesh_)
  : base_t{static_cast<solver::Mesh const&>(mesh_)}
  , mesh{mesh_}
{}

template<typename MT>
void CurvilinearMeshStatistics<MT>::calculate()
{
  base_t::calculate();

  // Calculate mean eta, averaged in the y direction
  auto const& invJ1 = MDView<Real** [1]>(
    mesh.invJ.data(), mesh.invJ.extent(0), mesh.invJ.extent(1));
  MeanVariance<Real, default_memory_pool> eta_statistics(invJ1);

  auto const& comm = mesh.comm();
  eta_statistics.calculate_local();
  eta_statistics.gather(0, comm.axis_comm[1]);
  if (comm.is_first(1)) {
    invJ_mean_y =
      Kokkos::create_mirror(Kokkos::HostSpace(), eta_statistics.mean);
    std::swap(invJ_mean_y, eta_statistics.mean);
  }
}

template<typename MT>
void CurvilinearMeshStatistics<MT>::write_to_tecplot(
  std::ostream&       out,
  FieldStatisticsType type) const
{
  if (type == FieldStatisticsType::XY) {
    base_t::write_to_tecplot_xy(out);
  }
  if (type == FieldStatisticsType::Y) {
    write_to_tecplot_y(out);
  }
}

template<typename MT>
void CurvilinearMeshStatistics<MT>::write_to_tecplot_y(std::ostream& out) const
{
  auto const& comm = mesh.comm();
  if (!comm.is_first(1) || !comm.is_first(2)) {
    mpipp::barrier(comm);
    return;
  }

  auto const nx = mesh.global_extent(0);
  auto const nz = mesh.global_extent(2);
  auto const Lx = 2 * Kokkos::numbers::pi_v<Real> / mesh.pex;
  for (int k = 0; k < nz; ++k) {
    for (int i = 0; i < nx; ++i) {
      auto const x = i * (Lx / nx);
      fmt::print(out, "{:.4e} ", x);
    }
    fmt::print(out, "\n");
  }
  for (int k = 0; k < nz; ++k) {
    for (int i = 0; i < nx; ++i) {
      if constexpr (std::is_same_v<MT, solver::BottomWaveMesh>) {
        auto const z = (gathered_zz(k) - 1) * invJ_mean_y(i, 0) + mesh.hbar;
        fmt::print(out, "{:.4e} ", z);
      } else if constexpr (std::is_same_v<MT, solver::TopWaveMesh>) {
        auto const z = gathered_zz(k) * invJ_mean_y(i, 0) - mesh.hbar;
        fmt::print(out, "{:.4e} ", z);
      } else {
        static_assert(std::is_same_v<MT, solver::BottomWaveMesh>
                        || std::is_same_v<MT, solver::TopWaveMesh>,
                      "Invalid MeshType");
      }
    }
    fmt::print(out, "\n");
  }

  mpipp::barrier(comm);
}

template<typename MT>
void CurvilinearMeshStatistics<MT>::cleanup()
{
  invJ_mean_y = {};
}

template<typename MT>
CurvilinearMeshStatistics<MT>::~CurvilinearMeshStatistics() = default;

// Explicit instantiation
template class CurvilinearMeshStatistics<solver::BottomWaveMesh>;
template class CurvilinearMeshStatistics<solver::TopWaveMesh>;

} // namespace alps::diagnostics
