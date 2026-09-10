#include "curvilinear_mesh.h"

#include <common/kokkos_abstraction/exec_policy.h>
#include <common/kokkos_abstraction/pool_space.h>
#include <spectral/spectral.h>

namespace alps::solver {

CurvilinearMesh::CurvilinearMesh(Grid const& grid_, Real Lz, int n_ghost)
  : Mesh(grid_, Lz, n_ghost)
  , ex(MDView<Real**, default_memory_pool>("eta_x", extent(0), extent(1)))
  , ey(MDView<Real**, default_memory_pool>("eta_y", extent(0), extent(1)))
  , J(MDView<Real**, default_memory_pool>("J", extent(0), extent(1)))
  , invJ(MDView<Real**, default_memory_pool>("J^{-1}", extent(0), extent(1)))
  , exr(MDView<Real**, default_memory_pool>("eta_x*J", extent(0), extent(1)))
  , eyr(MDView<Real**, default_memory_pool>("eta_y*J", extent(0), extent(1)))
  , et(MDView<Real**, default_memory_pool>("eta_t", extent(0), extent(1)))
{}

namespace {
template<typename MeshType>
void update_metric_coefficients_impl(MeshType const&                      mesh,
                                     MDView<Real const**> const&          eta,
                                     Kokkos::DefaultExecutionSpace const& space)
{
  static_assert(std::is_same_v<MeshType, BottomWaveMesh>
                  || std::is_same_v<MeshType, TopWaveMesh>,
                "Invalid MeshType");
  MDView<const Real** [1]> const eta_(eta.data(), eta.extent(0), eta.extent(1));

  MDView<Real** [1]> ex1(mesh.ex.data(), mesh.ex.extent(0), mesh.ex.extent(1));
  spectral::ddx(ex1, eta_, mesh.grid, space);

  MDView<Real** [1]> ey1(mesh.ey.data(), mesh.ey.extent(0), mesh.ey.extent(1));
  spectral::ddy(ey1, eta_, mesh.grid, space);

  const auto& invJ_ = mesh.invJ;
  const auto& J_    = mesh.J;
  const auto& ex_   = mesh.ex;
  const auto& ey_   = mesh.ey;
  const auto& exr_  = mesh.exr;
  const auto& eyr_  = mesh.eyr;
  const auto  hbar_ = mesh.hbar;
  auto policy = LoopPolicy<2>(space, {0, 0}, {mesh.extent(0), mesh.extent(1)});

  Kokkos::parallel_for(
    "metric exr eyr", policy, KOKKOS_LAMBDA(int i, int j) {
      invJ_(i, j) = std::is_same_v<MeshType, BottomWaveMesh>
                    ? hbar_ - eta(i, j)
                    : hbar_ + eta(i, j);
      J_(i, j)    = 1 / invJ_(i, j);
      exr_(i, j)  = ex_(i, j) / invJ_(i, j);
      eyr_(i, j)  = ey_(i, j) / invJ_(i, j);
    });

  space.fence();
}
} // anonymous namespace

void BottomWaveMesh::update_metric_coefficients(
  MDView<Real const**> const&          eta,
  Kokkos::DefaultExecutionSpace const& space) const
{
  update_metric_coefficients_impl(*this, eta, space);
}

void TopWaveMesh::update_metric_coefficients(
  MDView<Real const**> const&          eta,
  Kokkos::DefaultExecutionSpace const& space) const
{
  update_metric_coefficients_impl(*this, eta, space);
}

} // namespace alps::solver
