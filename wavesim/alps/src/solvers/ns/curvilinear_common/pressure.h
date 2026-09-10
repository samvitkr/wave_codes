#pragma once

#include <common/base/logging.h>
#include <common/base/macros.h>
#include <common/container/view_types.h>
#include <common/container/view_utils.h>
#include <common/device/device_traits.h>
#include <common/kokkos_abstraction/exec_policy.h>
#include <common/kokkos_abstraction/pool_space.h>
#include <common/real_type.h>
#include <common/runtime/async_utils.h>
#include <decomp/ghost_cell_exchange.h>
#include <solvers/mesh/mesh.h>
#include <solvers/ns/pressure_bctype.h>
#include <solvers/ns/pressure_eqn_coeffs.h>
#include <solvers/poisson/tridiagonal_solver.h>
#include <spectral/spectral.h>

#include <Kokkos_Profiling_ScopedRegion.hpp>
#include <mpipp/collectives.h>

namespace alps::solver {

struct PressureCurvilinearEqnOptions
{
  Real abs_tol{(Real)1e-6};
  Real rel_tol{(Real)1e-12};
  int  iterations{7};
};

template<typename MeshType>
class PressureCurvilinearEqn
{
 public:
  using solver_val_t =
    std::conditional_t<ALPS_POISSON_MIXED_PRECISION_IR, float, Real>;

  MeshType const* ptr_mesh;

  MDView<solver_val_t***> d, dl, du;

  std::unique_ptr<TridiagonalSolver<solver_val_t>> solver;

  PressureBCType top_bc_type{PressureBCType::NEUMANN};

  explicit PressureCurvilinearEqn(const MeshType& mesh);

  void initialize() const;

  void solve(const HaloView<Real***>&      pp,
             const MDView<Real***>&        div_u,
             Real                          dt,
             PressureCurvilinearEqnOptions options);

  ~PressureCurvilinearEqn();

 private:
  Logger logger;
};

namespace {
template<typename MT>
struct FunctorDzetaPhiDzeta
{
  MDView<Real***> fx;
  MDView<Real***> fy;
  MDView<Real***> fz;

  HaloView<Real const***> p;

  HaloView<Real const*> dzw;
  HaloView<Real const*> zw;
  HaloView<Real const*> dz;
  MDView<Real const**>  exr;
  MDView<Real const**>  eyr;
  MDView<Real const**>  J;

  int z_begin;
  int z_size;
  int stride;

  FunctorDzetaPhiDzeta(MDView<Real***> const&         fx_,
                       MDView<Real***> const&         fy_,
                       MDView<Real***> const&         fz_,
                       HaloView<Real const***> const& p_,
                       MT const&                      mesh,
                       int                            z_begin_,
                       int                            z_end_)
    : fx(fx_)
    , fy(fy_)
    , fz(fz_)
    , p(p_)
    , dzw(mesh.dzw)
    , zw(mesh.zw)
    , dz(mesh.dz)
    , exr(mesh.exr)
    , eyr(mesh.eyr)
    , J(mesh.J)
    , z_begin(z_begin_)
    , z_size([=] {
      if (z_end_ < z_begin_) {
        throw std::runtime_error("z_end < z_begin");
      }
      return z_end_ - z_begin_;
    }())
    , stride(static_cast<int>(p.stride(2)))
  {}

  KOKKOS_FUNCTION void operator()(GridPolicy<>::member_type const& team) const
  {
    int const j        = team.league_rank() / z_size;
    int const k        = team.league_rank() % z_size + z_begin;
    auto      inv_dzk1 = 1 / dz(k - 1);
    auto      inv_dzk  = 1 / dz(k);
    auto      alpha    = 1 / dzw(k - 1);
    auto      ratio    = dzw(k - 1) / 2 * inv_dzk;
    auto      ratio1   = dzw(k - 1) / 2 * inv_dzk1;

    auto const* p_ptr    = &p(0, j, k);
    auto        offset   = static_cast<int>(&fz(0, j, k) - fz.data());
    auto        offset_j = static_cast<int>(&exr(0, j) - exr.data());
    Kokkos::parallel_for(
      Kokkos::TeamThreadRange(team, fx.extent_int(0)), [&](int& i) {
        auto pw1 = itp2node(p_ptr[i], p_ptr[i - stride], ratio1);
        auto pz  = (p_ptr[i] - p_ptr[i - stride]) * inv_dzk1;
        auto pz1 = (p_ptr[i + stride] - p_ptr[i]) * inv_dzk;
        auto pw  = itp2node(p_ptr[i], p_ptr[i + stride], ratio);

        fz.data()[i + offset] = (pz1 - pz) * J(i, j);
        fx.data()[i + offset] =
          (MT::zeta_x(zw(k), exr.data()[i + offset_j]) * pw
           - MT::zeta_x(zw(k - 1), exr.data()[i + offset_j]) * pw1)
          * alpha;
        fy.data()[i + offset] =
          (MT::zeta_y(zw(k), eyr.data()[i + offset_j]) * pw
           - MT::zeta_y(zw(k - 1), eyr.data()[i + offset_j]) * pw1)
          * alpha;
      });
  }
};

template<typename MT>
struct FunctorPressureResidual
{
  MDView<Real***> residual;

  HaloView<Real const***> px;
  HaloView<Real const***> py;
  MDView<Real const***>   pzz;
  MDView<Real const***>   divU;

  HaloView<Real const*> dzw;
  HaloView<Real const*> dz;
  HaloView<Real const*> zw;
  MDView<Real const**>  exr;
  MDView<Real const**>  eyr;
  MDView<Real const**>  J;
  Real                  hbar;
  Real                  dt;

  int z_begin;
  int z_size;
  int stride;

  FunctorPressureResidual(MDView<Real***> const&       residual_,
                          HaloView<Real***> const&     px_,
                          HaloView<Real***> const&     py_,
                          MDView<Real***> const&       pzz_,
                          MDView<Real const***> const& divU_,
                          Real                         dt_,
                          MT const&                    mesh)
    : residual{residual_}
    , px{px_}
    , py{py_}
    , pzz{pzz_}
    , divU{divU_}
    , dzw{mesh.dzw}
    , dz{mesh.dz}
    , zw{mesh.zw}
    , exr{mesh.exr}
    , eyr{mesh.eyr}
    , J{mesh.J}
    , hbar{mesh.hbar}
    , dt{dt_}
    , z_begin{mesh.comm().is_first(2) ? 1 : 0}
    , z_size{[&mesh, nz = mesh.extent(2)] {
      int begin = mesh.comm().is_first(2) ? 1 : 0;
      int end   = mesh.comm().is_last(2) ? nz - 1 : nz;
      if (end < begin) throw std::runtime_error("z_end < z_begin");
      return end - begin;
    }()}
    , stride(static_cast<int>(residual.stride(2)))
  {}

  KOKKOS_FUNCTION void operator()(GridPolicy<>::member_type const& team) const
  {
    int const j        = team.league_rank() / z_size;
    int const k        = team.league_rank() % z_size + z_begin;
    auto      inv_dzk1 = 1 / dz(k - 1);
    auto      inv_dzk  = 1 / dz(k);
    auto      ratio    = dzw(k - 1) / 2 * inv_dzk;
    auto      ratio1   = dzw(k - 1) / 2 * inv_dzk1;
    auto      beta     = dzw(k - 1) * hbar;

    auto offset   = static_cast<int>(&residual(0, j, k) - residual.data());
    auto offset1  = static_cast<int>(&px(0, j, k) - px.data());
    auto offset_j = static_cast<int>(&exr(0, j) - exr.data());
    Kokkos::parallel_for(
      Kokkos::TeamThreadRange(team, residual.extent_int(0)), [&](int& i) {
        auto px1 = itp2node(
          px.data()[i + offset1], px.data()[i + stride + offset1], ratio);
        auto py1 = itp2node(
          py.data()[i + offset1], py.data()[i + stride + offset1], ratio);
        auto px0 = itp2node(
          px.data()[i + offset1], px.data()[i - stride + offset1], ratio1);
        auto py0 = itp2node(
          py.data()[i + offset1], py.data()[i - stride + offset1], ratio1);

        px0 = MT::zeta_x(zw(k), exr.data()[i + offset_j]) * px1
            - MT::zeta_x(zw(k - 1), exr.data()[i + offset_j]) * px0;
        py0 = MT::zeta_y(zw(k), eyr.data()[i + offset_j]) * py1
            - MT::zeta_y(zw(k - 1), eyr.data()[i + offset_j]) * py0;
        auto g = residual.data()[i + offset] * beta
               + (px0 + py0 + J(i, j) * pzz.data()[i + offset]) * hbar;

        residual.data()[i + offset] =
          (-g + (divU.data()[i + offset] / dt) * beta) * beta;
      });
  }

  void execute(Kokkos::DefaultExecutionSpace const& stream) const
  {
    const auto policy = [&, ny = residual.extent_int(1)] {
      if constexpr (is_hip_execution_space_v<decltype(stream)>) {
        return GridPolicy<>(stream, ny * z_size, 128);
      }
      return GridPolicy<>(stream, ny * z_size, Kokkos::AUTO());
    }();
    Kokkos::parallel_for("pressure residual final", policy, *this);
  }
};

template<typename EqnType>
void calc_laplace_p_nonlinear(const MDView<Real***>&               lap_p,
                              const HaloView<Real***>&             p,
                              const MDView<const Real***>&         div_u,
                              const Real                           dt,
                              const EqnType&                       eqn,
                              const Kokkos::DefaultExecutionSpace& space)
{
  const auto& mesh = *eqn.ptr_mesh;
  using MT         = std::decay_t<decltype(mesh)>;

  const auto& grid      = mesh.grid;
  const auto  nx        = local_extent(p, 0);
  const auto  ny        = local_extent(p, 1);
  const auto  nz        = local_extent(p, 2);
  auto const& J         = mesh.J;
  const auto& exr       = mesh.exr;
  const auto& eyr       = mesh.eyr;
  const auto& dz        = mesh.dz;
  const auto  is_top    = grid.comm().is_last(2);
  const auto  is_bottom = grid.comm().is_first(2);
  const auto  top_bc    = eqn.top_bc_type;

  // Temporary arrays
  HaloView<Real***, default_memory_pool> const px(
    Kokkos::view_alloc("px", Kokkos::WithoutInitializing),
    p.layout(),
    {begin(p, 0), begin(p, 1), begin(p, 2)});
  HaloView<Real***, default_memory_pool> const py(
    Kokkos::view_alloc("py", Kokkos::WithoutInitializing),
    p.layout(),
    {begin(p, 0), begin(p, 1), begin(p, 2)});
  MDView<Real***, default_memory_pool> const pzz(
    Kokkos::view_alloc("pzz", Kokkos::WithoutInitializing), lap_p.layout());
  const auto px_inner = create_inner_view(px).view();
  const auto py_inner = create_inner_view(py).view();
  const auto p_inner  = create_inner_view(p).view();

  const auto&    stream1   = space;
  const auto&    stream2   = get_next_stream();
  constexpr auto tile_size = [&]() -> Kokkos::Array<range_index_t, 3> {
    using Device = decltype(stream1);
    if constexpr (is_cuda_execution_space_v<Device>) return {16, 8, 1};
    if constexpr (is_hip_execution_space_v<Device>) return {64, 8, 1};
    return {0, 0, 0};
  }();

  // Update the ghost cells of pressure field
  auto reqs_p = async_update_halo_z(grid.x_pencil(), p, 1);

  /* Calculate ∂(ζ_x Φ)/∂ζ and ∂(ζ_y Φ)/∂ζ at cell centers */
  auto const DzetaPhiDzeta_policy = [&]() {
    if constexpr (is_hip_execution_space_v<Kokkos::DefaultExecutionSpace>) {
      return GridPolicy<>(stream1, ny * (nz - 2), 128);
    }
    return GridPolicy<>(stream1, ny * (nz - 2), Kokkos::AUTO());
  }();
  Kokkos::parallel_for(
    "d(zeta_i Phi)/dzeta",
    DzetaPhiDzeta_policy,
    FunctorDzetaPhiDzeta(px_inner, py_inner, pzz, p, mesh, 1, nz - 1));
  // Wait for the ghost cell update before using the boundary values
  reqs_p.waitall();
  if (!is_bottom) {
    auto policy = GridPolicy<>(stream2, ny, Kokkos::AUTO());
    Kokkos::parallel_for(
      "d(zeta_i Phi)/dzeta 0",
      policy,
      FunctorDzetaPhiDzeta(px_inner, py_inner, pzz, p, mesh, 0, 1));
  }
  if (!is_top) {
    auto policy = GridPolicy<>(stream2, ny, Kokkos::AUTO());
    Kokkos::parallel_for(
      "d(zeta_i Phi)/dzeta nz-1",
      policy,
      FunctorDzetaPhiDzeta(px_inner, py_inner, pzz, p, mesh, nz - 1, nz));
  }
  if (is_top && top_bc == PressureBCType::DIRICHLET) {
    auto const& zz = mesh.zz;
    Kokkos::parallel_for(
      "d(zeta_i Phi)/dzeta",
      LoopPolicy<2>(stream2, {0, 0}, {nx, ny}),
      KOKKOS_LAMBDA(int i, int j) {
        auto dz0    = dz(nz - 2);
        auto dz1    = dz(nz - 3);
        auto alpha  = dz0 / (dz0 + dz1);
        auto coeff2 = (1 + alpha) / dz0;
        auto coeff1 = -1 / alpha / dz1;
        auto coeff0 = alpha / dz1;
        px_inner(i, j, nz - 1) =
          MT::zeta_x(zz(nz - 3), exr(i, j)) * p(i, j, nz - 3) * coeff0
          + MT::zeta_x(zz(nz - 2), exr(i, j)) * p(i, j, nz - 2) * coeff1
          + MT::zeta_x(zz(nz - 1), exr(i, j)) * p(i, j, nz - 1) * coeff2;
        py_inner(i, j, nz - 1) =
          MT::zeta_y(zz(nz - 3), eyr(i, j)) * p(i, j, nz - 3) * coeff0
          + MT::zeta_y(zz(nz - 2), eyr(i, j)) * p(i, j, nz - 2) * coeff1
          + MT::zeta_y(zz(nz - 1), eyr(i, j)) * p(i, j, nz - 1) * coeff2;
      });
  }
  auto event_1 = get_device_event();
  enqueue(event_1, stream2);
  wait_for(
    event_1,
    stream1); // let subsequent computation wait for the above two ops finish

  /* Calculate ∂Φ/∂ξ+∂(ζ_x Φ)/∂ζ */
  spectral::ddx(lap_p, p_inner, grid, stream1);
  Kokkos::parallel_for(
    "d(Phi)/dx",
    LoopPolicy<3>(stream1, {0, 0, 0}, {nx, ny, nz}, tile_size),
    KOKKOS_LAMBDA(int i, int j, int k) {
      px_inner(i, j, k) += lap_p(i, j, k);
    });
  spectral::dealias(px_inner, grid, stream1);
  spectral::ddx(lap_p, px_inner, grid, stream1); // ∂[∂Φ/∂ξ+∂(ζ_x Φ)/∂ζ]/∂ξ
  /* Calculate ∂Φ/∂ψ+∂(ζ_y Φ)/∂ζ */
  stream1.fence(); // ddy and dealias cannot overlap with the above computation
  spectral::ddy_and_add(py_inner, p_inner, grid, stream2);
  spectral::dealias(py_inner, grid, stream2);
  // ∂[∂Φ/∂ψ+∂(ζ_y Φ)/∂ζ]/∂ψ
  spectral::ddy_and_add(lap_p, py_inner, grid, stream2);

  // halo update of px may overlap with above computation
  update_halo_z(grid.x_pencil(), px, 2);
  stream2.fence();
  // For Dirichlet velocity boundary conditions, u=\hat{u} at the boundary imply
  // the pressure corrections at the boundaries are zero.
  if (is_bottom) {
    Kokkos::deep_copy(stream1, subview(px_inner, ALL, ALL, 0), 0);
    Kokkos::deep_copy(stream1, subview(py_inner, ALL, ALL, 0), 0);
  }
  // For Dirichlet pressure BC at the top (free surface), dp/dx and dp/dy are
  // retained.
  if (is_top && top_bc == PressureBCType::NEUMANN) {
    Kokkos::deep_copy(stream1, subview(px_inner, ALL, ALL, nz - 1), 0);
    Kokkos::deep_copy(stream1, subview(py_inner, ALL, ALL, nz - 1), 0);
  }
  // halo update of py may overlap with the subsequent computation
  auto reqs_y = async_update_halo_z(grid.x_pencil(), py, 3);

  /* Calculate ∂(∂(J Φ)/∂ζ)/∂ζ at boundaries */
  auto const policy2d = LoopPolicy<2>(stream2, {0, 0}, {nx, ny});
  if (is_bottom) {
    Kokkos::parallel_for(
      "dPhi/dzeta bottom", policy2d, KOKKOS_LAMBDA(int i, int j) {
        pzz(i, j, 1) = ((p(i, j, 2) - p(i, j, 1)) / dz(1) - 0) * J(i, j);
      });
  }
  if (is_top && top_bc == PressureBCType::NEUMANN) {
    auto const k = nz - 2;
    Kokkos::parallel_for(
      "dPhi/dzeta top", policy2d, KOKKOS_LAMBDA(int i, int j) {
        pzz(i, j, k) =
          (0 - (p(i, j, k) - p(i, j, k - 1)) / dz(k - 1)) * J(i, j);
      });
  }
  if (is_top && top_bc == PressureBCType::DIRICHLET) {
    Kokkos::parallel_for(
      "dPhi/dzeta top", policy2d, KOKKOS_LAMBDA(int i, int j) {
        auto dz0    = dz(nz - 2);
        auto dz1    = dz(nz - 3);
        auto beta   = dz0 / (dz0 + dz1);
        auto coeff2 = (1 + beta) / dz0;
        auto coeff1 = -1 / beta / dz1;
        auto coeff0 = beta / dz1;
        auto pz_top = p(i, j, nz - 3) * coeff0 + p(i, j, nz - 2) * coeff1
                    + p(i, j, nz - 1) * coeff2;
        pzz(i, j, nz - 2) =
          (pz_top - (p(i, j, nz - 2) - p(i, j, nz - 3)) / dz(nz - 3)) * J(i, j);
      });
  }
  /* dealias J∂(∂Φ/∂ζ)/∂ζ */
  spectral::dealias(pzz, grid, stream2);
  enqueue(event_1, stream2);

  reqs_y.waitall(); // wait for the halo update of py finish
  /* Calculate ∂[ζ_x(∂Φ/∂ξ+∂(ζ_x Φ)/∂ζ)]/∂ζ and ∂[ζ_y(∂Φ/∂ψ+∂(ζ_y Φ)/∂ζ)]/∂ζ */
  // calculate only the inner points not involving halo points
  wait_for(event_1, stream1);
  FunctorPressureResidual(lap_p, px, py, pzz, div_u, dt, mesh).execute(stream1);

  stream1.fence();
}

template<class CurveMesh>
void set_pressure_eqn_rhs_bc(const MDView<Real***>&                   sigma,
                             const PressureCurvilinearEqn<CurveMesh>& peqn,
                             const Kokkos::DefaultExecutionSpace&     stream)
{
  auto const& mesh      = *peqn.ptr_mesh;
  const auto& grid      = mesh.grid;
  const auto  is_top    = grid.comm().is_last(2);
  const auto  is_bottom = grid.comm().is_first(2);

  const auto  ends = grid.extents();
  const auto  nz   = ends[2];
  const auto& dz   = mesh.dz;
  if (is_bottom) {
    Kokkos::parallel_for(
      "peqn rhs bottom bc",
      LoopPolicy<2>(stream, {0, 0}, {ends[0], ends[1]}),
      KOKKOS_LAMBDA(int i, int j) {
        auto alpha     = 2 + dz(1) / dz(0);
        sigma(i, j, 0) = sigma(i, j, 1) / alpha;
      });
  }

  if (is_top && peqn.top_bc_type == PressureBCType::NEUMANN) {
    Kokkos::parallel_for(
      "peqn rhs top bc",
      LoopPolicy<2>(stream, {0, 0}, {ends[0], ends[1]}),
      KOKKOS_LAMBDA(int i, int j) {
        auto alpha          = -(2 + dz(nz - 3) / dz(nz - 2));
        sigma(i, j, nz - 1) = sigma(i, j, nz - 2) / alpha;
      });
  }
  if (is_top && peqn.top_bc_type == PressureBCType::DIRICHLET) {
    // Set the top boundary to 0
    Kokkos::deep_copy(stream, subview(sigma, ALL, ALL, nz - 1), 0);
  }
}

template<typename MeshType>
void set_pressure_eqn_rhs(const MDView<Real***>&                  sigma,
                          const HaloView<Real***>&                p,
                          const MDView<const Real***>&            div_u,
                          Real                                    dt,
                          const PressureCurvilinearEqn<MeshType>& peqn,
                          const Kokkos::DefaultExecutionSpace&    stream)
{
  auto const region = Kokkos::Profiling::ScopedRegion("pressure NL residual");

  calc_laplace_p_nonlinear(sigma, p, div_u, dt, peqn, stream);

  set_pressure_eqn_rhs_bc(sigma, peqn, stream);

  spectral::dealias(sigma, peqn.ptr_mesh->grid, stream);
}
} // anonymous namespace

template<typename MeshType>
PressureCurvilinearEqn<MeshType>::PressureCurvilinearEqn(const MeshType& mesh)
  : ptr_mesh{&mesh}
  , d(MDView<solver_val_t***, default_memory_pool>(
      Kokkos::view_alloc("d", Kokkos::WithoutInitializing),
      mesh.grid
        .template get_r2c_xy_output_layout<solver_val_t,
                                           Kokkos::DefaultExecutionSpace>()))
  , dl(MDView<solver_val_t***, default_memory_pool>(
      Kokkos::view_alloc("dl", Kokkos::WithoutInitializing),
      d.layout()))
  , du(MDView<solver_val_t***, default_memory_pool>(
      Kokkos::view_alloc("du", Kokkos::WithoutInitializing),
      d.layout()))
  , solver{create_tridiagonal_solver<solver_val_t>(
      (mesh.global_extent(1) / 3) * 2,
      d.extent_int(1),
      d.extent_int(2),
      mesh.grid.comm().axis_comm[2])}
  , logger{get_logger("p_eqn")}
{}

template<typename MeshType>
void PressureCurvilinearEqn<MeshType>::initialize() const
{
  const auto stream = get_next_stream();

  set_pressure_eqn_coefficients(d, dl, du, top_bc_type, *ptr_mesh, stream);
  stream.fence();

  solver->setup(d, dl, du);
}

template<typename MeshType>
PressureCurvilinearEqn<MeshType>::~PressureCurvilinearEqn() = default;

namespace {
Real estimated_next_residual(std::vector<Real> const& history)
{
  if (history.size() < 2) return history.back();

  auto n  = history.size();
  auto e0 = history[n - 1];
  auto e1 = history[n - 2];

  // assuming {e_n}=lambda * {e_(n-1)}
  return (e0 / e1) * e0;
}

enum class ConvergenceStatus
{
  DIVERGED,
  NOT_CONVERGING,
  STALLED,
  CONVERGED
};

ConvergenceStatus
check_convergence(std::vector<Real> const& history, Real abs_tol, Real rel_tol)
{
  if (history.empty()) return ConvergenceStatus::NOT_CONVERGING;

  auto const last  = history.back();
  auto const first = history.front();

  if (Kokkos::isnan(last) || last < 0 // gpu may produce negative norm
      || last > first * 10) {
    return ConvergenceStatus::DIVERGED;
  }
  if (last <= abs_tol || last < first * rel_tol) {
    return ConvergenceStatus::CONVERGED;
  }

  constexpr double stall_threshold = 0.85;
  constexpr int    stall_window    = 2;
  if (history.size() > stall_window) {
    auto const prev = *(history.crbegin() + stall_window);
    if (double(last / prev) > stall_threshold) {
      return ConvergenceStatus::STALLED;
    }
  }

  return ConvergenceStatus::NOT_CONVERGING;
}

template<typename EqnValT>
struct MPEquation
{
  MDView<Real***>    residual;
  MDView<EqnValT***> rhs;
  MDView<EqnValT***> sigma_y;
  HaloView<Real***>  pp;

  HaloView<Real const*> dzw;
  Real                  hbar;

  int z_begin;
  int z_end;

  MPEquation(HaloView<Real***> const& pp_,
             MDView<Real***> const&   residual_,
             Mesh const&              mesh)
    : residual{residual_}
    , rhs{[&]() -> MDView<EqnValT***> {
      if constexpr (std::is_same_v<EqnValT, Real>) {
        return residual_;
      } else {
        return MDView<EqnValT***, default_memory_pool>(
          Kokkos::view_alloc("p rhs", Kokkos::WithoutInitializing),
          residual_.layout());
      }
      ALPS_UNREACHABLE(decltype(rhs){});
    }()}
    , sigma_y{MDView<EqnValT***, default_memory_pool>(
        Kokkos::view_alloc("sigma_y", Kokkos::WithoutInitializing),
        mesh.grid
          .get_r2c_xy_output_layout<EqnValT, Kokkos::DefaultExecutionSpace>())}
    , pp{pp_}
    , dzw{mesh.dzw}
    , hbar{mesh.hbar}
    , z_begin{mesh.comm().is_first(2) ? 1 : 0}
    , z_end{mesh.comm().is_last(2) ? mesh.extent(2) - 1 : mesh.extent(2)}
  {}

  template<typename EqnType>
  void solve(EqnType const&                       eqn,
             Kokkos::DefaultExecutionSpace const& stream) const
  {
    auto const& grid = eqn.ptr_mesh->grid;

    // if the solver overwrites coefficients, the coefficients need to be reset
    auto const copy_stream = get_next_stream();
    if (eqn.solver->is_coeff_overwritten) {
      set_pressure_eqn_coefficients(
        eqn.d, eqn.dl, eqn.du, eqn.top_bc_type, *eqn.ptr_mesh, copy_stream);
    }

    spectral::fft_r2c_xy(sigma_y, rhs, grid, stream);

    if (eqn.top_bc_type == PressureBCType::NEUMANN) {
      // Set the mean pressure, (kx, ky) = (0, 0)), to zero
      if (grid.comm().is_first(1) && grid.comm().is_last(2)) {
        auto sigma_y_0 = Kokkos::subview(sigma_y, 0, 0, sigma_y.extent(2) - 1);
        Kokkos::deep_copy(stream, sigma_y_0, 0);
      }
    }

    stream.fence();
    copy_stream.fence();
    mpipp::barrier(eqn.solver->comm_);

    eqn.solver->solve(sigma_y, eqn.d, eqn.dl, eqn.du);

    spectral::fft_c2r_xy(rhs, sigma_y, true, grid, stream);
  }

  template<int multiplier>
  KOKKOS_FORCEINLINE_FUNCTION void
  assign_impl(GridPolicy<>::member_type const& team, Real& l_max) const
  {
    int  k     = team.league_rank() / multiplier + z_begin;
    int  j0    = team.league_rank() % multiplier;
    auto coeff = dzw(k - 1) * hbar;
    Real s     = 0;
    Kokkos::parallel_reduce(
      Kokkos::TeamThreadRange(team, rhs.extent_int(1) / multiplier),
      [&](int& j, Real& inner_lmax) {
        Real t_max{};
        Kokkos::parallel_reduce(
          Kokkos::ThreadVectorRange(team, rhs.extent_int(0)),
          [=](int& i, Real& v_max) {
            auto val = residual(i, j * multiplier + j0, k);
            if constexpr (!std::is_same_v<EqnValT, Real>) {
              rhs(i, j * multiplier + j0, k) = static_cast<EqnValT>(val);
            }
            // RHS residual is divided by (dzw*H)^2 to estimate the divergence
            // residual
            v_max = Kokkos::max(Kokkos::abs(val) / coeff / coeff, v_max);
          },
          Kokkos::Max<Real>(t_max));
        inner_lmax = Kokkos::max(t_max, inner_lmax);
      },
      Kokkos::Max<Real>(s));
    l_max = Kokkos::max(l_max, s);
  }

  KOKKOS_FUNCTION void operator()(std::integral_constant<int, 8> /*tag*/,
                                  GridPolicy<>::member_type const& team,
                                  Real&                            l_max) const
  {
    assign_impl<8>(team, l_max);
  }

  KOKKOS_FUNCTION void operator()(std::integral_constant<int, 4> /*tag*/,
                                  GridPolicy<>::member_type const& team,
                                  Real&                            l_max) const
  {
    assign_impl<4>(team, l_max);
  }

  KOKKOS_FUNCTION void operator()(std::integral_constant<int, 2> /*tag*/,
                                  GridPolicy<>::member_type const& team,
                                  Real&                            l_max) const
  {
    assign_impl<2>(team, l_max);
  }

  Real assign_residual_local(Kokkos::DefaultExecutionSpace const& stream) const
  {
    auto constexpr vector_len =
      is_cuda_execution_space_v<Kokkos::DefaultExecutionSpace>  ? 32
      : is_hip_execution_space_v<Kokkos::DefaultExecutionSpace> ? 64
                                                                : 1;
    Real nrm_local{};
    if (rhs.extent_int(1) % 8 == 0) {
      auto const policy = GridPolicy<std::integral_constant<int, 8>>(
        stream, (z_end - z_begin) * 8, Kokkos::AUTO(), vector_len);
      parallel_reduce(policy, *this, Kokkos::Max<Real>(nrm_local));
    } else if (rhs.extent_int(1) % 4 == 0) {
      auto const policy = GridPolicy<std::integral_constant<int, 4>>(
        stream, (z_end - z_begin) * 4, Kokkos::AUTO(), vector_len);
      parallel_reduce(policy, *this, Kokkos::Max<Real>(nrm_local));
    } else if (rhs.extent_int(1) % 2 == 0) {
      auto const policy = GridPolicy<std::integral_constant<int, 2>>(
        stream, (z_end - z_begin) * 2, Kokkos::AUTO(), vector_len);
      parallel_reduce(policy, *this, Kokkos::Max<Real>(nrm_local));
    }
    return nrm_local;
  }

  // assign (and cast to low precision) residual to rhs
  // return the maximum residual
  Real assign_residual(MPIComm3D const&                     comm,
                       Kokkos::DefaultExecutionSpace const& stream) const
  {
    Real nrm_local = assign_residual_local(stream);
    if (comm.is_first(2)) {
      Kokkos::deep_copy(stream,
                        subview(rhs, Kokkos::ALL, Kokkos::ALL, 0),
                        subview(residual, Kokkos::ALL, Kokkos::ALL, 0));
    }
    if (comm.is_last(2)) {
      Kokkos::deep_copy(stream,
                        subview(rhs, Kokkos::ALL, Kokkos::ALL, z_end),
                        subview(residual, Kokkos::ALL, Kokkos::ALL, z_end));
    }

    Real nrm{};
    mpipp::allreduce(nrm_local, nrm, mpipp::max<Real>(), comm);
    stream.fence();
    return nrm;
  }

  void add_correction(Kokkos::DefaultExecutionSpace const& stream) const
  {
    constexpr auto tile_size = []() -> Kokkos::Array<range_index_t, 3> {
      using Device = Kokkos::DefaultExecutionSpace;
      if constexpr (is_cuda_execution_space_v<Device>) return {32, 4, 1};
      if constexpr (is_hip_execution_space_v<Device>) return {64, 8, 1};
      return {0, 0, 0};
    }();
    auto const& dp = this->rhs;
    auto const& p  = this->pp;
    Kokkos::parallel_for(
      "add correction",
      LoopPolicy<3>(stream,
                    {0, 0, z_begin},
                    {rhs.extent_int(0), rhs.extent_int(1), z_end},
                    tile_size),
      KOKKOS_LAMBDA(int i, int j, int k) {
        p(i, j, k) += static_cast<Real>(dp(i, j, k));
      });
  }
};

void extrapolate_zero_neumann_boundary(
  MDView<Real**> const&                f0,
  MDView<Real**> const&                f1,
  MDView<Real**> const&                f2,
  Real const                           dz0,
  Real const                           dz1,
  Kokkos::DefaultExecutionSpace const& stream)
{
  Kokkos::parallel_for(
    LoopPolicy<2>(stream, {0, 0}, {f0.extent(0), f0.extent(1)}),
    KOKKOS_LAMBDA(int i, int j) {
      auto alpha  = dz0 + dz1;
      auto beta   = alpha + dz0;
      auto coeff1 = alpha * alpha / dz1 / beta;
      auto coeff2 = -dz0 * dz0 / dz1 / beta;
      f0(i, j)    = coeff1 * f1(i, j) + coeff2 * f2(i, j);
    });
}
} // namespace

template<class CurveMesh>
void PressureCurvilinearEqn<CurveMesh>::solve(
  const HaloView<Real***>&            pp,
  const MDView<Real***>&              div_u,
  const Real                          dt,
  const PressureCurvilinearEqnOptions options)
{
  const auto& mesh      = *ptr_mesh;
  const auto  is_top    = mesh.comm().is_last(2);
  const auto  is_bottom = mesh.comm().is_first(2);
  const auto  nz        = local_extent(pp, 2);

  const auto residual = MDView<Real***, default_memory_pool>(
    Kokkos::view_alloc("residual", Kokkos::WithoutInitializing),
    div_u.layout());

  MPEquation<solver_val_t> const mp_eqn(pp, residual, mesh);

  auto           stream    = get_next_stream();
  constexpr auto tile_size = []() -> Kokkos::Array<range_index_t, 3> {
    using Device = Kokkos::DefaultExecutionSpace;
    if constexpr (is_cuda_execution_space_v<Device>) return {32, 4, 1};
    if constexpr (is_hip_execution_space_v<Device>) return {64, 8, 1};
    return {0, 0, 0};
  }();

  // Calculate Φ=J^{-1}p
  const auto& invJ = mesh.invJ;
  Kokkos::parallel_for(
    "scale Phi",
    LoopPolicy<3>(stream, local_begins(pp), local_ends(pp), tile_size),
    KOKKOS_LAMBDA(int i, int j, int k) { pp(i, j, k) *= invJ(i, j); });

  //==================
  // Iteration start
  //==================
  auto const rank = mesh.comm().rank();

  std::vector<Real> converge_history{0};
  converge_history.reserve(options.iterations + 1);
  ConvergenceStatus state        = ConvergenceStatus::NOT_CONVERGING;
  Real              nrm_residual = 0;
  int               n_solves     = 0;
  stream.fence();
  for (int it = 0;
       it < options.iterations && state != ConvergenceStatus::CONVERGED;
       ++it) {
    Kokkos::Profiling::ScopedRegion region(fmt::format("p iter {}", it));

    set_pressure_eqn_rhs(residual, pp, div_u, dt, *this, stream);

    nrm_residual            = mp_eqn.assign_residual(mesh.comm(), stream);
    converge_history.back() = nrm_residual;

    if (rank == 0) {
      logger->debug(
        "Iteration {:3d}, ‖r‖={:11.3e}", it, converge_history.back());
    }

    state = check_convergence(
      converge_history, options.abs_tol / dt, options.rel_tol);
    if (state == ConvergenceStatus::DIVERGED) {
      if (rank == 0) {
        logger->error("Pressure solve diverges.");
      }
      throw std::runtime_error("Solution diverging.");
    }
    if (state == ConvergenceStatus::CONVERGED) {
      break;
    }
    if (state == ConvergenceStatus::STALLED) {
      if (rank == 0) {
        logger->warn("Pressure solve have stalled.");
      }
      break;
    }

    mp_eqn.solve(*this, stream);

    mp_eqn.add_correction(stream);
    // Set boundary of pressure to be dp/dz=0
    if (is_top && top_bc_type == PressureBCType::NEUMANN) {
      extrapolate_zero_neumann_boundary(subview(pp, ALL, ALL, nz - 1).view(),
                                        subview(pp, ALL, ALL, nz - 2).view(),
                                        subview(pp, ALL, ALL, nz - 3).view(),
                                        mesh.dz_h(nz - 2),
                                        mesh.dz_h(nz - 3),
                                        stream);
    }
    if (is_bottom) {
      extrapolate_zero_neumann_boundary(subview(pp, ALL, ALL, 0).view(),
                                        subview(pp, ALL, ALL, 1).view(),
                                        subview(pp, ALL, ALL, 2).view(),
                                        mesh.dz_h(0),
                                        mesh.dz_h(1),
                                        stream);
    }
    stream.fence();

    ++n_solves;

    converge_history.push_back(estimated_next_residual(converge_history));
    state = check_convergence(
      converge_history, options.abs_tol / dt, options.rel_tol);
  }

  if (rank != 0) return; // only rank 0 logs the convergence information
  logger->log(state == ConvergenceStatus::CONVERGED ? spdlog::level::info
                                                    : spdlog::level::warn,
              "Solve {} after {} iterations with {}residual {:.4e}.",
              state == ConvergenceStatus::CONVERGED ? "converged"
                                                    : "NOT converging",
              n_solves,
              nrm_residual != converge_history.back() ? "estimated " : "",
              converge_history.back());
}

} // namespace alps::solver
