#include "dynamic_smagorinsky_functors.h"

#include <common/base/macros.h>
#include <common/container/matrix_field.h>
#include <common/container/vector_field.h>
#include <common/container/view_utils.h>
#include <common/device/device_traits.h>
#include <solvers/mesh/mesh.h>

namespace alps::solver::detail {

template<typename T>
KOKKOS_FORCEINLINE_FUNCTION auto Smag(T S11, T S12, T S13, T S22, T S23, T S33)
{
  return Kokkos::sqrt(2
                      * (S11 * S11 + S22 * S22 + S33 * S33
                         + 2 * (S12 * S12 + S13 * S13 + S23 * S23)));
}

template<typename FT>
Lij_functor<FT>::Lij_functor(Vector3Field<Real***> const&    u_vec,
                             Vector3Field<FT***> const&      u_test_,
                             SymmTensor33Field<FT***> const& Lij_)
  : u(u_vec.x)
  , v(u_vec.y)
  , w(u_vec.z)
  , u_test(create_inner_view(u_test_.x).view())
  , v_test(create_inner_view(u_test_.y).view())
  , w_test(create_inner_view(u_test_.z).view())
  , Lxx(create_inner_view(Lij_.xx).view())
  , Lxy(create_inner_view(Lij_.xy).view())
  , Lxz(create_inner_view(Lij_.xz).view())
  , Lyy(create_inner_view(Lij_.yy).view())
  , Lyz(create_inner_view(Lij_.yz).view())
  , Lzz(create_inner_view(Lij_.zz).view())
{
  if (!check_same_layout_and_offset(u, v, w)) {
    throw std::runtime_error("Mismatched extents in u, v, w");
  }
  if (!check_same_layout_and_offset(
        u_test, v_test, w_test, Lxx, Lxy, Lxz, Lyy, Lyz, Lzz)) {
    throw std::runtime_error(
      "Mismatched extents in u_test, v_test, w_test, Lij");
  }
}

template<typename FT>
KOKKOS_FUNCTION void Lij_functor<FT>::operator()(Lij_Step1_tag const& /*tag*/,
                                                 int i,
                                                 int j,
                                                 int k) const
{
  auto offset = &u_test(i, j, k) - u_test.data();
  auto u_ijk  = u(i, j, k);
  auto v_ijk  = v(i, j, k);
  auto w_ijk  = w(i, j, k);
  auto w_ijk1 = w(i, j, k - 1);

  u_test.data()[offset] = static_cast<FT>(u_ijk);
  v_test.data()[offset] = static_cast<FT>(v_ijk);
  Lxx.data()[offset]    = static_cast<FT>(u_ijk * u_ijk);
  Lxy.data()[offset]    = static_cast<FT>(u_ijk * v_ijk);
  Lyy.data()[offset]    = static_cast<FT>(v_ijk * v_ijk);
  const auto w_c        = itp2center(w_ijk1, w_ijk);
  w_test.data()[offset] = static_cast<FT>(w_c);
  Lxz.data()[offset]    = static_cast<FT>(u_ijk * w_c);
  Lyz.data()[offset]    = static_cast<FT>(v_ijk * w_c);
  Lzz.data()[offset]    = static_cast<FT>(w_c * w_c);
}

template<typename FT>
KOKKOS_FUNCTION void Lij_functor<FT>::operator()(Lij_Step2_tag const& /*tag*/,
                                                 int i,
                                                 int j,
                                                 int k) const
{
  auto       offset = &Lxx(i, j, k) - Lxx.data();
  auto const ut     = u_test.data()[offset];
  auto const vt     = v_test.data()[offset];
  auto const wt     = w_test.data()[offset];
  Lxx.data()[offset] -= ut * ut;
  Lxy.data()[offset] -= ut * vt;
  Lyy.data()[offset] -= vt * vt;
  Lxz.data()[offset] -= ut * wt;
  Lyz.data()[offset] -= vt * wt;
  Lzz.data()[offset] -= wt * wt;
}

template<typename FT>
void Lij_functor<FT>::execute(LoopPolicy<3, Lij_Step1_tag> const& policy) const
{
  Kokkos::parallel_for("Lij_uiuj", policy, *this);
}
template<typename FT>
void Lij_functor<FT>::execute(LoopPolicy<3, Lij_Step2_tag> const& policy) const
{
  Kokkos::parallel_for("Lij", policy, *this);
}

template<typename FT>
KOKKOS_FUNCTION void Mij_functor<FT>::operator()(Mij_Step1_tag const& /*tag*/,
                                                 int i,
                                                 int j,
                                                 int k) const
{
  auto offset_hv = &S_mag(i, j, k) - S_mag.data();
  auto offset    = &St_xx(i, j, k) - St_xx.data();

  const auto S11 = Sij.xx.data()[offset_hv];
  const auto S22 = Sij.yy.data()[offset_hv];
  const auto S33 = Sij.zz.data()[offset_hv];

  auto S               = S11 * S11 + S22 * S22 + S33 * S33;
  St_xx.data()[offset] = static_cast<FT>(S11);
  St_yy.data()[offset] = static_cast<FT>(S22);
  St_zz.data()[offset] = static_cast<FT>(S33);

  const auto S12    = Sij.xy.data()[offset_hv];
  auto       stride = S_mag.stride(2);
  const auto S13_c =
    itp2center(Sij.xz.data()[offset_hv - stride], Sij.xz.data()[offset_hv]);
  const auto S23_c =
    itp2center(Sij.yz.data()[offset_hv - stride], Sij.yz.data()[offset_hv]);
  St_xy.data()[offset] = static_cast<FT>(S12);
  St_xz.data()[offset] = static_cast<FT>(S13_c);
  St_yz.data()[offset] = static_cast<FT>(S23_c);

  S += 2 * (S13_c * S13_c + S23_c * S23_c + S12 * S12);
  S = Kokkos::sqrt(2 * S);

  S_mag.data()[offset_hv] = S;
  Mxx.data()[offset]      = static_cast<FT>(S * S11);
  Mxy.data()[offset]      = static_cast<FT>(S * S12);
  Mxz.data()[offset]      = static_cast<FT>(S * S13_c);
  Myy.data()[offset]      = static_cast<FT>(S * S22);
  Myz.data()[offset]      = static_cast<FT>(S * S23_c);
  Mzz.data()[offset]      = static_cast<FT>(S * S33);
}

template<typename FT>
KOKKOS_FUNCTION void Mij_functor<FT>::operator()(Mij_Step2_tag const& /*tag*/,
                                                 int i,
                                                 int j,
                                                 int k) const
{
  auto       offset = &St_xx(i, j, k) - St_xx.data();
  auto const Sxx    = St_xx.data()[offset];
  auto const Sxy    = St_xy.data()[offset];
  auto const Sxz    = St_xz.data()[offset];
  auto const Syy    = St_yy.data()[offset];
  auto const Syz    = St_yz.data()[offset];
  auto const Szz    = St_zz.data()[offset];

  auto S               = Smag(Sxx, Sxy, Sxz, Syy, Syz, Szz);
  St_xx.data()[offset] = S; // store |S| for later use
  S *= ratio2;
  Mxx.data()[offset] -= S * Sxx;
  Mxy.data()[offset] -= S * Sxy;
  Mxz.data()[offset] -= S * Sxz;
  Myy.data()[offset] -= S * Syy;
  Myz.data()[offset] -= S * Syz;
  Mzz.data()[offset] -= S * Szz;
}

template<typename FT>
Mij_functor<FT>::Mij_functor(SymmTensor33Field<Real***> const& Sij_,
                             SymmTensor33Field<FT***> const&   Sij_test_,
                             SymmTensor33Field<FT***> const&   Mij_,
                             HaloView<Real***> const&          S_mag_,
                             FT test_filter_ratio)
  : Sij{Sij_}
  , St_xx{create_inner_view(Sij_test_.xx).view()}
  , St_xy{create_inner_view(Sij_test_.xy).view()}
  , St_xz{create_inner_view(Sij_test_.xz).view()}
  , St_yy{create_inner_view(Sij_test_.yy).view()}
  , St_yz{create_inner_view(Sij_test_.yz).view()}
  , St_zz{create_inner_view(Sij_test_.zz).view()}
  , Mxx{create_inner_view(Mij_.xx).view()}
  , Mxy{create_inner_view(Mij_.xy).view()}
  , Mxz{create_inner_view(Mij_.xz).view()}
  , Myy{create_inner_view(Mij_.yy).view()}
  , Myz{create_inner_view(Mij_.yz).view()}
  , Mzz{create_inner_view(Mij_.zz).view()}
  , S_mag{S_mag_}
  , ratio2{test_filter_ratio * test_filter_ratio}
{
  if (!check_same_layout_and_offset(
        Sij.xx, Sij.xy, Sij.xz, Sij.yy, Sij.yz, Sij.zz, S_mag)) {
    throw std::runtime_error("Mismatched extents in Sij and S_mag");
  }
  if (!check_same_layout_and_offset(St_xx,
                                    St_xy,
                                    St_xz,
                                    St_yy,
                                    St_yz,
                                    St_zz,
                                    Mxx,
                                    Mxy,
                                    Mxz,
                                    Myy,
                                    Myz,
                                    Mzz)) {
    throw std::runtime_error("Mismatched extents in Sij_test and Mij");
  }
}

template<typename FT>
void Mij_functor<FT>::execute(LoopPolicy<3, Mij_Step1_tag> const& policy) const
{
  Kokkos::parallel_for("S*Sij", policy, *this);
}
template<typename FT>
void Mij_functor<FT>::execute(LoopPolicy<3, Mij_Step2_tag> const& policy) const
{
  Kokkos::parallel_for("Mij", policy, *this);
}

template<typename FT>
LijMij_functor<FT>::LijMij_functor(SymmTensor33Field<FT***> const& Lij_,
                                   SymmTensor33Field<FT***> const& Mij_,
                                   MDView<Real* [2]> const&        local_sum_,
                                   int                             z_begin_,
                                   int                             z_end_)
  : Lij{Lij_}
  , Mij{Mij_}
  , local_sum{local_sum_}
  , z_begin{z_begin_}
  , z_end{z_end_}
  , nx{local_extent(Lij_, 0)}
  , ny{local_extent(Lij_, 1)}
{
  if (!check_same_layout_and_offset(Lij.xx,
                                    Lij.xy,
                                    Lij.xz,
                                    Lij.yy,
                                    Lij.yz,
                                    Lij.zz,
                                    Mij.xx,
                                    Mij.xy,
                                    Mij.xz,
                                    Mij.yy,
                                    Mij.yz,
                                    Mij.zz)) {
    throw std::runtime_error("Mismatched extents in Lij, Mij");
  }
}

template<typename FT>
KOKKOS_FUNCTION void LijMij_functor<FT>::operator()(member_t const& team) const
{
  const auto n_y_blocks = (ny + y_block_size - 1) / y_block_size;
  const auto k          = team.league_rank() / n_y_blocks + z_begin;
  const auto y_begin    = (team.league_rank() % n_y_blocks) * y_block_size;
  const auto y_end      = Kokkos::min(y_begin + y_block_size, ny);

  Kokkos::complex<Real> team_sum{0, 0};
  Kokkos::parallel_reduce(
    Kokkos::TeamThreadRange(team, y_begin, y_end),
    KOKKOS_TR_LAMBDA(int j, Kokkos::complex<Real>& j_sum) {
      Kokkos::complex<FT> i_sum{0, 0};
      auto                offset = &Lij.xx(0, j, k) - Lij.xx.data();
      parallel_reduce(
        Kokkos::ThreadVectorRange(team, nx),
        KOKKOS_TR_LAMBDA(int i, Kokkos::complex<FT>& threadSum) {
#define M_AijBij(Aij, Bij)                                      \
  Aij.xx.data()[offset + i] * Bij.xx.data()[offset + i]         \
    + Aij.yy.data()[offset + i] * Bij.yy.data()[offset + i]     \
    + Aij.zz.data()[offset + i] * Bij.zz.data()[offset + i]     \
    + (Aij.xy.data()[offset + i] * Bij.xy.data()[offset + i]    \
       + Aij.xz.data()[offset + i] * Bij.xz.data()[offset + i]  \
       + Aij.yz.data()[offset + i] * Bij.yz.data()[offset + i]) \
        * 2
          const auto LijMij = M_AijBij(Lij, Mij);
          const auto MijMij = M_AijBij(Mij, Mij);
#undef M_AijBij
          threadSum += Kokkos::complex<FT>(LijMij, MijMij) / nx;
        },
        i_sum);
      Kokkos::single(Kokkos::PerThread(team), [&] { j_sum += i_sum / ny; });
    },
    team_sum);

  Kokkos::single(Kokkos::PerTeam(team), [this, k, team_sum] {
    Kokkos::atomic_add(&local_sum(k, 0), team_sum.real());
    Kokkos::atomic_add(&local_sum(k, 1), team_sum.imag());
  });
}

template<typename FT>
void LijMij_functor<FT>::execute(
  Kokkos::DefaultExecutionSpace const& space) const
{
  const auto n_y_blocks = (ny + y_block_size - 1) / y_block_size;
  auto const policy     = [&space, teams = (z_end - z_begin) * n_y_blocks] {
    using Device = Kokkos::DefaultExecutionSpace;
    if constexpr (is_cuda_execution_space_v<Device>) {
      return GridPolicy<>(space, teams, Kokkos::AUTO(), 16);
    }
    if constexpr (is_hip_execution_space_v<Device>) {
      return GridPolicy<>(space, teams, 4, 64);
    }
    return GridPolicy<>(space, teams, Kokkos::AUTO());
  }();
  Kokkos::parallel_for("local LM MM", policy, *this);
}

template<typename FT>
KOKKOS_FUNCTION void SijMag_functor<FT>::operator()(int i, int j, int k) const
{
  const auto S11 = Sij.xx(i, j, k);
  const auto S22 = Sij.yy(i, j, k);
  const auto S33 = Sij.zz(i, j, k);

  auto S = S11 * S11 + S22 * S22 + S33 * S33;

  const auto S12   = Sij.xy(i, j, k);
  const auto S13_c = itp2center(Sij.xz(i, j, k - 1), Sij.xz(i, j, k));
  const auto S23_c = itp2center(Sij.yz(i, j, k - 1), Sij.yz(i, j, k));

  S += 2 * (S13_c * S13_c + S23_c * S23_c + S12 * S12);

  S_mag(i, j, k) = Kokkos::sqrt(2 * S);
}

template<typename FT>
SijMag_functor<FT>::SijMag_functor(SymmTensor33Field<Real***> const& Sij_,
                                   HaloView<Real***> const&          S_mag_)
  : Sij{Sij_}
  , S_mag{S_mag_}
{}

template<typename FT>
void SijMag_functor<FT>::execute(LoopPolicy<3> const& policy) const
{
  Kokkos::parallel_for("SijMag", policy, *this);
}

#if (ALPS_SGS_USE_LOWER_PRECISION)
template struct Lij_functor<float>;
template struct Mij_functor<float>;
template struct LijMij_functor<float>;
template struct SijMag_functor<float>;
#else
template struct Lij_functor<Real>;
template struct Mij_functor<Real>;
template struct LijMij_functor<Real>;
template struct SijMag_functor<Real>;
#endif

} // namespace alps::solver::detail
