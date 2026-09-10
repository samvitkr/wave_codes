#include "dynamic_smagorinsky_scalar_functors.h"

#include <common/base/macros.h>
#include <common/container/matrix_field.h>
#include <common/container/vector_field.h>
#include <common/container/view_utils.h>
#include <common/device/device_traits.h>
#include <solvers/mesh/mesh.h>

namespace alps::solver::detail {

template<typename FT>
Ki_functor<FT>::Ki_functor(HaloView<Real const***> const& f_,
                           Vector3Field<Real***> const&   u_vec,
                           MDView<FT***> const&           f_test_,
                           Vector3Field<FT***> const&     u_test_,
                           Vector3Field<FT***> const&     Ki_)
  : f(f_)
  , u(u_vec.x)
  , v(u_vec.y)
  , w(u_vec.z)
  , f_test(f_test_)
  , u_test(u_test_)
  , Ki(Ki_)
{}

template<typename FT>
KOKKOS_FUNCTION void Ki_functor<FT>::operator()(Ki_Step1_tag const& /*tag*/,
                                                int i,
                                                int j,
                                                int k) const
{
  u_test.x(i, j, k) = static_cast<FT>(u(i, j, k));
  u_test.y(i, j, k) = static_cast<FT>(v(i, j, k));
  Ki.x(i, j, k)     = static_cast<FT>(f(i, j, k) * u(i, j, k));
  Ki.y(i, j, k)     = static_cast<FT>(f(i, j, k) * v(i, j, k));
  const auto w_c    = itp2center(w(i, j, k - 1), w(i, j, k));
  u_test.z(i, j, k) = static_cast<FT>(w_c);
  Ki.z(i, j, k)     = static_cast<FT>(f(i, j, k) * w_c);
  f_test(i, j, k)   = static_cast<FT>(f(i, j, k));
}

template<typename FT>
KOKKOS_FUNCTION void Ki_functor<FT>::operator()(Ki_Step2_tag const& /*tag*/,
                                                int i,
                                                int j,
                                                int k) const
{
  Ki.x(i, j, k) -= f_test(i, j, k) * u_test.x(i, j, k);
  Ki.y(i, j, k) -= f_test(i, j, k) * u_test.y(i, j, k);
  Ki.z(i, j, k) -= f_test(i, j, k) * u_test.z(i, j, k);
}

template<typename FT>
void Ki_functor<FT>::execute(LoopPolicy<3, Ki_Step1_tag> const& policy) const
{
  Kokkos::parallel_for("Ki_fui", policy, *this);
}
template<typename FT>
void Ki_functor<FT>::execute(LoopPolicy<3, Ki_Step2_tag> const& policy) const
{
  Kokkos::parallel_for("Ki", policy, *this);
}

template<typename FT>
KOKKOS_FUNCTION void Xi_functor<FT>::operator()(Xi_Step1_tag const& /*tag*/,
                                                int i,
                                                int j,
                                                int k) const
{
  const auto S1   = theta_i.x(i, j, k);
  const auto S2   = theta_i.y(i, j, k);
  const auto S3_c = itp2center(theta_i.z(i, j, k - 1), theta_i.z(i, j, k));

  auto S                  = S_mag(i, j, k);
  theta_i_test.x(i, j, k) = static_cast<FT>(S1);
  theta_i_test.y(i, j, k) = static_cast<FT>(S2);
  theta_i_test.z(i, j, k) = static_cast<FT>(S3_c);

  Xi.x(i, j, k) = static_cast<FT>(S * S1);
  Xi.y(i, j, k) = static_cast<FT>(S * S2);
  Xi.z(i, j, k) = static_cast<FT>(S * S3_c);
}

template<typename FT>
KOKKOS_FUNCTION void Xi_functor<FT>::operator()(Xi_Step2_tag const& /*tag*/,
                                                int i,
                                                int j,
                                                int k) const
{
  const auto S = S_mag_test(i, j, k);
  Xi.x(i, j, k) -= ratio2 * S * theta_i_test.x(i, j, k);
  Xi.y(i, j, k) -= ratio2 * S * theta_i_test.y(i, j, k);
  Xi.z(i, j, k) -= ratio2 * S * theta_i_test.z(i, j, k);
}

template<typename FT>
Xi_functor<FT>::Xi_functor(Vector3Field<Real***> const&   theta_i_,
                           Vector3Field<FT***> const&     theta_i_test_,
                           Vector3Field<FT***> const&     Xi_,
                           HaloView<Real const***> const& S_mag_,
                           HaloView<FT const***> const&   S_mag_test_,
                           FT                             test_filter_ratio)
  : theta_i{theta_i_}
  , theta_i_test{theta_i_test_}
  , Xi{Xi_}
  , S_mag{create_inner_view(S_mag_).view()}
  , S_mag_test{create_inner_view(S_mag_test_).view()}
  , ratio2{test_filter_ratio * test_filter_ratio}
{}

template<typename FT>
void Xi_functor<FT>::execute(LoopPolicy<3, Xi_Step1_tag> const& policy) const
{
  Kokkos::parallel_for("S*theta_i", policy, *this);
}
template<typename FT>
void Xi_functor<FT>::execute(LoopPolicy<3, Xi_Step2_tag> const& policy) const
{
  Kokkos::parallel_for("Xi", policy, *this);
}

template<typename FT>
KiXi_functor<FT>::KiXi_functor(Vector3Field<FT***> const& Ki_,
                               Vector3Field<FT***> const& Xi_,
                               MDView<Real* [2]> const&   local_sum_,
                               int                        z_begin_,
                               int                        z_end_)
  : Ki{Ki_}
  , Xi{Xi_}
  , local_sum{local_sum_}
  , z_begin{z_begin_}
  , z_end{z_end_}
  , nx{local_extent(Ki_, 0)}
  , ny{local_extent(Xi_, 1)}
{}

template<typename FT>
KOKKOS_FUNCTION void KiXi_functor<FT>::operator()(member_t const& team) const
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
      parallel_reduce(
        Kokkos::ThreadVectorRange(team, nx),
        KOKKOS_TR_LAMBDA(int i, Kokkos::complex<FT>& threadSum) {
          const auto KiXi = Ki.x(i, j, k) * Xi.x(i, j, k)
                          + Ki.y(i, j, k) * Xi.y(i, j, k)
                          + Ki.z(i, j, k) * Xi.z(i, j, k);
          const auto XiXi = Xi.x(i, j, k) * Xi.x(i, j, k)
                          + Xi.y(i, j, k) * Xi.y(i, j, k)
                          + Xi.z(i, j, k) * Xi.z(i, j, k);
          threadSum += Kokkos::complex<FT>(KiXi, XiXi) / nx;
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
void KiXi_functor<FT>::execute(Kokkos::DefaultExecutionSpace const& space) const
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
  Kokkos::parallel_for("local KX XX", policy, *this);
}

#if (ALPS_SGS_USE_LOWER_PRECISION)
template struct Ki_functor<float>;
template struct Xi_functor<float>;
template struct KiXi_functor<float>;
#else
template struct Ki_functor<Real>;
template struct Xi_functor<Real>;
template struct KiXi_functor<Real>;
#endif

} // namespace alps::solver::detail
