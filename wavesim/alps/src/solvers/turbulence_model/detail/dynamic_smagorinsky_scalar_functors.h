#pragma once

#include <common/container/matrix_field.h>
#include <common/container/vector_field.h>
#include <common/kokkos_abstraction/exec_policy.h>
#include <common/real_type.h>

namespace alps::solver::detail {

struct Ki_Step1_tag
{};
struct Ki_Step2_tag
{};

template<typename FT>
struct Ki_functor
{
  HaloView<Real const***> f;
  HaloView<Real const***> u;
  HaloView<Real const***> v;
  HaloView<Real const***> w;
  MDView<FT***>           f_test;
  Vector3Field<FT***>     u_test;
  Vector3Field<FT***>     Ki;

  KOKKOS_FUNCTION void
  operator()(Ki_Step1_tag const& tag, int i, int j, int k) const;

  KOKKOS_FUNCTION void
  operator()(Ki_Step2_tag const& tag, int i, int j, int k) const;

  Ki_functor(HaloView<Real const***> const& f_,
             Vector3Field<Real***> const&   u_vec,
             MDView<FT***> const&           f_test_,
             Vector3Field<FT***> const&     u_test_,
             Vector3Field<FT***> const&     Ki_);

  void execute(LoopPolicy<3, Ki_Step1_tag> const& policy) const;
  void execute(LoopPolicy<3, Ki_Step2_tag> const& policy) const;
};

struct Xi_Step1_tag
{};
struct Xi_Step2_tag
{};

template<typename FT>
struct Xi_functor
{
  Vector3Field<Real***> theta_i;
  Vector3Field<FT***>   theta_i_test;
  Vector3Field<FT***>   Xi;
  MDView<Real const***> S_mag;
  MDView<FT const***>   S_mag_test;
  FT                    ratio2;

  KOKKOS_FUNCTION void
  operator()(Xi_Step1_tag const& tag, int i, int j, int k) const;

  KOKKOS_FUNCTION void
  operator()(Xi_Step2_tag const& tag, int i, int j, int k) const;

  Xi_functor(Vector3Field<Real***> const&   theta_i_,
             Vector3Field<FT***> const&     theta_i_test_,
             Vector3Field<FT***> const&     Xi_,
             HaloView<Real const***> const& S_mag_,
             HaloView<FT const***> const&   S_mag_test_,
             FT                             test_filter_ratio);

  void execute(LoopPolicy<3, Xi_Step1_tag> const& policy) const;
  void execute(LoopPolicy<3, Xi_Step2_tag> const& policy) const;
};

template<class FT>
struct KiXi_functor
{
  static constexpr int y_block_size = 16;

  using member_t = GridPolicy<>::member_type;
  Vector3Field<FT***> Ki;
  Vector3Field<FT***> Xi;
  MDView<Real* [2]>   local_sum;
  int                 z_begin;
  int                 z_end;
  int                 nx;
  int                 ny;

  KOKKOS_FUNCTION void operator()(member_t const& team) const;

  KiXi_functor(Vector3Field<FT***> const& Ki_,
               Vector3Field<FT***> const& Xi_,
               MDView<Real* [2]> const&   local_sum_,
               int                        z_begin_,
               int                        z_end_);

  void execute(Kokkos::DefaultExecutionSpace const& space) const;
};

} // namespace alps::solver::detail
