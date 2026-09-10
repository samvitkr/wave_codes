#pragma once

#include <common/container/matrix_field.h>
#include <common/container/vector_field.h>
#include <common/kokkos_abstraction/exec_policy.h>
#include <common/real_type.h>

namespace alps::solver::detail {

struct Lij_Step1_tag
{};
struct Lij_Step2_tag
{};

template<typename FT>
struct Lij_functor
{
  HaloView<Real const***> u;
  HaloView<Real const***> v;
  HaloView<Real const***> w;
  MDView<FT***>           u_test;
  MDView<FT***>           v_test;
  MDView<FT***>           w_test;
  MDView<FT***>           Lxx;
  MDView<FT***>           Lxy;
  MDView<FT***>           Lxz;
  MDView<FT***>           Lyy;
  MDView<FT***>           Lyz;
  MDView<FT***>           Lzz;

  KOKKOS_FUNCTION void
  operator()(Lij_Step1_tag const& tag, int i, int j, int k) const;

  KOKKOS_FUNCTION void
  operator()(Lij_Step2_tag const& tag, int i, int j, int k) const;

  Lij_functor(Vector3Field<Real***> const&    u_vec,
              Vector3Field<FT***> const&      u_test_,
              SymmTensor33Field<FT***> const& Lij_);

  void execute(LoopPolicy<3, Lij_Step1_tag> const& policy) const;
  void execute(LoopPolicy<3, Lij_Step2_tag> const& policy) const;
};

struct Mij_Step1_tag
{};
struct Mij_Step2_tag
{};

template<typename FT>
struct Mij_functor
{
  SymmTensor33Field<Real***> Sij;
  MDView<FT***>              St_xx;
  MDView<FT***>              St_xy;
  MDView<FT***>              St_xz;
  MDView<FT***>              St_yy;
  MDView<FT***>              St_yz;
  MDView<FT***>              St_zz;
  MDView<FT***>              Mxx;
  MDView<FT***>              Mxy;
  MDView<FT***>              Mxz;
  MDView<FT***>              Myy;
  MDView<FT***>              Myz;
  MDView<FT***>              Mzz;
  HaloView<Real***>          S_mag;
  FT                         ratio2;

  KOKKOS_FUNCTION void
  operator()(Mij_Step1_tag const& tag, int i, int j, int k) const;

  KOKKOS_FUNCTION void
  operator()(Mij_Step2_tag const& tag, int i, int j, int k) const;

  Mij_functor(SymmTensor33Field<Real***> const& Sij_,
              SymmTensor33Field<FT***> const&   Sij_test_,
              SymmTensor33Field<FT***> const&   Mij_,
              HaloView<Real***> const&          S_mag_,
              FT                                test_filter_ratio);

  void execute(LoopPolicy<3, Mij_Step1_tag> const& policy) const;
  void execute(LoopPolicy<3, Mij_Step2_tag> const& policy) const;
};

template<class FT>
struct LijMij_functor
{
  static constexpr int y_block_size = 16;

  using member_t = GridPolicy<>::member_type;
  SymmTensor33Field<FT***> Lij;
  SymmTensor33Field<FT***> Mij;
  MDView<Real* [2]>        local_sum;
  int                      z_begin;
  int                      z_end;
  int                      nx;
  int                      ny;

  KOKKOS_FUNCTION void operator()(member_t const& team) const;

  LijMij_functor(SymmTensor33Field<FT***> const& Lij_,
                 SymmTensor33Field<FT***> const& Mij_,
                 MDView<Real* [2]> const&        local_sum_,
                 int                             z_begin_,
                 int                             z_end_);

  void execute(Kokkos::DefaultExecutionSpace const& space) const;
};

template<typename FT>
struct SijMag_functor
{
  SymmTensor33Field<Real***> Sij;
  HaloView<Real***>          S_mag;

  KOKKOS_FUNCTION void operator()(int i, int j, int k) const;

  SijMag_functor(SymmTensor33Field<Real***> const& Sij_,
                 HaloView<Real***> const&          S_mag_);

  void execute(LoopPolicy<3> const& policy) const;
};

template<typename FT>
struct Ki_functor
{
  HaloView<Real const***> f;
  HaloView<Real const***> u;
  HaloView<Real const***> v;
  HaloView<Real const***> w;
  HaloView<Real const***> f_test;
  Vector3Field<FT***>     u_test;
  Vector3Field<FT***>     Ki;

  KOKKOS_FUNCTION void
  operator()(Lij_Step1_tag const& tag, int i, int j, int k) const;

  KOKKOS_FUNCTION void
  operator()(Lij_Step2_tag const& tag, int i, int j, int k) const;

  Ki_functor(HaloView<Real const***> const& f_,
             Vector3Field<Real***> const&   u_vec,
             HaloView<Real const***> const& f_test_,
             Vector3Field<FT***> const&     u_test_,
             Vector3Field<FT***> const&     Ki_);

  void execute(LoopPolicy<3, Lij_Step1_tag> const& policy) const;
  void execute(LoopPolicy<3, Lij_Step2_tag> const& policy) const;
};

} // namespace alps::solver::detail
