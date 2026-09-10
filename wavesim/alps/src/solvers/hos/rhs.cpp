#include "rhs.h"

#include <common/base/macros.h>
#include <common/container/view_utils.h>
#include <common/kokkos_abstraction/exec_policy.h>
#include <common/kokkos_abstraction/pool_space.h>
#include <common/math.h>
#include <common/runtime/async_utils.h>
#include <spectral/spectral.h>

#include <Kokkos_Core.hpp>

#include <limits>

namespace alps::solver::hos {

namespace detail {

template<bool surface_tension_enabled>
struct dfdt_rhs_functor
{
  using member_t = GridPolicy<>::member_type;

  MDView<Real** [2]>       F_t;
  MDView<Real const**>     eta;
  MDView<Real const** [2]> F_x;
  MDView<Real const** [2]> F_y;
  MDView<Real const** [2]> tmp;
  MDView<Real const**>     ws;
  MDView<Real const**>     pa;
  MDView<Real const**>     p_st;
  Real                     Fr2;

  int nx;

  KOKKOS_FUNCTION void operator()(member_t const& team) const
  {
    const int j = team.league_rank();
    Kokkos::parallel_for(
      Kokkos::TeamVectorRange(team, nx), KOKKOS_TR_LAMBDA(int i) {
        // eta_t
        F_t(i, j, 0) =
          -(F_x(i, j, 1) * F_x(i, j, 0) + F_y(i, j, 1) * F_y(i, j, 0))
          + tmp(i, j, 0) * ws(i, j);
        // vps_t
        auto p_surface = pa(i, j) + eta(i, j) / Fr2;
        if constexpr (surface_tension_enabled) {
          p_surface += p_st(i, j);
        }
        auto uu = F_x(i, j, 1) * F_x(i, j, 1) + F_y(i, j, 1) * F_y(i, j, 1)
                - tmp(i, j, 0) * tmp(i, j, 1);
        F_t(i, j, 1) = -(p_surface + uu / 2);
      });
  }

  dfdt_rhs_functor(MDView<Real** [2]>       F_t_,
                   MDView<Real const** [2]> F_,
                   MDView<Real const** [2]> F_x_,
                   MDView<Real const** [2]> F_y_,
                   MDView<Real const** [2]> tmp_,
                   MDView<Real const**>     ws_,
                   MDView<Real const**>     pa_,
                   MDView<Real const**>     p_st_,
                   Real                     Fr2_)
    : F_t{F_t_}
    , eta{subview(F_, ALL, ALL, 0)}
    , F_x{F_x_}
    , F_y{F_y_}
    , tmp{tmp_}
    , ws{ws_}
    , pa{pa_}
    , p_st{p_st_}
    , Fr2{Fr2_}
    , nx{local_extent(F_, 0)}
  {}
};

} // namespace detail

void calc_evolution_rhs(MDView<Real** [2]>                   F_t,
                        MDView<const Real** [2]>             F,
                        MDView<const Real**>                 ws,
                        MDView<const Real**>                 pa,
                        Real                                 Fr2,
                        Real                                 We,
                        Grid const&                          grid,
                        Kokkos::DefaultExecutionSpace const& space)
{
  auto ends = local_extents(F);

  auto         stream2 = get_next_stream();
  GridPolicy<> policy2(stream2, ends[1], Kokkos::AUTO);

  MDView<Real** [2][3], default_memory_pool> tmp0(
    Kokkos::view_alloc("tmp", Kokkos::WithoutInitializing),
    F.extent(0),
    F.extent(1));
  auto F_x = subview(tmp0, ALL, ALL, ALL, 0);
  auto F_y = subview(tmp0, ALL, ALL, ALL, 1);
  auto tmp = subview(tmp0, ALL, ALL, ALL, 2);

  auto event = get_device_event();
  spectral::ddx(F_x, F, grid, stream2);
  enqueue(event, stream2);
  spectral::ddy(F_y, F, grid, space);
  wait_for(event, space); // sync ddx with space

  using Kokkos::parallel_for;
  GridPolicy<> policy1(space, ends[1], Kokkos::AUTO);
  const int    nx = ends[0];
  parallel_for(
    "scaling", policy1, KOKKOS_LAMBDA(decltype(policy1)::member_type team) {
      const int j = team.league_rank();
      parallel_for(
        Kokkos::TeamVectorRange(team, nx), KOKKOS_TR_LAMBDA(int i) {
          tmp(i, j, 0) =
            1 + F_x(i, j, 0) * F_x(i, j, 0) + F_y(i, j, 0) * F_y(i, j, 0);
        });
    });
  parallel_for(
    "ws^2", policy2, KOKKOS_LAMBDA(decltype(policy2)::member_type team) {
      const int j = team.league_rank();
      Kokkos::parallel_for(
        Kokkos::TeamVectorRange(team, nx),
        KOKKOS_TR_LAMBDA(int i) { tmp(i, j, 1) = ws(i, j) * ws(i, j); });
    });
  enqueue(event, stream2);
  wait_for(event, space);
  spectral::dealias(tmp, grid, space);

  if (We > std::numeric_limits<Real>::epsilon() && !std::isinf(We)) {
    MDView<Real**, default_memory_pool> p_st(
      Kokkos::view_alloc("p tension", Kokkos::WithoutInitializing),
      pa.layout());
    calc_surface_tension(p_st,
                         subview(F_x, ALL, ALL, index_range(0, 1)),
                         subview(F_y, ALL, ALL, index_range(0, 1)),
                         We,
                         grid,
                         space);
    parallel_for(
      "dFdt rhs",
      policy1,
      detail::dfdt_rhs_functor<true>(F_t, F, F_x, F_y, tmp, ws, pa, p_st, Fr2));
    space.fence(); // wait before p_st is deallocated
  } else {
    parallel_for(
      "dFdt rhs",
      policy1,
      detail::dfdt_rhs_functor<false>(F_t, F, F_x, F_y, tmp, ws, pa, {}, Fr2));
    space.fence();
  }
}

void calc_surface_tension(MDView<Real**>                       p_st,
                          MDView<const Real** [1]>             eta_x,
                          MDView<const Real** [1]>             eta_y,
                          Real const                           We,
                          Grid const&                          grid,
                          Kokkos::DefaultExecutionSpace const& space)
{
  MDView<Real** [3], default_memory_pool> tmp(
    Kokkos::view_alloc("eta_2nd", Kokkos::WithoutInitializing),
    eta_x.extent(0),
    eta_x.extent(1));

  auto event   = get_device_event();
  auto stream2 = get_next_stream();
  enqueue(event, space);
  wait_for(event, stream2); // sync stream2 with external stream

  spectral::ddx(subview(tmp, ALL, ALL, index_range(0, 1)), eta_x, grid, space);
  spectral::ddy(
    subview(tmp, ALL, ALL, index_range(1, 2)), eta_y, grid, stream2);
  enqueue(event, stream2);
  spectral::ddx(subview(tmp, ALL, ALL, index_range(2, 3)), eta_y, grid, space);

  wait_for(event, space); // sync ddy with external stream
  LoopPolicy<2> policy(space, {0, 0}, {eta_x.extent(0), eta_x.extent(1)});
  Kokkos::parallel_for(
    "surface tension", policy, KOKKOS_LAMBDA(int i, int j) {
      auto s =
        1 + eta_x(i, j, 0) * eta_x(i, j, 0) + eta_y(i, j, 0) * eta_y(i, j, 0);
      auto t = -(2 / We)
             * ((1 + eta_y(i, j, 0) * eta_y(i, j, 0)) * tmp(i, j, 0)
                - 2 * eta_x(i, j, 0) * eta_y(i, j, 0) * tmp(i, j, 2)
                + (1 + eta_x(i, j, 0) * eta_x(i, j, 0)) * tmp(i, j, 1))
             / 2;
      p_st(i, j) = t * Kokkos::rsqrt(s * s * s);
    });

  space.fence();
}

} // namespace alps::solver::hos
