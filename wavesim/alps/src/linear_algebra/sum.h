//
// Created by xuananqing on 3/29/23.
//

#pragma once

#include <common/base/macros.h>
#include <common/kokkos_abstraction/exec_policy.h>

#include <Kokkos_Core.hpp>

namespace alps {

template<typename T, typename MemSpace>
struct SumXY
{
  using policy_t = GridPolicy<>;
  using member_t = policy_t::member_type;
  Kokkos::View<T const***, Kokkos::LayoutLeft, MemSpace> f;
  Kokkos::View<T*, Kokkos::LayoutLeft, MemSpace>         sum;
  int                                                    nx{};
  int                                                    ny{};

  SumXY(Kokkos::View<T const***, Kokkos::LayoutLeft, MemSpace> const& f_,
        Kokkos::View<T*, Kokkos::LayoutLeft, MemSpace> const&         sum_)
    : f(f_)
    , sum(sum_)
    , nx{f_.extent_int(0)}
    , ny{f_.extent_int(1)}
  {}

  KOKKOS_FUNCTION void operator()(member_t const& team) const
  {
    const auto k = team.league_rank();

    T team_sum{0};
    Kokkos::parallel_reduce(
      Kokkos::TeamThreadRange(team, ny),
      KOKKOS_TR_LAMBDA(int j, T& j_sum) {
        T i_sum{0};
        Kokkos::parallel_reduce(
          Kokkos::ThreadVectorRange(team, nx),
          KOKKOS_TR_LAMBDA(int i, T& threadSum) { threadSum += f(i, j, k); },
          i_sum);
        Kokkos::single(Kokkos::PerThread(team), [&] { j_sum += i_sum; });
      },
      team_sum);

    Kokkos::single(Kokkos::PerTeam(team),
                   [this, k, team_sum] { sum(k) = team_sum; });
  }

  void execute(GridPolicy<> policy) const
  {
    Kokkos::parallel_for("SumXY " + f.label(), policy, *this);
  }
};
} // namespace alps
