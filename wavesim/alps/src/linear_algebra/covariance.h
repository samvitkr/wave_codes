//
// Created by xuananqing on 4/4/23.
//

#pragma once

#include <common/kokkos_abstraction/exec_policy.h>

#include <Kokkos_Core.hpp>
#include <mpipp/comm.h>

namespace alps {

/// A class for computing the mean and variance of a distributed 3D field along
/// the y-direction (second axis)
template<typename T, typename MemSpace>
struct Covariance
{
  using policy_t = GridPolicy<>;
  using member_t = policy_t::member_type;

  Kokkos::View<T const***, Kokkos::LayoutLeft, MemSpace> f;
  Kokkos::View<T const***, Kokkos::LayoutLeft, MemSpace> g;

  Kokkos::View<T**, Kokkos::LayoutLeft, Kokkos::HostSpace> mean_f;
  Kokkos::View<T**, Kokkos::LayoutLeft, Kokkos::HostSpace> mean_g;
  Kokkos::View<T**, Kokkos::LayoutLeft, Kokkos::HostSpace> M2;

  Kokkos::View<T**, Kokkos::LayoutLeft, MemSpace> mean_f_d;
  Kokkos::View<T**, Kokkos::LayoutLeft, MemSpace> mean_g_d;
  Kokkos::View<T**, Kokkos::LayoutLeft, MemSpace> M2_d;

  Covariance(Kokkos::View<T const***, Kokkos::LayoutLeft, MemSpace> const& f_,
             Kokkos::View<T const***, Kokkos::LayoutLeft, MemSpace> const& g_);

  Covariance(
    Kokkos::View<T**, Kokkos::LayoutLeft, Kokkos::HostSpace> const& cov_,
    Kokkos::View<T const***, Kokkos::LayoutLeft, MemSpace> const&   f_,
    Kokkos::View<T const***, Kokkos::LayoutLeft, MemSpace> const&   g_);

  KOKKOS_FUNCTION void operator()(member_t const& team) const;

  void calculate_local();

  void gather(int root, const mpipp::communicator& comm);
};

} // namespace alps
