#pragma once

#include <common/container/view_types.h>
#include <common/real_type.h>
#include <io/hdf5.h>
#include <spectral/spectral_fwd.h>

#include <filesystem>
#include <utility>

namespace alps::solver::hos {

struct HOSField
{
  HOSField(Grid const& spectral_grid);

  MDView<Real** [2]> value;
  MDView<Real**>     pa;

  Grid const& grid;

  mutable double time;

  Real kx0;
  Real ky0;

  /// get surface elevation
  auto eta() const
  {
    return Kokkos::subview(value, Kokkos::ALL, Kokkos::ALL, 0);
  }

  /// get surface velocity potential
  auto vps() const
  {
    return Kokkos::subview(value, Kokkos::ALL, Kokkos::ALL, 1);
  }

  /// get surface elevation with a 3rd dimension
  auto eta_3() const
  {
    return Kokkos::subview(value, Kokkos::ALL, Kokkos::ALL, std::pair(0, 1));
  }

  /// get surface velocity potential with a 3rd dimension
  auto vps_3() const
  {
    return Kokkos::subview(value, Kokkos::ALL, Kokkos::ALL, std::pair(1, 2));
  }

  void save(HighFive::File& h5_file) const;

  void save(std::filesystem::path filename) const;
};

} // namespace alps::solver::hos
