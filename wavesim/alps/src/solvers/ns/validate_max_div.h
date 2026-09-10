//
// Created by xuanx004 on 6/30/24.
//

#pragma once

#include <common/real_type.h>
#include <solvers/mesh/mesh_fwd.h>

#include <Kokkos_Core.hpp>
#include <mpipp/collectives.h>

#include <array>
#include <utility>

namespace alps::solver {

template<typename FieldType>
std::pair<Real, std::array<int, 3>> validate_max_div(FieldType const& flow)
{
  auto [div_nrm_local, max_div_indices] =
    max_div(flow.u, flow.mesh, CenterPt());
  if (Kokkos::isnan(div_nrm_local) || div_nrm_local < 0) {
    throw std::runtime_error("Result NaN or 0");
  }
  if (div_nrm_local > std::numeric_limits<decltype(div_nrm_local)>::max() / 2) {
    throw std::runtime_error("Solution diverging");
  }

  // Obtain the maximum divergence and the rank where it occurs
  const auto&                comm = flow.mesh.comm();
  const auto                 rank = comm.rank();
  const std::pair<Real, int> local_max_and_rank(div_nrm_local, rank);
  std::pair<Real, int>       div_nrm;
  mpipp::allreduce(local_max_and_rank,
                   div_nrm,
                   mpipp::max_with_location<decltype(div_nrm)>(),
                   comm);

  // Allreduce the indices of the maximum divergence
  if (rank == div_nrm.second) {
    const auto offsets = flow.mesh.partition().offsets;
    max_div_indices[0] += offsets[0];
    max_div_indices[1] += offsets[1];
    max_div_indices[2] += offsets[2];
  }
  mpipp::bcast(nonstd::span(max_div_indices), div_nrm.second, comm);
  return {div_nrm.first, max_div_indices};
}

} // namespace alps::solver
