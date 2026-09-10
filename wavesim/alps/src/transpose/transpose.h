#pragma once

#include "transposer_base.h"
#include "transposer_pool.h"

#include <Kokkos_Core.hpp>

namespace alps::transpose {

template<typename InView, typename OutView, typename ExecSpace>
void transpose(OutView&&                  output,
               InView&&                   input,
               int                        input_n0_global,
               int                        input_n1_global,
               TransposeOps               op,
               mpipp::communicator const& comm,
               ExecSpace const&           space)
{
  static_assert(Kokkos::is_view_v<std::decay_t<InView>>,
                "input must be a view");
  static_assert(Kokkos::is_view_v<std::decay_t<OutView>>,
                "output must be a view");

  using T    = typename std::decay_t<InView>::non_const_value_type;
  int   nz   = input.extent_int(2); // only for plan creation
  auto& pool = TransposerPool<T, ExecSpace>::get_instance();

  auto& transposer =
    pool.get_transposer(comm, input_n0_global, input_n1_global, nz);
  transposer.execute(
    std::forward<OutView>(output), std::forward<InView>(input), op, space);
}

} // namespace alps::transpose
