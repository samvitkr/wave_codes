#pragma once

#include <common/container/view_types.h>
#include <common/container/view_utils.h>
#include <common/kokkos_abstraction/exec_policy.h>

#include <Kokkos_Core.hpp>

#include <utility>

template<typename T, typename ES>
struct TransposeResultVerifier
{
  using mem_space = typename ES::memory_space;

  int                           global_nx_, global_ny_, num_ranks_;
  alps::MDView<T***, mem_space> global_reference_;

  TransposeResultVerifier(int nx, int ny, int nz, int num_ranks)
    : global_nx_(nx)
    , global_ny_(ny)
    , num_ranks_(num_ranks)
    , global_reference_("reference array (global)", nx, ny, nz)
  {}

  void initialize_reference() const
  {
    auto const stream = ES();
    // initialize global data
    auto const& ref = global_reference_;
    Kokkos::parallel_for(
      alps::LoopPolicy<3, ES>(stream, alps::begins(ref), alps::ends(ref)),
      KOKKOS_LAMBDA(int i, int j, int k) {
        ref(i, j, k) = i * 10000 + j * 100 + k;
      });
    stream.fence();
  }

  alps::MDView<T const***, mem_space> get_input_view(int rank) const
  {
    auto block  = global_ny_ / num_ranks_;
    auto offset = rank * block;
    auto input  = alps::MDView<T***, mem_space>(
      "input view (local)", global_nx_, block, global_reference_.extent(2));
    Kokkos::deep_copy(input,
                      subview(global_reference_,
                              Kokkos::ALL,
                              std::pair(offset, offset + block),
                              Kokkos::ALL));
    return input;
  }

  template<typename Op>
  bool compare(alps::MDView<T const***, mem_space> const& result,
               int                                        rank,
               Op const& /*op*/,
               T add_value = 0) const
  {
    auto& ref    = global_reference_;
    auto  offset = rank * (global_nx_ / num_ranks_);
    // check result
    alps::LoopPolicy<3, ES> policy(
      ES(),
      {0, 0, 0},
      {global_ny_,
       global_nx_ / num_ranks_,
       Kokkos::min(result.extent_int(2), global_reference_.extent_int(2))});
    bool OVResult{true};
    Kokkos::parallel_reduce(
      policy,
      KOKKOS_LAMBDA(int i, int j, int k, bool& equal) {
        T before_op_ref = ref(j + offset, i, k);
        T ref_val       = add_value;
        Op::apply(ref_val, before_op_ref);
        equal = equal && (ref_val == result(i, j, k));
      },
      Kokkos::LAnd<bool>(OVResult));

    ES().fence();
    return OVResult;
  }
};
