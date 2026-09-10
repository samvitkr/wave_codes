#pragma once

#include "pencil_plan.h"
#include <transpose/transpose.h>

#include <type_traits>
#include <utility>

namespace alps {

/// Transpose from X-pencil to Y-pencil. (asynchronous)
template<class InType, class OutType, class ExecSpace>
void transpose_xy(OutType&&         out,
                  InType&&          in,
                  const PencilPlan& pencil,
                  const ExecSpace&  space)
{
  static_assert(Kokkos::is_view_v<std::decay_t<InType>>, "in must be a view");
  static_assert(Kokkos::is_view_v<std::decay_t<OutType>>, "out must be a view");

  transpose::transpose(std::forward<OutType>(out),
                       std::forward<InType>(in),
                       pencil.global_extent(0, Pencil::X),
                       pencil.global_extent(1, Pencil::X),
                       transpose::TransposeOpAssign(),
                       pencil.comm().axis_comm[1],
                       space);
}

/// Transpose from Y-pencil to X-pencil. (asynchronous)
template<class InType, class OutType, class ExecSpace>
void transpose_yx(OutType&&         out,
                  InType&&          in,
                  const PencilPlan& pencil,
                  const ExecSpace&  space)
{
  static_assert(Kokkos::is_view_v<std::decay_t<InType>>, "in must be a view");
  static_assert(Kokkos::is_view_v<std::decay_t<OutType>>, "out must be a view");

  transpose::transpose(std::forward<OutType>(out),
                       std::forward<InType>(in),
                       pencil.global_extent(0, Pencil::Y),
                       pencil.global_extent(1, Pencil::Y),
                       transpose::TransposeOpAssign(),
                       pencil.comm().axis_comm[1],
                       space);
}

/// Transpose from Y-pencil to X-pencil and add to the output. (asynchronous)
template<class InType, class OutType, class ExecSpace>
void transpose_yx_and_add(OutType&&         out,
                          InType&&          in,
                          const PencilPlan& pencil,
                          const ExecSpace&  space)
{
  static_assert(Kokkos::is_view_v<std::decay_t<InType>>, "in must be a view");
  static_assert(Kokkos::is_view_v<std::decay_t<OutType>>, "out must be a view");

  transpose::transpose(std::forward<OutType>(out),
                       std::forward<InType>(in),
                       pencil.global_extent(0, Pencil::Y),
                       pencil.global_extent(1, Pencil::Y),
                       transpose::TransposeOpAdd(),
                       pencil.comm().axis_comm[1],
                       space);
}
} // namespace alps
