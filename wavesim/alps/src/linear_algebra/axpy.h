#pragma once

#include <common/container/view_types.h>

#include <Kokkos_Core.hpp>

namespace alps {
namespace detail {
template<class T, class alpha_t, class MemSpace>
struct axpy_kernel
{
  template<class DataType>
  using vector_t =
    MDView<DataType,
           MemSpace,
           Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::Restrict>>;

  axpy_kernel(alpha_t alpha, const vector_t<const T*>& x, const vector_t<T*>& y)
    : alpha_{alpha}
    , x_{x}
    , y_{y}
  {}

  KOKKOS_INLINE_FUNCTION void operator()(int i) const
  {
    y_(i) += alpha_ * x_(i);
  }

  alpha_t            alpha_;
  vector_t<const T*> x_;
  vector_t<T*>       y_;
};
} // namespace detail

template<class policy_t, class alpha_t, class input_t, class output_t>
void axpy_with_policy(alpha_t         alpha,
                      const input_t&  x,
                      const output_t& y,
                      const policy_t& policy)
{
  using T = typename input_t::non_const_value_type;
  using vector_t =
    MDView<T*,
           typename input_t::memory_space,
           Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::Restrict>>;
  vector_t x_(x.data(), x.span());
  vector_t y_(y.data(), y.span());
  Kokkos::parallel_for(
    policy,
    detail::axpy_kernel<T, alpha_t, typename input_t::memory_space>(
      alpha, x_, y_));
}

} // namespace alps
