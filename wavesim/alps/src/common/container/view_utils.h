#pragma once

#include "view_types.h"

#include <Kokkos_Core.hpp>
#include <Kokkos_OffsetView.hpp>
#include <Kokkos_Pair.hpp>

namespace alps {

using Kokkos::ALL; // NOLINT

template<class I, typename = std::enable_if_t<std::is_integral_v<I>>>
constexpr auto index_range(I i1, I i2)
{
  return Kokkos::make_pair(i1, i2);
}

using range_index_t = std::int64_t;

namespace detail {
template<class ViewType, unsigned... Is>
Kokkos::Array<int, ViewType::rank>
extents_helper(ViewType const& view,
               std::integer_sequence<unsigned, Is...> /*unused*/)
{
  return {static_cast<int>(view.extent(Is))...};
}

template<class ViewType, size_t... Is>
Kokkos::Array<range_index_t, ViewType::rank>
ends_helper(ViewType const& view, std::index_sequence<Is...> /*unused*/)
{
  if constexpr (Kokkos::Experimental::is_offset_view_v<ViewType>) {
    return {(view.begin(Is) + static_cast<range_index_t>(view.extent(Is)))...};
  }
  return {static_cast<range_index_t>(view.extent(Is))...};
}

template<class OutType, class ViewType, size_t... Is>
Kokkos::Array<OutType, ViewType::rank>
local_ends_helper(ViewType const& view, std::index_sequence<Is...> /*unused*/)
{
  if constexpr (Kokkos::Experimental::is_offset_view_v<ViewType>) {
    return {(static_cast<OutType>(view.begin(Is) * 2)
             + static_cast<OutType>(view.extent(Is)))...};
  }
  return {static_cast<OutType>(view.extent(Is))...};
}
} // namespace detail

template<class ViewType>
std::enable_if_t<is_view_or_offset_view_v<ViewType>,
                 Kokkos::Array<int, ViewType::rank>>
extents(const ViewType& view)
{
  return detail::extents_helper(
    view, std::make_integer_sequence<unsigned, ViewType::rank>{});
}

template<class ViewType>
std::enable_if_t<is_view_or_offset_view_v<ViewType>, int>
extent(const ViewType& view, int dim)
{
  return static_cast<int>(view.extent(dim));
}

template<class ViewType>
constexpr std::enable_if_t<is_view_or_offset_view_v<ViewType>, int>
begin(const ViewType& view, int dim)
{
  if constexpr (Kokkos::Experimental::is_offset_view<ViewType>::value) {
    return int(view.begin(dim));
  } else {
    (void)dim;
    return 0;
  }
  return 0;
}

template<class ViewType>
constexpr std::enable_if_t<is_view_or_offset_view_v<ViewType>,
                           Kokkos::Array<range_index_t, ViewType::rank>>
begins(const ViewType& view)
{
  if constexpr (Kokkos::Experimental::is_offset_view_v<ViewType>) {
    return view.begins();
  }
  return {0};
}

template<class ViewType>
std::enable_if_t<is_view_or_offset_view_v<ViewType>, int>
end(const ViewType& view, int dim)
{
  if constexpr (Kokkos::Experimental::is_offset_view<ViewType>::value) {
    return static_cast<int>(view.begin(dim) + view.extent(dim));
  } else {
    return static_cast<int>(view.extent(dim));
  }
  return 0;
}

template<class ViewType>
std::enable_if_t<is_view_or_offset_view_v<ViewType>,
                 Kokkos::Array<range_index_t, ViewType::rank>>
ends(const ViewType& view)
{
  return detail::ends_helper(view, std::make_index_sequence<ViewType::rank>{});
}

template<class ViewType>
std::enable_if_t<is_view_or_offset_view_v<ViewType>, int>
local_extent(const ViewType& view, int dim)
{
  if constexpr (Kokkos::Experimental::is_offset_view<ViewType>::value) {
    return static_cast<int>(view.extent(dim) + view.begin(dim) * 2);
  } else {
    return static_cast<int>(view.extent(dim));
  }
  return 0;
}

template<class ViewType>
std::enable_if_t<is_view_or_offset_view_v<ViewType>,
                 Kokkos::Array<int, ViewType::rank>>
local_extents(const ViewType& view)
{
  return detail::local_ends_helper<int>(
    view, std::make_index_sequence<ViewType::rank>{});
}

template<class ViewType>
constexpr std::enable_if_t<is_view_or_offset_view_v<ViewType>,
                           Kokkos::Array<range_index_t, ViewType::rank>>
local_begins(const ViewType& /*view*/)
{
  return {0};
}

template<class ViewType>
constexpr std::enable_if_t<is_view_or_offset_view_v<ViewType>, int>
local_begin(const ViewType& /*view*/, int /*dim*/)
{
  return 0;
}

template<class ViewType>
std::enable_if_t<is_view_or_offset_view_v<ViewType>, int>
local_end(const ViewType& view, int dim)
{
  return local_extent(view, dim);
}

template<class ViewType>
std::enable_if_t<is_view_or_offset_view_v<ViewType>,
                 Kokkos::Array<range_index_t, ViewType::rank>>
local_ends(const ViewType& view)
{
  return detail::local_ends_helper<range_index_t>(
    view, std::make_index_sequence<ViewType::rank>{});
}

namespace detail {
template<class ViewType, unsigned... Is>
decltype(auto)
create_inner_view_helper(ViewType const& v,
                         std::integer_sequence<unsigned, Is...> /*unused*/)
{
  return Kokkos::Experimental::subview(
    v, index_range((int64_t)0, (int64_t)v.extent(Is) + v.begin(Is) * 2)...);
}

template<class ViewType, std::size_t... Is, class... Args>
decltype(auto) trailing_subview_helper(ViewType const& v,
                                       std::index_sequence<Is...> /*unused*/,
                                       Args... args)
{
  return subview(v, ((void)Is, Kokkos::ALL)..., args...);
}

template<typename ViewType, typename... CP, size_t... Is>
auto create_haloview_helper(
  Kokkos::Impl::ViewCtorProp<CP...> const&  arg_prop,
  Kokkos::Array<int, ViewType::rank> const& local_extents,
  Kokkos::Array<int, ViewType::rank> const& n_ghosts,
  std::index_sequence<Is...> /*unused*/)
{
  return ViewType(arg_prop,
                  typename ViewType::traits::array_layout(
                    (local_extents[Is] + n_ghosts[Is] * 2)...),
                  {-n_ghosts[Is]...});
}
} // namespace detail

/**
 * @brief Create a HaloView with allocation properties (created with
 * Kokkos::view_alloc or Kokkos::view_wrap), local extents, and number of ghost
 * cells. The resulting HaloView will have index ranges [ -n_ghosts[d],
 * local_extents[d] + n_ghosts[d] ] in each dimension d.
 */
template<typename D, typename... P, typename... CP>
HaloView<D, P...>
create_haloview(Kokkos::Impl::ViewCtorProp<CP...> const&    arg_prop,
                Kokkos::Array<int, HaloView<D, P...>::rank> local_extents,
                Kokkos::Array<int, HaloView<D, P...>::rank> n_ghosts)
{
  return detail::create_haloview_helper<HaloView<D, P...>>(
    arg_prop,
    local_extents,
    n_ghosts,
    std::make_index_sequence<HaloView<D, P...>::rank>{});
}

template<class D, class... P>
Kokkos::Experimental::OffsetView<D, Kokkos::LayoutLeft, P...>
create_inner_view(const HaloView<D, P...>& src)
{
  if (!src.is_allocated()) {
    return src;
  }
  auto constexpr rank = Kokkos::Experimental::OffsetView<D, P...>::rank;
  return detail::create_inner_view_helper(
    src, std::make_integer_sequence<unsigned, rank>{});
}

template<class ViewType,
         class... Args,
         typename = std::enable_if_t<is_view_or_offset_view_v<ViewType>>>
decltype(auto) trailing_subview(ViewType const& v, Args... args)
{
  static_assert(sizeof...(Args) <= ViewType::rank,
                "Number of subview arguments exceeds the rank");
  static_assert(sizeof...(Args) >= 1, "At least one argument is required");

  return detail::trailing_subview_helper(
    v, std::make_index_sequence<ViewType::rank - sizeof...(Args)>(), args...);
}

namespace detail {
template<typename T, unsigned... Is, typename... Args>
bool check_same_layout_and_offset_helper(
  std::integer_sequence<unsigned, Is...> /*unused*/,
  T const& first,
  Args const&... args)
{
  auto check_extent = [&](int dim) {
    return ((first.extent(dim) == args.extent(dim)) && ...);
  };
  auto check_stride = [&](int dim) {
    return ((first.stride(dim) == args.stride(dim)) && ...);
  };
  auto check_begin = [&](int dim) {
    return ((begin(first, dim) == begin(args, dim)) && ...);
  };

  return (check_extent(Is) && ...) && (check_stride(Is) && ...)
      && (check_begin(Is) && ...);
}
} // namespace detail

/**
 * @brief Checks if the given HaloViews have the same layout and offset.
 *
 * This function compares the extent, stride, and begin values of the HaloViews.
 */
template<typename T, typename... Args>
bool check_same_layout_and_offset(T const& first, Args const&... args)
{
  static_assert(sizeof...(Args) >= 1,
                "At least one additional argument is required");
  static_assert(is_view_or_offset_view_v<T>
                  && (is_view_or_offset_view_v<Args> && ...),
                "All arguments must be views or offset views");
  static_assert(((T::rank == Args::rank) && ...),
                "All arguments must have the same rank");

  return detail::check_same_layout_and_offset_helper(
    std::make_integer_sequence<unsigned, T::rank>{}, first, args...);
}
} // namespace alps
