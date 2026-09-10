#pragma once

#include <Kokkos_OffsetView.hpp>

#include <type_traits>

namespace alps {

/// @brief Alias for a Kokkos::View with LayoutLeft
template<class DataType, class... Properties>
using MDView = Kokkos::View<DataType, Kokkos::LayoutLeft, Properties...>;

/// @brief Alias for a Kokkos::OffsetView with LayoutLeft, used for blocks with
/// halo points, wherein it is assumed that the total number of points is
/// N+2*n_ghosts, and the indices -n_ghost..-1 and N..N+n_ghost-1 (inclusive)
/// are halo points, and the indices 0..N-1 are the interior points.
template<class DataType, class... Properties>
using HaloView =
  Kokkos::Experimental::OffsetView<DataType, Kokkos::LayoutLeft, Properties...>;

template<class>
struct is_mdview : public std::false_type
{};

template<class D, class... P>
struct is_mdview<MDView<D, P...>> : public std::true_type
{};

template<class D, class... P>
struct is_mdview<const MDView<D, P...>> : public std::true_type
{};

template<class>
struct is_haloview : public std::false_type
{};

template<class D, class... P>
struct is_haloview<HaloView<D, P...>> : public std::true_type
{};

template<class D, class... P>
struct is_haloview<const HaloView<D, P...>> : public std::true_type
{};

template<class T>
inline constexpr bool is_view_or_offset_view_v =
  Kokkos::is_view_v<T> || Kokkos::Experimental::is_offset_view_v<T>;

} // namespace alps
