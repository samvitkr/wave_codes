#pragma once

#include "view_types.h"
#include "view_utils.h"

#include <Kokkos_Core.hpp>

#include <string>

namespace alps {

template<typename DataType, typename... Properties>
class Vector3Field;

namespace detail {
template<typename... P, typename Label>
auto prepend_ctor_prop_label(
  Kokkos::Impl::ViewCtorProp<P...> const& ctor_prop,
  Label const&                            suffix,
  std::enable_if_t<Kokkos::Impl::is_view_label<Label>::value, int>* /*tag*/ = 0)
{
  static_assert(Kokkos::Impl::ViewCtorProp<P...>::has_label,
                "ViewCtorProp must have a label to prepend suffix");

  auto old_label =
    Kokkos::Impl::get_property<Kokkos::Impl::LabelTag>(ctor_prop);
  Kokkos::Impl::ViewCtorProp<P...> new_ctor_prop(ctor_prop);
  static_cast<Kokkos::Impl::ViewCtorProp<void, std::string>&>(new_ctor_prop)
    .value = old_label + suffix;
  return new_ctor_prop;
}
} // namespace detail

/// @brief 3D field of 3D vectors
/**
 * The vector field is stored in SoA format, and its components are accesses by
 * the member x, y, z.
 */
template<typename DataType, typename... Properties>
class Vector3Field
{
 public:
  using view_t = HaloView<DataType, Properties...>;

  using traits = typename view_t::traits;

  using data_type           = typename view_t::data_type;
  using const_data_type     = typename view_t::const_data_type;
  using non_const_data_type = typename view_t::non_const_data_type;

  using value_type           = typename view_t::value_type;
  using const_value_type     = typename view_t::const_value_type;
  using non_const_value_type = typename view_t::non_const_value_type;

  using array_layout = typename view_t::array_layout;
  using dimension    = typename view_t::dimension;

  using specialize = typename view_t::specialize;

  static constexpr unsigned rank         = dimension::rank;
  static constexpr unsigned rank_dynamic = dimension::rank_dynamic;

  using device_type       = typename view_t::device_type;
  using execution_space   = typename view_t::execution_space;
  using host_mirror_space = typename view_t::host_mirror_space;
  using memory_space      = typename view_t::memory_space;
  using memory_traits     = typename view_t::memory_traits;

  using begins_type = typename view_t::begins_type;

  view_t x, y, z;

  /// @brief Construct an empty Vector3Field
  Vector3Field() = default;

  /**
   * @brief Construct a Vector3Field with the given label, size of the local
   * region, and number of ghost cells
   *
   * @param label The label of the field
   * @param local_extents The size of the local region
   * @param n_ghost The number of ghost cells on each side of the axes
   */
  template<
    typename Label,
    typename = std::enable_if_t<Kokkos::Impl::is_view_label<Label>::value>>
  Vector3Field(Label const&             label,
               Kokkos::Array<int, rank> local_extents,
               Kokkos::Array<int, rank> n_ghost = {})
    : Vector3Field(Kokkos::Impl::ViewCtorProp<std::string>(label),
                   local_extents,
                   n_ghost)
  {}

  /**
   * @brief Construct a Vector3Field with the given allocation properties, size
   * of the local region, and number of ghost cells
   *
   * @param arg_prop The allocation properties (created with
   * Kokkos::view_alloc)
   * @param local_extents The size of the local region
   * @param n_ghost The number of ghost cells on each side of the axes
   */
  template<
    typename... P,
    typename = std::enable_if_t<!Kokkos::Impl::ViewCtorProp<P...>::has_pointer>>
  Vector3Field(Kokkos::Impl::ViewCtorProp<P...> const& arg_prop,
               Kokkos::Array<int, rank>                local_extents,
               Kokkos::Array<int, rank>                n_ghost = {})
    : x(create_haloview<DataType, Properties...>(
        detail::prepend_ctor_prop_label(arg_prop, "_x"),
        local_extents,
        n_ghost))
    , y(create_haloview<DataType, Properties...>(
        detail::prepend_ctor_prop_label(arg_prop, "_y"),
        local_extents,
        n_ghost))
    , z(create_haloview<DataType, Properties...>(
        detail::prepend_ctor_prop_label(arg_prop, "_z"),
        local_extents,
        n_ghost))
  {}

  // NOLINTBEGIN(*explicit-constructor)
  // Copy constructor for compatible views
  template<typename RT, typename... RP>
  Vector3Field(const Vector3Field<RT, RP...>& other)
    : x{other.x}
    , y{other.y}
    , z{other.z}
  {}
  // NOLINTEND(*explicit-constructor)

  // Copy assignment for compatible views
  template<typename RT, typename... RP>
  Vector3Field& operator=(const Vector3Field<RT, RP...>& other)
  {
    x = other.x;
    y = other.y;
    z = other.z;
    return *this;
  }

  /**
   * @brief Construct a Vector3Field from given x, y, and z arrays
   * @note Arrays are shallow copied, and the resulting Vector3Field will
   *  reference the same data as the input arrays
   */
  template<typename RT1,
           typename RT2,
           typename RT3,
           typename... RP1,
           typename... RP2,
           typename... RP3>
  Vector3Field(HaloView<RT1, RP1...> const& x_component,
               HaloView<RT2, RP2...> const& y_component,
               HaloView<RT3, RP3...> const& z_component)
    : x{x_component}
    , y{y_component}
    , z{z_component}
  {
    if (!check_same_layout_and_offset(x, y, z)) {
      throw std::runtime_error(
        "Incompatible layouts or offsets when constructing Vector3Field from "
        "given components");
    }
  }
};

template<typename DataType, typename... Properties>
int local_extent(const Vector3Field<DataType, Properties...>& vec, int dim)
{
  return local_extent(vec.x, dim);
}

template<typename DataType, typename... Properties>
auto local_extents(const Vector3Field<DataType, Properties...>& vec)
{
  return local_extents(vec.x);
}
} // namespace alps

namespace Kokkos {
template<class ExecSpace, class DT, class... DP>
inline void deep_copy(
  const ExecSpace&                                  exec_space,
  const alps::Vector3Field<DT, DP...>&              dst,
  typename ViewTraits<DT, DP...>::const_value_type& value,
  std::enable_if_t<
    Kokkos::is_execution_space<ExecSpace>::value
    && std::is_void_v<typename ViewTraits<DT, DP...>::specialize>>* = nullptr)
{
  deep_copy(exec_space, dst.x, value);
  deep_copy(exec_space, dst.y, value);
  deep_copy(exec_space, dst.z, value);
}

template<class ExecSpace, class DT, class... DP, class ST, class... SP>
inline void deep_copy(
  const ExecSpace&                     exec_space,
  const alps::Vector3Field<DT, DP...>& dst,
  const alps::Vector3Field<ST, SP...>& src,
  std::enable_if_t<
    (Kokkos::is_execution_space<ExecSpace>::value
     && std::is_void_v<typename ViewTraits<DT, DP...>::specialize>
     && std::is_void_v<typename ViewTraits<ST, SP...>::specialize>)>* = nullptr)
{
  deep_copy(exec_space, dst.x, src.x);
  deep_copy(exec_space, dst.y, src.y);
  deep_copy(exec_space, dst.z, src.z);
}

template<class DT, class... DP>
inline void deep_copy(
  const alps::Vector3Field<DT, DP...>&              dst,
  typename ViewTraits<DT, DP...>::const_value_type& value,
  std::enable_if_t<
    std::is_void_v<typename ViewTraits<DT, DP...>::specialize>>* = nullptr)
{
  deep_copy(dst.x, value);
  deep_copy(dst.y, value);
  deep_copy(dst.z, value);
}

template<class DT, class... DP, class ST, class... SP>
inline void deep_copy(
  const alps::Vector3Field<DT, DP...>& dst,
  const alps::Vector3Field<ST, SP...>& src,
  std::enable_if_t<
    (std::is_void_v<typename ViewTraits<DT, DP...>::specialize>
     && std::is_void_v<typename ViewTraits<ST, SP...>::specialize>)>* = nullptr)
{
  deep_copy(dst.x, src.x);
  deep_copy(dst.y, src.y);
  deep_copy(dst.z, src.z);
}
} // namespace Kokkos
