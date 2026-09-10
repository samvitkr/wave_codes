#pragma once

#include "vector_field.h"
#include "view_types.h"
#include "view_utils.h"

#include <Kokkos_Core.hpp>

#include <string>

namespace alps {

/// @brief 3D field of 3x3 matrices
/**
 * The tensor field is stored in SoA format, and its components are accesses by
 * the member xx, xy, xz, yx, yy, yz, zx, zy, zz.
 */
template<typename DataType, typename... Properties>
class Tensor33Field
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

  view_t xx, xy, xz;
  view_t yx, yy, yz;
  view_t zx, zy, zz;

  /// @brief Construct an empty Tensor33Field
  Tensor33Field() = default;

  /** @brief Construct a Tensor33Field with the given label, size of the local
   * region, and number of ghost cells
   *
   * @param label The label of the field
   * @param local_extents The size of the local region
   * @param n_ghost The number of ghost cells on each side of the axes
   */
  template<
    typename Label,
    typename = std::enable_if_t<Kokkos::Impl::is_view_label<Label>::value>>
  Tensor33Field(Label const&             label,
                Kokkos::Array<int, rank> local_extents,
                Kokkos::Array<int, rank> n_ghost = {})
    : Tensor33Field(Kokkos::Impl::ViewCtorProp<std::string>(label),
                    local_extents,
                    n_ghost)
  {}

  /**
   * @brief Construct a Tensor33Field with the given allocation properties, size
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
  Tensor33Field(Kokkos::Impl::ViewCtorProp<P...> const& arg_prop,
                Kokkos::Array<int, rank>                local_extents,
                Kokkos::Array<int, rank>                n_ghost = {})
    : xx(create_haloview<DataType, Properties...>(
        detail::prepend_ctor_prop_label(arg_prop, "_xx"),
        local_extents,
        n_ghost))
    , xy(create_haloview<DataType, Properties...>(
        detail::prepend_ctor_prop_label(arg_prop, "_xy"),
        local_extents,
        n_ghost))
    , xz(create_haloview<DataType, Properties...>(
        detail::prepend_ctor_prop_label(arg_prop, "_xz"),
        local_extents,
        n_ghost))
    , yx(create_haloview<DataType, Properties...>(
        detail::prepend_ctor_prop_label(arg_prop, "_yx"),
        local_extents,
        n_ghost))
    , yy(create_haloview<DataType, Properties...>(
        detail::prepend_ctor_prop_label(arg_prop, "_yy"),
        local_extents,
        n_ghost))
    , yz(create_haloview<DataType, Properties...>(
        detail::prepend_ctor_prop_label(arg_prop, "_yz"),
        local_extents,
        n_ghost))
    , zx(create_haloview<DataType, Properties...>(
        detail::prepend_ctor_prop_label(arg_prop, "_zx"),
        local_extents,
        n_ghost))
    , zy(create_haloview<DataType, Properties...>(
        detail::prepend_ctor_prop_label(arg_prop, "_zy"),
        local_extents,
        n_ghost))
    , zz(create_haloview<DataType, Properties...>(
        detail::prepend_ctor_prop_label(arg_prop, "_zz"),
        local_extents,
        n_ghost))
  {}

  // NOLINTBEGIN(*explicit-constructor)
  // Copy constructor for compatible views
  template<typename RT, typename... RMS>
  Tensor33Field(const Tensor33Field<RT, RMS...>& other)
    : xx{other.xx}
    , xy{other.xy}
    , xz{other.xz}
    , yx{other.yx}
    , yy{other.yy}
    , yz{other.yz}
    , zx{other.zx}
    , zy{other.zy}
    , zz{other.zz}
  {}
  // NOLINTEND(*explicit-constructor)

  // Copy assignment for compatible views
  template<typename RT, typename... RMS>
  Tensor33Field& operator=(const Tensor33Field<RT, RMS...>& other)
  {
    xx = other.xx;
    xy = other.xy;
    xz = other.xz;
    yx = other.yx;
    yy = other.yy;
    yz = other.yz;
    zx = other.zx;
    zy = other.zy;
    zz = other.zz;
    return *this;
  }

  /// @brief Construct a Tensor33Field from given x, y, and z Vector3Field
  template<typename RT1,
           typename RT2,
           typename RT3,
           typename... RP1,
           typename... RP2,
           typename... RP3>
  Tensor33Field(const Vector3Field<RT1, RP1...>& x,
                const Vector3Field<RT2, RP2...>& y,
                const Vector3Field<RT3, RP3...>& z)
    : xx{x.x}
    , xy{x.y}
    , xz{x.z}
    , yx{y.x}
    , yy{y.y}
    , yz{y.z}
    , zx{z.x}
    , zy{z.y}
    , zz{z.z}
  {
    if (!check_same_layout_and_offset(xx, xy, xz, yx, yy, yz, zx, zy, zz)) {
      throw std::runtime_error(
        "Incompatible layouts or offsets when constructing Tensor33Field from "
        "Vector3Field components");
    }
  }
};

/// @brief 3D field of symmetric 3x3 matrices
/**
 * The tensor field is stored in SoA format, and its components are accesses by
 * the member xx, xy, xz, yx, yy, yz, zx, zy, zz. Symmetric components reference
 * the same data.
 */
template<typename DataType, typename... Properties>
class SymmTensor33Field
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

  view_t xx, xy, xz;
  view_t yx, yy, yz;
  view_t zx, zy, zz;

  /// @brief Construct an empty SymmTensor33Field
  SymmTensor33Field() = default;

  /// @brief Construct a SymmTensor33Field with the given label, size of the
  /// local region, and number of ghost cells
  /**
   * @param label The label of the field
   * @param local_extents The size of the local region
   * @param n_ghost The number of ghost cells on each side of the axes
   */
  template<
    typename Label,
    typename = std::enable_if_t<Kokkos::Impl::is_view_label<Label>::value>>
  SymmTensor33Field(Label const&             label,
                    Kokkos::Array<int, rank> local_extents,
                    Kokkos::Array<int, rank> n_ghost = {})
    : SymmTensor33Field(Kokkos::Impl::ViewCtorProp<std::string>(label),
                        local_extents,
                        n_ghost)
  {}

  /**
   * @brief Construct a SymmTensor33Field with the given allocation
   * properties, size of the local region, and number of ghost cells
   *
   * @param arg_prop The allocation properties (created with
   * Kokkos::view_alloc)
   * @param local_extents The size of the local region
   * @param n_ghost The number of ghost cells on each side of the axes
   */
  template<
    typename... P,
    typename = std::enable_if_t<!Kokkos::Impl::ViewCtorProp<P...>::has_pointer>>
  SymmTensor33Field(Kokkos::Impl::ViewCtorProp<P...> const& arg_prop,
                    Kokkos::Array<int, rank>                local_extents,
                    Kokkos::Array<int, rank>                n_ghost = {})
    : xx(create_haloview<DataType, Properties...>(
        detail::prepend_ctor_prop_label(arg_prop, "_xx"),
        local_extents,
        n_ghost))
    , xy(create_haloview<DataType, Properties...>(
        detail::prepend_ctor_prop_label(arg_prop, "_xy"),
        local_extents,
        n_ghost))
    , xz(create_haloview<DataType, Properties...>(
        detail::prepend_ctor_prop_label(arg_prop, "_xz"),
        local_extents,
        n_ghost))
    , yx(xy)
    , yy(create_haloview<DataType, Properties...>(
        detail::prepend_ctor_prop_label(arg_prop, "_yy"),
        local_extents,
        n_ghost))
    , yz(create_haloview<DataType, Properties...>(
        detail::prepend_ctor_prop_label(arg_prop, "_yz"),
        local_extents,
        n_ghost))
    , zx(xz)
    , zy(yz)
    , zz(create_haloview<DataType, Properties...>(
        detail::prepend_ctor_prop_label(arg_prop, "_zz"),
        local_extents,
        n_ghost))
  {}

  // NOLINTBEGIN(*explicit-constructor)
  // Copy constructor for compatible views
  template<typename RT, typename... RP>
  SymmTensor33Field(const SymmTensor33Field<RT, RP...>& other)
    : xx{other.xx}
    , xy{other.xy}
    , xz{other.xz}
    , yx{xy}
    , yy{other.yy}
    , yz{other.yz}
    , zx{xz}
    , zy{yz}
    , zz{other.zz}
  {}
  // NOLINTEND(*explicit-constructor)

  // Copy assignment for compatible views
  template<typename RT, typename... RP>
  SymmTensor33Field& operator=(const SymmTensor33Field<RT, RP...>& other)
  {
    xx = other.xx;
    xy = other.xy;
    xz = other.xz;
    yx = xy;
    yy = other.yy;
    yz = other.yz;
    zx = xz;
    zy = yz;
    zz = other.zz;
    return *this;
  }
};

template<typename DataType, typename... Properties>
int local_extent(const Tensor33Field<DataType, Properties...>& vec, int dim)
{
  return local_extent(vec.xx, dim);
}

template<typename DataType, typename... Properties>
auto local_extents(const Tensor33Field<DataType, Properties...>& vec)
{
  return local_extents(vec.xx);
}

template<typename DataType, typename... Properties>
int local_extent(const SymmTensor33Field<DataType, Properties...>& vec, int dim)
{
  return local_extent(vec.xx, dim);
}

template<typename DataType, typename... Properties>
auto local_extents(const SymmTensor33Field<DataType, Properties...>& vec)
{
  return local_extents(vec.xx);
}
} // namespace alps

namespace Kokkos {
template<class ExecSpace, class DT, class... DP>
inline void deep_copy(
  const ExecSpace&                                  exec_space,
  const alps::Tensor33Field<DT, DP...>&             dst,
  typename ViewTraits<DT, DP...>::const_value_type& value,
  std::enable_if_t<
    Kokkos::is_execution_space<ExecSpace>::value
    && std::is_void_v<typename ViewTraits<DT, DP...>::specialize>>* = nullptr)
{
  deep_copy(exec_space, dst.xx, value);
  deep_copy(exec_space, dst.xy, value);
  deep_copy(exec_space, dst.xz, value);
  deep_copy(exec_space, dst.yx, value);
  deep_copy(exec_space, dst.yy, value);
  deep_copy(exec_space, dst.yz, value);
  deep_copy(exec_space, dst.zx, value);
  deep_copy(exec_space, dst.zy, value);
  deep_copy(exec_space, dst.zz, value);
}

template<class ExecSpace, class DT, class... DP, class ST, class... SP>
inline void deep_copy(
  const ExecSpace&                      exec_space,
  const alps::Tensor33Field<DT, DP...>& dst,
  const alps::Tensor33Field<ST, SP...>& src,
  std::enable_if_t<
    (Kokkos::is_execution_space<ExecSpace>::value
     && std::is_void_v<typename ViewTraits<DT, DP...>::specialize>
     && std::is_void_v<typename ViewTraits<ST, SP...>::specialize>)>* = nullptr)
{
  deep_copy(exec_space, dst.xx, src.xx);
  deep_copy(exec_space, dst.xy, src.xy);
  deep_copy(exec_space, dst.xz, src.xz);
  deep_copy(exec_space, dst.yx, src.yx);
  deep_copy(exec_space, dst.yy, src.yy);
  deep_copy(exec_space, dst.yz, src.yz);
  deep_copy(exec_space, dst.zx, src.zx);
  deep_copy(exec_space, dst.zy, src.zy);
  deep_copy(exec_space, dst.zz, src.zz);
}

template<class DT, class... DP>
inline void deep_copy(
  const alps::Tensor33Field<DT, DP...>&             dst,
  typename ViewTraits<DT, DP...>::const_value_type& value,
  std::enable_if_t<
    std::is_void_v<typename ViewTraits<DT, DP...>::specialize>>* = nullptr)
{
  deep_copy(dst.xx, value);
  deep_copy(dst.xy, value);
  deep_copy(dst.xz, value);
  deep_copy(dst.yx, value);
  deep_copy(dst.yy, value);
  deep_copy(dst.yz, value);
  deep_copy(dst.zx, value);
  deep_copy(dst.zy, value);
  deep_copy(dst.zz, value);
}

template<class DT, class... DP, class ST, class... SP>
inline void deep_copy(
  const alps::Tensor33Field<DT, DP...>& dst,
  const alps::Tensor33Field<ST, SP...>& src,
  std::enable_if_t<
    (std::is_void_v<typename ViewTraits<DT, DP...>::specialize>
     && std::is_void_v<typename ViewTraits<ST, SP...>::specialize>)>* = nullptr)
{
  deep_copy(dst.xx, src.xx);
  deep_copy(dst.xy, src.xy);
  deep_copy(dst.xz, src.xz);
  deep_copy(dst.yx, src.yx);
  deep_copy(dst.yy, src.yy);
  deep_copy(dst.yz, src.yz);
  deep_copy(dst.zx, src.zx);
  deep_copy(dst.zy, src.zy);
  deep_copy(dst.zz, src.zz);
}

template<class ExecSpace, class DT, class... DP>
inline void deep_copy(
  const ExecSpace&                                  exec_space,
  const alps::SymmTensor33Field<DT, DP...>&         dst,
  typename ViewTraits<DT, DP...>::const_value_type& value,
  std::enable_if_t<
    Kokkos::is_execution_space<ExecSpace>::value
    && std::is_void_v<typename ViewTraits<DT, DP...>::specialize>>* = nullptr)
{
  deep_copy(exec_space, dst.xx, value);
  deep_copy(exec_space, dst.xy, value);
  deep_copy(exec_space, dst.xz, value);
  deep_copy(exec_space, dst.yy, value);
  deep_copy(exec_space, dst.yz, value);
  deep_copy(exec_space, dst.zz, value);
}

template<class ExecSpace, class DT, class... DP, class ST, class... SP>
inline void deep_copy(
  const ExecSpace&                          exec_space,
  const alps::SymmTensor33Field<DT, DP...>& dst,
  const alps::SymmTensor33Field<ST, SP...>& src,
  std::enable_if_t<
    (Kokkos::is_execution_space<ExecSpace>::value
     && std::is_void_v<typename ViewTraits<DT, DP...>::specialize>
     && std::is_void_v<typename ViewTraits<ST, SP...>::specialize>)>* = nullptr)
{
  deep_copy(exec_space, dst.xx, src.xx);
  deep_copy(exec_space, dst.xy, src.xy);
  deep_copy(exec_space, dst.xz, src.xz);
  deep_copy(exec_space, dst.yy, src.yy);
  deep_copy(exec_space, dst.yz, src.yz);
  deep_copy(exec_space, dst.zz, src.zz);
}

template<class DT, class... DP>
inline void deep_copy(
  const alps::SymmTensor33Field<DT, DP...>&         dst,
  typename ViewTraits<DT, DP...>::const_value_type& value,
  std::enable_if_t<
    std::is_void_v<typename ViewTraits<DT, DP...>::specialize>>* = nullptr)
{
  deep_copy(dst.xx, value);
  deep_copy(dst.xy, value);
  deep_copy(dst.xz, value);
  deep_copy(dst.yy, value);
  deep_copy(dst.yz, value);
  deep_copy(dst.zz, value);
}

template<class DT, class... DP, class ST, class... SP>
inline void deep_copy(
  const alps::SymmTensor33Field<DT, DP...>& dst,
  const alps::SymmTensor33Field<ST, SP...>& src,
  std::enable_if_t<
    (std::is_void_v<typename ViewTraits<DT, DP...>::specialize>
     && std::is_void_v<typename ViewTraits<ST, SP...>::specialize>)>* = nullptr)
{
  deep_copy(dst.xx, src.xx);
  deep_copy(dst.xy, src.xy);
  deep_copy(dst.xz, src.xz);
  deep_copy(dst.yy, src.yy);
  deep_copy(dst.yz, src.yz);
  deep_copy(dst.zz, src.zz);
}
} // namespace Kokkos
