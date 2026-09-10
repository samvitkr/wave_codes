#pragma once

#include "datatype.h"

#include <mpipp/config.h>

#include <complex>
#include <cstddef>
#include <type_traits>
#include <utility>

//--- forward declarations of Kokkos types----------------------------
namespace Kokkos {
template<class Scalar>
class complex;
} // namespace Kokkos

namespace mpipp {

// Trait to check if a type has raw_handle() returning MPI_Op
template<typename T, typename = void>
struct has_mpi_op_raw_handle : std::false_type
{};

template<typename T>
struct has_mpi_op_raw_handle<
  T,
  std::void_t<decltype(std::declval<const T&>().raw_handle())>>
  : std::is_same<decltype(std::declval<const T&>().raw_handle()), MPI_Op>
{};

template<typename T>
inline constexpr bool has_mpi_op_raw_handle_v = has_mpi_op_raw_handle<T>::value;

// Predefined operations
template<typename DataType>
class max
{
  static_assert(!std::is_same_v<DataType, DataType>,
                "MPI_MAX is not valid for this type");
};

template<typename DataType>
class min
{
  static_assert(!std::is_same_v<DataType, DataType>,
                "MPI_MIN is not valid for this type");
};

template<typename DataType>
class sum
{
  static_assert(!std::is_same_v<DataType, DataType>,
                "MPI_SUM is not valid for this type");
};

template<typename DataType>
using plus = sum<DataType>; // for backward compatibility

template<typename DataType>
class product
{
  static_assert(!std::is_same_v<DataType, DataType>,
                "MPI_PROD is not valid for this type");
};

template<typename DataType>
using multiplies = product<DataType>; // for backward compatibility

template<typename DataType>
class logical_and
{
  static_assert(!std::is_same_v<DataType, DataType>,
                "MPI_LAND is not valid for this type");
};

template<typename DataType>
class logical_or
{
  static_assert(!std::is_same_v<DataType, DataType>,
                "MPI_LOR is not valid for this type");
};

template<typename DataType>
class logical_xor
{
  static_assert(!std::is_same_v<DataType, DataType>,
                "MPI_LXOR is not valid for this type");
};

template<typename DataType>
class bit_and
{
  static_assert(!std::is_same_v<DataType, DataType>,
                "MPI_BAND is not valid for this type");
};

template<typename DataType>
class bit_or
{
  static_assert(!std::is_same_v<DataType, DataType>,
                "MPI_BOR is not valid for this type");
};

template<typename DataType>
class bit_xor
{
  static_assert(!std::is_same_v<DataType, DataType>,
                "MPI_BXOR is not valid for this type");
};

template<typename DataType>
class max_with_location
{
  static_assert(!std::is_same_v<DataType, DataType>,
                "MPI_MAXLOC is not valid for this type");
};

template<typename DataType>
class min_with_location
{
  static_assert(!std::is_same_v<DataType, DataType>,
                "MPI_MINLOC is not valid for this type");
};

template<typename DataType>
class replace;

template<typename DataType>
class no_op;

// -------------------------------------------------------------------

template<typename T, typename Op>
struct is_accumulate_compatible : std::false_type
{};

template<typename T, typename Op>
struct is_get_accumulate_compatible : std::false_type
{};

template<typename T, typename Op>
struct is_fetch_and_op_compatible : std::false_type
{};

template<typename T, typename Op>
inline constexpr bool is_accumulate_compatible_v =
  is_accumulate_compatible<T, Op>::value;

template<typename T, typename Op>
inline constexpr bool is_get_accumulate_compatible_v =
  is_get_accumulate_compatible<T, Op>::value;

template<typename T, typename Op>
inline constexpr bool is_fetch_and_op_compatible_v =
  is_fetch_and_op_compatible<T, Op>::value;

// -------------------------------------------------------------------

#define MPIPP_OP_INTRINSIC(type, op_class, op_value, acc, getacc, fetch) \
  template<>                                                             \
  class op_class<type>                                                   \
  {                                                                      \
   public:                                                               \
    MPI_Op raw_handle() const                                            \
    {                                                                    \
      return op_value;                                                   \
    }                                                                    \
  };                                                                     \
  template<>                                                             \
  struct is_accumulate_compatible<type, op_class<type>>                  \
    : std::bool_constant<acc>                                            \
  {};                                                                    \
  template<>                                                             \
  struct is_get_accumulate_compatible<type, op_class<type>>              \
    : std::bool_constant<getacc>                                         \
  {};                                                                    \
  template<>                                                             \
  struct is_fetch_and_op_compatible<type, op_class<type>>                \
    : std::bool_constant<fetch>                                          \
  {}

// Predefined intrinsic operations as defined in section 5.9.2 of the MPI 3.0
// --------------------------------------------------------------------------
// MPI_MAX
// C integer
MPIPP_OP_INTRINSIC(short, max, MPI_MAX, true, true, true);
MPIPP_OP_INTRINSIC(int, max, MPI_MAX, true, true, true);
MPIPP_OP_INTRINSIC(long, max, MPI_MAX, true, true, true);
MPIPP_OP_INTRINSIC(long long, max, MPI_MAX, true, true, true);
MPIPP_OP_INTRINSIC(signed char, max, MPI_MAX, true, true, true);
MPIPP_OP_INTRINSIC(unsigned char, max, MPI_MAX, true, true, true);
MPIPP_OP_INTRINSIC(unsigned short, max, MPI_MAX, true, true, true);
MPIPP_OP_INTRINSIC(unsigned int, max, MPI_MAX, true, true, true);
MPIPP_OP_INTRINSIC(unsigned long, max, MPI_MAX, true, true, true);
MPIPP_OP_INTRINSIC(unsigned long long, max, MPI_MAX, true, true, true);
// floating point
MPIPP_OP_INTRINSIC(float, max, MPI_MAX, true, true, true);
MPIPP_OP_INTRINSIC(double, max, MPI_MAX, true, true, true);
MPIPP_OP_INTRINSIC(long double, max, MPI_MAX, true, true, true);

// MPI_MIN
// C integer
MPIPP_OP_INTRINSIC(short, min, MPI_MIN, true, true, true);
MPIPP_OP_INTRINSIC(int, min, MPI_MIN, true, true, true);
MPIPP_OP_INTRINSIC(long, min, MPI_MIN, true, true, true);
MPIPP_OP_INTRINSIC(long long, min, MPI_MIN, true, true, true);
MPIPP_OP_INTRINSIC(signed char, min, MPI_MIN, true, true, true);
MPIPP_OP_INTRINSIC(unsigned char, min, MPI_MIN, true, true, true);
MPIPP_OP_INTRINSIC(unsigned short, min, MPI_MIN, true, true, true);
MPIPP_OP_INTRINSIC(unsigned int, min, MPI_MIN, true, true, true);
MPIPP_OP_INTRINSIC(unsigned long, min, MPI_MIN, true, true, true);
MPIPP_OP_INTRINSIC(unsigned long long, min, MPI_MIN, true, true, true);
// floating point
MPIPP_OP_INTRINSIC(float, min, MPI_MIN, true, true, true);
MPIPP_OP_INTRINSIC(double, min, MPI_MIN, true, true, true);
MPIPP_OP_INTRINSIC(long double, min, MPI_MIN, true, true, true);

// MPI_SUM
// C integer
MPIPP_OP_INTRINSIC(short, sum, MPI_SUM, true, true, true);
MPIPP_OP_INTRINSIC(int, sum, MPI_SUM, true, true, true);
MPIPP_OP_INTRINSIC(long, sum, MPI_SUM, true, true, true);
MPIPP_OP_INTRINSIC(long long, sum, MPI_SUM, true, true, true);
MPIPP_OP_INTRINSIC(signed char, sum, MPI_SUM, true, true, true);
MPIPP_OP_INTRINSIC(unsigned char, sum, MPI_SUM, true, true, true);
MPIPP_OP_INTRINSIC(unsigned short, sum, MPI_SUM, true, true, true);
MPIPP_OP_INTRINSIC(unsigned int, sum, MPI_SUM, true, true, true);
MPIPP_OP_INTRINSIC(unsigned long, sum, MPI_SUM, true, true, true);
MPIPP_OP_INTRINSIC(unsigned long long, sum, MPI_SUM, true, true, true);
// floating point
MPIPP_OP_INTRINSIC(float, sum, MPI_SUM, true, true, true);
MPIPP_OP_INTRINSIC(double, sum, MPI_SUM, true, true, true);
MPIPP_OP_INTRINSIC(long double, sum, MPI_SUM, true, true, true);
// complex
MPIPP_OP_INTRINSIC(std::complex<float>, sum, MPI_SUM, true, true, true);
MPIPP_OP_INTRINSIC(std::complex<double>, sum, MPI_SUM, true, true, true);
MPIPP_OP_INTRINSIC(std::complex<long double>, sum, MPI_SUM, true, true, true);
MPIPP_OP_INTRINSIC(Kokkos::complex<float>, sum, MPI_SUM, true, true, true);
MPIPP_OP_INTRINSIC(Kokkos::complex<double>, sum, MPI_SUM, true, true, true);
MPIPP_OP_INTRINSIC(Kokkos::complex<long double>,
                   sum,
                   MPI_SUM,
                   true,
                   true,
                   true);

// MPI_PROD
// C integer
MPIPP_OP_INTRINSIC(short, product, MPI_PROD, true, true, true);
MPIPP_OP_INTRINSIC(int, product, MPI_PROD, true, true, true);
MPIPP_OP_INTRINSIC(long, product, MPI_PROD, true, true, true);
MPIPP_OP_INTRINSIC(long long, product, MPI_PROD, true, true, true);
MPIPP_OP_INTRINSIC(signed char, product, MPI_PROD, true, true, true);
MPIPP_OP_INTRINSIC(unsigned char, product, MPI_PROD, true, true, true);
MPIPP_OP_INTRINSIC(unsigned short, product, MPI_PROD, true, true, true);
MPIPP_OP_INTRINSIC(unsigned int, product, MPI_PROD, true, true, true);
MPIPP_OP_INTRINSIC(unsigned long, product, MPI_PROD, true, true, true);
MPIPP_OP_INTRINSIC(unsigned long long, product, MPI_PROD, true, true, true);
// floating point
MPIPP_OP_INTRINSIC(float, product, MPI_PROD, true, true, true);
MPIPP_OP_INTRINSIC(double, product, MPI_PROD, true, true, true);
MPIPP_OP_INTRINSIC(long double, product, MPI_PROD, true, true, true);
// complex
MPIPP_OP_INTRINSIC(std::complex<float>, product, MPI_PROD, true, true, true);
MPIPP_OP_INTRINSIC(std::complex<double>, product, MPI_PROD, true, true, true);
MPIPP_OP_INTRINSIC(std::complex<long double>,
                   product,
                   MPI_PROD,
                   true,
                   true,
                   true);
MPIPP_OP_INTRINSIC(Kokkos::complex<float>, product, MPI_PROD, true, true, true);
MPIPP_OP_INTRINSIC(Kokkos::complex<double>,
                   product,
                   MPI_PROD,
                   true,
                   true,
                   true);
MPIPP_OP_INTRINSIC(Kokkos::complex<long double>,
                   product,
                   MPI_PROD,
                   true,
                   true,
                   true);

MPIPP_OP_INTRINSIC(bool, logical_and, MPI_LAND, true, true, true);
MPIPP_OP_INTRINSIC(bool, logical_or, MPI_LOR, true, true, true);
MPIPP_OP_INTRINSIC(bool, logical_xor, MPI_LXOR, true, true, true);

MPIPP_OP_INTRINSIC(std::byte, bit_and, MPI_BAND, true, true, true);
MPIPP_OP_INTRINSIC(std::byte, bit_or, MPI_BOR, true, true, true);
MPIPP_OP_INTRINSIC(std::byte, bit_xor, MPI_BXOR, true, true, true);

#define MPIPP_LIST2(A, B) A, B
MPIPP_OP_INTRINSIC(std::pair<MPIPP_LIST2(float, int)>,
                   max_with_location,
                   MPI_MAXLOC,
                   true,
                   true,
                   false);
MPIPP_OP_INTRINSIC(std::pair<MPIPP_LIST2(double, int)>,
                   max_with_location,
                   MPI_MAXLOC,
                   true,
                   true,
                   false);
MPIPP_OP_INTRINSIC(std::pair<MPIPP_LIST2(long, int)>,
                   max_with_location,
                   MPI_MAXLOC,
                   true,
                   true,
                   false);
MPIPP_OP_INTRINSIC(std::pair<MPIPP_LIST2(int, int)>,
                   max_with_location,
                   MPI_MAXLOC,
                   true,
                   true,
                   false);
MPIPP_OP_INTRINSIC(std::pair<MPIPP_LIST2(short, int)>,
                   max_with_location,
                   MPI_MAXLOC,
                   true,
                   true,
                   false);
MPIPP_OP_INTRINSIC(std::pair<MPIPP_LIST2(long double, int)>,
                   max_with_location,
                   MPI_MAXLOC,
                   true,
                   true,
                   false);

MPIPP_OP_INTRINSIC(Kokkos::pair<MPIPP_LIST2(float, int)>,
                   max_with_location,
                   MPI_MAXLOC,
                   true,
                   true,
                   false);
MPIPP_OP_INTRINSIC(Kokkos::pair<MPIPP_LIST2(double, int)>,
                   max_with_location,
                   MPI_MAXLOC,
                   true,
                   true,
                   false);
MPIPP_OP_INTRINSIC(Kokkos::pair<MPIPP_LIST2(long, int)>,
                   max_with_location,
                   MPI_MAXLOC,
                   true,
                   true,
                   false);
MPIPP_OP_INTRINSIC(Kokkos::pair<MPIPP_LIST2(int, int)>,
                   max_with_location,
                   MPI_MAXLOC,
                   true,
                   true,
                   false);
MPIPP_OP_INTRINSIC(Kokkos::pair<MPIPP_LIST2(short, int)>,
                   max_with_location,
                   MPI_MAXLOC,
                   true,
                   true,
                   false);
MPIPP_OP_INTRINSIC(Kokkos::pair<MPIPP_LIST2(long double, int)>,
                   max_with_location,
                   MPI_MAXLOC,
                   true,
                   true,
                   false);

MPIPP_OP_INTRINSIC(std::pair<MPIPP_LIST2(float, int)>,
                   min_with_location,
                   MPI_MINLOC,
                   true,
                   true,
                   false);
MPIPP_OP_INTRINSIC(std::pair<MPIPP_LIST2(double, int)>,
                   min_with_location,
                   MPI_MINLOC,
                   true,
                   true,
                   false);
MPIPP_OP_INTRINSIC(std::pair<MPIPP_LIST2(long, int)>,
                   min_with_location,
                   MPI_MINLOC,
                   true,
                   true,
                   false);
MPIPP_OP_INTRINSIC(std::pair<MPIPP_LIST2(int, int)>,
                   min_with_location,
                   MPI_MINLOC,
                   true,
                   true,
                   false);
MPIPP_OP_INTRINSIC(std::pair<MPIPP_LIST2(short, int)>,
                   min_with_location,
                   MPI_MINLOC,
                   true,
                   true,
                   false);
MPIPP_OP_INTRINSIC(std::pair<MPIPP_LIST2(long double, int)>,
                   min_with_location,
                   MPI_MINLOC,
                   true,
                   true,
                   false);

MPIPP_OP_INTRINSIC(Kokkos::pair<MPIPP_LIST2(float, int)>,
                   min_with_location,
                   MPI_MINLOC,
                   true,
                   true,
                   false);
MPIPP_OP_INTRINSIC(Kokkos::pair<MPIPP_LIST2(double, int)>,
                   min_with_location,
                   MPI_MINLOC,
                   true,
                   true,
                   false);
MPIPP_OP_INTRINSIC(Kokkos::pair<MPIPP_LIST2(long, int)>,
                   min_with_location,
                   MPI_MINLOC,
                   true,
                   true,
                   false);
MPIPP_OP_INTRINSIC(Kokkos::pair<MPIPP_LIST2(int, int)>,
                   min_with_location,
                   MPI_MINLOC,
                   true,
                   true,
                   false);
MPIPP_OP_INTRINSIC(Kokkos::pair<MPIPP_LIST2(short, int)>,
                   min_with_location,
                   MPI_MINLOC,
                   true,
                   true,
                   false);
MPIPP_OP_INTRINSIC(Kokkos::pair<MPIPP_LIST2(long double, int)>,
                   min_with_location,
                   MPI_MINLOC,
                   true,
                   true,
                   false);

#undef MPIPP_LIST2

#undef MPIPP_OP_INTRINSIC

// MPI_REPLACE
template<typename DataType>
class replace
{
  static_assert(mpipp::is_builtin_type_v<DataType>,
                "MPI_REPLACE is only valid for predefined MPI types");

 public:
  MPI_Op raw_handle() const { return MPI_REPLACE; }
};

template<typename DataType>
struct is_accumulate_compatible<DataType, mpipp::replace<DataType>>
  : std::true_type
{};
template<typename DataType>
struct is_get_accumulate_compatible<DataType, mpipp::replace<DataType>>
  : std::true_type
{};
template<typename DataType>
struct is_fetch_and_op_compatible<DataType, mpipp::replace<DataType>>
  : std::true_type
{};

// MPI_NO_OP
template<typename DataType>
class no_op
{
  static_assert(mpipp::is_builtin_type_v<DataType>,
                "MPI_NO_OP is only valid for predefined MPI types");

 public:
  MPI_Op raw_handle() const { return MPI_NO_OP; }
};

template<typename DataType>
struct is_accumulate_compatible<DataType, mpipp::no_op<DataType>>
  : std::false_type
{};
template<typename DataType>
struct is_get_accumulate_compatible<DataType, mpipp::no_op<DataType>>
  : std::true_type
{};
template<typename DataType>
struct is_fetch_and_op_compatible<DataType, mpipp::no_op<DataType>>
  : std::true_type
{};

} // namespace mpipp
