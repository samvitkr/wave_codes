#pragma once

#include <mpipp/config.h>

#include <nonstd/span.hpp>

#include <complex>
#include <cstddef>
#include <type_traits>
#include <utility>

//--- forward declarations of Kokkos types----------------------------
namespace Kokkos {
template<class Scalar>
class complex;

template<class T, class U>
struct pair;
} // namespace Kokkos

namespace Kokkos::Experimental::Impl {
/// @brief templated struct for determining if half_t is an alias to float.
/// @tparam T The type to specialize on.
template<class T>
struct is_float16;

/// @brief templated struct for determining if bhalf_t is an alias to float.
/// @tparam T The type to specialize on.
template<class T>
struct is_bfloat16;
} // namespace Kokkos::Experimental::Impl

namespace mpipp {

template<class T, typename Enable = void>
struct is_builtin_type : std::false_type
{};

template<class T, typename Enable = void>
struct is_derived_type : std::false_type
{};

template<class T, typename Enable = void>
struct is_mpi_type : std::disjunction<is_builtin_type<T>, is_derived_type<T>>
{};

template<class T, typename Enable = void>
struct is_cas_compatible : std::false_type
{};

template<class T>
struct is_builtin_type<T, std::enable_if_t<std::is_enum_v<T>>>
  : is_builtin_type<std::underlying_type_t<T>>
{};

template<class T>
std::enable_if_t<std::is_enum_v<T>, MPI_Datatype>
get_type(std::enable_if_t<std::is_enum_v<T>, int> = 0)
{
  return get_type<std::underlying_type_t<T>>();
}

//--------------------------------------------------------------------

#define MPIPP_PREDEF_TYPE(cxx_type, mpi_type, cas_compatible) \
  template<>                                                  \
  struct is_builtin_type<cxx_type> : std::true_type           \
  {};                                                         \
  template<>                                                  \
  struct is_cas_compatible<cxx_type> : cas_compatible         \
  {};                                                         \
  template<class T>                                           \
  inline MPI_Datatype get_type(                               \
    std::enable_if_t<std::is_same_v<T, cxx_type>, int> = 0)   \
  {                                                           \
    return mpi_type;                                          \
  }

// List of predefined MPI datatype
// https://www.mpi-forum.org/docs/mpi-3.0/mpi30-report.pdf
MPIPP_PREDEF_TYPE(std::nullptr_t, MPI_DATATYPE_NULL, std::false_type)
// Grouping types as defined in section 5.9.2 of the MPI 3.0 standard
MPIPP_PREDEF_TYPE(char, MPI_CHAR, std::false_type)
MPIPP_PREDEF_TYPE(wchar_t, MPI_WCHAR, std::false_type)
// C integer
MPIPP_PREDEF_TYPE(short, MPI_SHORT, std::true_type)
MPIPP_PREDEF_TYPE(int, MPI_INT, std::true_type)
MPIPP_PREDEF_TYPE(long, MPI_LONG, std::true_type)
MPIPP_PREDEF_TYPE(long long, MPI_LONG_LONG, std::true_type)
MPIPP_PREDEF_TYPE(signed char, MPI_SIGNED_CHAR, std::true_type)
MPIPP_PREDEF_TYPE(unsigned char, MPI_UNSIGNED_CHAR, std::true_type)
MPIPP_PREDEF_TYPE(unsigned short, MPI_UNSIGNED_SHORT, std::true_type)
MPIPP_PREDEF_TYPE(unsigned int, MPI_UNSIGNED, std::true_type)
MPIPP_PREDEF_TYPE(unsigned long, MPI_UNSIGNED_LONG, std::true_type)
MPIPP_PREDEF_TYPE(unsigned long long, MPI_UNSIGNED_LONG_LONG, std::true_type)
// floating point
MPIPP_PREDEF_TYPE(float, MPI_FLOAT, std::false_type)
MPIPP_PREDEF_TYPE(double, MPI_DOUBLE, std::false_type)
MPIPP_PREDEF_TYPE(long double, MPI_LONG_DOUBLE, std::false_type)
// logical
MPIPP_PREDEF_TYPE(bool, MPI_CXX_BOOL, std::true_type)
// byte
MPIPP_PREDEF_TYPE(std::byte, MPI_BYTE, std::true_type)
// complex
MPIPP_PREDEF_TYPE(std::complex<float>, MPI_CXX_FLOAT_COMPLEX, std::false_type)
MPIPP_PREDEF_TYPE(std::complex<double>, MPI_CXX_DOUBLE_COMPLEX, std::false_type)
MPIPP_PREDEF_TYPE(std::complex<long double>,
                  MPI_CXX_LONG_DOUBLE_COMPLEX,
                  std::false_type)
MPIPP_PREDEF_TYPE(Kokkos::complex<float>,
                  MPI_CXX_FLOAT_COMPLEX,
                  std::false_type)
MPIPP_PREDEF_TYPE(Kokkos::complex<double>,
                  MPI_CXX_DOUBLE_COMPLEX,
                  std::false_type)
MPIPP_PREDEF_TYPE(Kokkos::complex<long double>,
                  MPI_CXX_LONG_DOUBLE_COMPLEX,
                  std::false_type)

template<class T>
struct is_builtin_type<
  T,
  std::enable_if_t<Kokkos::Experimental::Impl::is_float16<T>::value>>
  : std::true_type
{};
template<class T>
struct is_cas_compatible<
  T,
  std::enable_if_t<Kokkos::Experimental::Impl::is_float16<T>::value>>
  : std::true_type
{};
template<class T>
struct is_builtin_type<
  T,
  std::enable_if_t<Kokkos::Experimental::Impl::is_bfloat16<T>::value>>
  : std::true_type
{};
template<class T>
struct is_cas_compatible<
  T,
  std::enable_if_t<Kokkos::Experimental::Impl::is_bfloat16<T>::value>>
  : std::true_type
{};
template<class T>
inline MPI_Datatype get_type(
  std::enable_if_t<Kokkos::Experimental::Impl::is_float16<T>::value, int> = 0)
{
  return MPI_UINT16_T;
}
template<class T>
inline MPI_Datatype get_type(
  std::enable_if_t<Kokkos::Experimental::Impl::is_bfloat16<T>::value, int> = 0)
{
  return MPI_UINT16_T;
}

// INTERNAL ONLY
#define MPIPP_LIST2(A, B) A, B
MPIPP_PREDEF_TYPE(std::pair<MPIPP_LIST2(float, int)>,
                  MPI_FLOAT_INT,
                  std::false_type)
MPIPP_PREDEF_TYPE(std::pair<MPIPP_LIST2(double, int)>,
                  MPI_DOUBLE_INT,
                  std::false_type)
MPIPP_PREDEF_TYPE(std::pair<MPIPP_LIST2(long, int)>,
                  MPI_LONG_INT,
                  std::false_type)
MPIPP_PREDEF_TYPE(std::pair<MPIPP_LIST2(int, int)>, MPI_2INT, std::false_type)
MPIPP_PREDEF_TYPE(std::pair<MPIPP_LIST2(short, int)>,
                  MPI_SHORT_INT,
                  std::false_type)
MPIPP_PREDEF_TYPE(std::pair<MPIPP_LIST2(long double, int)>,
                  MPI_LONG_DOUBLE_INT,
                  std::false_type)

MPIPP_PREDEF_TYPE(Kokkos::pair<MPIPP_LIST2(float, int)>,
                  MPI_FLOAT_INT,
                  std::false_type)
MPIPP_PREDEF_TYPE(Kokkos::pair<MPIPP_LIST2(double, int)>,
                  MPI_DOUBLE_INT,
                  std::false_type)
MPIPP_PREDEF_TYPE(Kokkos::pair<MPIPP_LIST2(long, int)>,
                  MPI_LONG_INT,
                  std::false_type)
MPIPP_PREDEF_TYPE(Kokkos::pair<MPIPP_LIST2(int, int)>,
                  MPI_2INT,
                  std::false_type)
MPIPP_PREDEF_TYPE(Kokkos::pair<MPIPP_LIST2(short, int)>,
                  MPI_SHORT_INT,
                  std::false_type)
MPIPP_PREDEF_TYPE(Kokkos::pair<MPIPP_LIST2(long double, int)>,
                  MPI_LONG_DOUBLE_INT,
                  std::false_type)
#undef MPIPP_LIST2

#undef MPIPP_PREDEF_TYPE

//--------------------------------------------------------------------

template<class T>
inline constexpr auto is_builtin_type_v = is_builtin_type<T>::value;

template<class T>
inline constexpr auto is_derived_type_v = is_derived_type<T>::value;

template<class T>
inline constexpr auto is_mpi_type_v = is_mpi_type<T>::value;

template<class T>
inline constexpr auto is_cas_compatible_v = is_cas_compatible<T>::value;

} // namespace mpipp
