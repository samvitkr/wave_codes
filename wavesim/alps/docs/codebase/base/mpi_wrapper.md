# MPI wrapper

A thin wrapper around MPI functions is provided in the directory [`mpipp`](#dir_mpipp) or [`namespace mpipp`](#namespace_mpipp). The intent of the wrapper is to provide easier-to-use function signatures of MPI functions, because MPI functions often take a large number of repetitive arguments. Wrapper functions can also determine the correct MPI data type for many intrinsic types, such as `int`, `double`, `float`, `std::complex<double>`, etc.

Example:
```cpp
auto comm = mpipp::COMM_WORLD();
auto rank = comm.rank();

// reduce a local sum to the global sum
// local_sum_h and global_sum are MDViews
mpipp::allreduce(
  local_sum_h.data(), // data() returns the pointer to the underlying data
  local_sum_h.extent_int(0), // number of elements should be int type
  global_sum.data(),
  mpipp::plus<Real>(),
  comm);
```
Original MPI functions can still be used if a wrapper function is not provided. The MPI data type can be obtained by calling `mpipp::get_type<T>()`, which is useful when the type `T` may be changed at compile time. For example, `alps::Real` type may be `double` or `float` depending on the compile-time configuration. Using `mpipp::get_type<Real>()` ensures that the correct MPI data type is used.
```cpp
MPI_Allreduce(
  local_sum_h.data(),
  local_sum_h.extent_int(0),
  mpipp::get_type<Real>(),
  global_sum.data(),
  local_sum_h.extent_int(0),
  mpipp::get_type<Real>(),
  MPI_SUM,
  comm);
```

Refer to [`namespace mpipp`](#namespace_mpipp) for full documentation.