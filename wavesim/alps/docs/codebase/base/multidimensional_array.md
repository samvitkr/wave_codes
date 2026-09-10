# Multi-dimensional array

## `MDView` and `HaloView`

{type}`alps::MDView` and {type}`alps::HaloView` are aliases of [`Kokkos::View`](https://kokkos.org/kokkos-core-wiki/API/core/view/view.html) and [`Kokkos::OffsetView`](https://kokkos.org/kokkos-core-wiki/API/containers/Offset-View.html), respectively, which are used to store multi-dimensional arrays.

```{eval-rst}
.. doxygentypedef:: alps::MDView

.. doxygentypedef:: alps::HaloView

```

Both {type}`alps::MDView` and {type}`alps::HaloView` use the **column-major layout** (or Fortran ordering), i.e. the fastest-varying index is the first axis. These are specified by the template parameter `Kokkos::LayoutLeft` in the above `using` definitions.

The difference between {type}`alps::MDView` and {type}`alps::HaloView` is that the latter allows non-zero starting indices in each dimension. `alps::HaloView` is used to represent a data block with halo regions. Specifically, this codebase follows the convention that the number of halo points is symmetric in each dimension, and the negative indices correspond to the lower halo region. Therefore, the distinction between the local and halo regions implicitly uses the following layout:
```
 -n_halo, ..., -1, 0, 1, ..., local_extent-1, local_extent, ..., local_extent+n_halo-1
|<- halo region ->|<- --- local region -- -->|<------------- halo region ------------>|
```
where `n_halo` is the number of halo points in each dimension, and `local_extent` is the extent of the local region in each dimension.

{type}`alps::MDView` is used when it is not necessary to specify the halo regions.

```{note}
Kokkos::View switches the layout between row-major and column-major depending on the default execution space. We fix the layout to column-major by specifying `Kokkos::LayoutLeft` for {type}`alps::MDView` and {type}`alps::HaloView` to ensure a consistent programming interface.
```

In [`common/container/view_utils.h`](#file_common_container_view_utils.h), some utility functions are provided to query the extents of the local and halo regions of {type}`alps::MDView` and {type}`alps::HaloView` object.
```cpp
#include <common/container/view_types.h>
#include <common/container/view_utils.h>

// a HaloView with 2 halo points in the last dimension
auto const hv = alps::HaloView<double***>("hv", Kokkos::LayoutLeft(64, 32, 17), {0, 0, -1}); // indices start from (0, 0, -1)
auto const mv = alps::MDView<double***>("v", 16, 16, 9);

// get full extents of the view
auto e = alps::extents(hv); // a Kokkos::Array<int, 3> with values {64, 32, 17}
auto e = alps::extents(mv); // a Kokkos::Array<int, 3> with values {16, 16, 9}
auto e0 = alps::extent(hv, 2); // 17
auto e1 = alps::extent(mv, 2); // 9

// get local extents of the view
auto e = alps::local_extents(hv); // a Kokkos::Array<int, 3> with values {64, 32, 15}
auto e = alps::local_extents(mv); // a Kokkos::Array<int, 3> with values {16, 16, 9} (same as extents)
auto e0 = alps::local_extent(hv, 2); // 15

// get the begin and end indices
// useful for defining the loop region
auto v0 = alps::begins(hv); // {0, 0, -1}
auto v1 = alps::ends(hv); // {64, 32, 16}, the end indices are exclusive
auto v1 = alps::ends(mv); // {16, 16, 9}, the end indices are exclusive
auto e0 = alps::begin(hv, 2); // -1
auto e1 = alps::end(hv, 2); // 16

// get the begin and end indices corresponding to the local region
auto v0 = alps::local_begins(hv); // {0, 0, 0}
auto v1 = alps::local_ends(hv); // {64, 32, 15}, the end indices are exclusive
auto v0 = alps::local_begins(mv); // {0, 0, 0}
auto v1 = alps::local_ends(mv); // {16, 16, 9}, the end indices are exclusive
auto e0 = alps::local_begin(hv, 2); // 0
auto e1 = alps::local_end(hv, 2); // 15

// create a subview of the local region of a HaloView 
auto const inner_view = alps:create_inner_view(hv); // a HaloView with extents {64, 32, 15} and begin indices {0, 0, 0}
auto const inner = inner_view.view(); // convert to a MDView
```

## `Vector3Field`, `Tensor33Field` and `SymmTensor33Field`

{class}`alps::Vector3Field`, {class}`alps::Tensor33Field` and {class}`alps::SymmTensor33Field` are SoA (structure of arrays) data structures, which wrap several {type}`alps::HaloView`s to represent the vector, 3x3 matrix and symmetric matrix defined on a 3D domain, respectively. These are often used to store the flow variables in the solver layer.

```cpp
#include <common/container/vector_field.h>
#include <common/container/matrix_field.h>

// construction by specifying the extents
auto const vec_u = alps::Vector3Field<double***>("vec_u", {64, 32, 15}, {0, 0, 1});
auto const u = vec_u.x;
auto const v = vec_u.y;
auto const w = vec_u.z;
```

These classes provide other methods to construct and access, see the documentation and implementation linked below for more details.
- [`common/container/vector_field.h`](#file_common_container_vector_field.h)
- [`common/container/matrix_field.h`](#file_common_container_matrix_field.h)

```{eval-rst}
.. doxygenstruct:: alps::Vector3Field
   :private-members:
   :members:
   :undoc-members:

.. doxygenstruct:: alps::Tensor33Field
   :private-members:
   :members:
   :undoc-members:

.. doxygenstruct:: alps::SymmTensor33Field
   :private-members:
   :members:
   :undoc-members:
```