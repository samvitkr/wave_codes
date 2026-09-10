Domain decomposition
======================

## Description of the domain decomposition

The domain decomposition is described by a composition of classes, including {class}`alps::MPIComm3D`, {cpp:class}`alps::BlockPartition`, {cpp:class}`alps::PencilPlan`.

The {cpp:class}`alps::MPIComm3D` class wraps a MPI communicator with a three-dimensional Cartesian topology created from [`MPI_Cart_create`](https://docs.open-mpi.org/en/main/man-openmpi/man3/MPI_Cart_create.3.html). The class also stores some information about the processor grid for quick access. The relation between the linear MPI rank and the Cartesian rank follows the Fortran order, i.e. the axis 0 is the fastest-varying axis, and axis 2 is the slowest-varying axis. This class is trivially copyable and movable.

```{note}
Internally, the Cartesian rank of the MPI communicator uses C-order. The class {cpp:class}`MPIComm3D` reverses the order of the user input when interacting with the MPI API.
```

The {cpp:class}`alps::BlockPartition` class describes a three-dimensional sub-block inside a computational grid. The class also wraps the {cpp:class}`alps::MPIComm3D` processor grid that the computational grid is distributed on. The class is trivially copyable and movable.

The {cpp:class}`alps::PencilPlan` class describes the pencil decomposition of a three-dimensional computational grid. The class includes two {cpp:class}`alps::BlockPartition` objects, {cpp:member}`alps::PencilPlan::x_pencil` and {cpp:member}`alps::PencilPlan::y_pencil`. The underlying {cpp:class}`alps::MPIComm3D` communicator is shared between the two {cpp:class}`alps::BlockPartition` objects while the global extents of the two pencils are transposed in the x-y plane. The number of processors in the first dimension (axis 0) of the {cpp:class}`alps::MPIComm3D` communicator is always 1.

The examples below show the construction of the above classes:
```cpp
// Create a MPI communicator with a 3D Cartesian topology with 1x1x8 processors and periodic in the 0- and 1-axis.
auto const comm = alps::MPIComm3D(mpipp::COMM_WORLD(), {1, 1, 8}, {1, 1, 0});

// Create a pencil decomposition of a 128x128x64 computational grid using the above communicator.
// 64 is distributed on the 8 processors in the 2-axis.
auto const pencil = alps::PencilPlan(comm, {128, 128, 64});
```
A frequently used pattern is to query whether the current MPI process is at the boundary:
```cpp
auto const is_top = pencil.comm().is_last(2); // test whether is the last processor in the 2-axis
auto const is_bottom = pencil.comm().is_first(2); // test whether is the first processor in the 2-axis
```
The following code snippet shows how to get the neighbor MPI rank, which is useful for exchanging ghost cells (see @ref ghost_cell_exchange.h):
```cpp
auto const upper_id = comm.next_proc_on_axis[2];
auto const lower_id = comm.prev_proc_on_axis[2];
```

## Transposition

The transposition of data between the x-pencil and y-pencil distribution can be performed by calling {cpp:func}`alps::transpose_xy` and {cpp:func}`alps::transpose_yx`.

> __Warn__ The transposition functions are always asynchronous, so one must pay attention to the execution order and synchronization.

Example usage:
```cpp
// Below, `in` and `out` are rank-3 alps::MDView
auto const stream = alps::get_next_stream();

alps::transpose_xy(out, in, pencil, stream); // queue the transpose operation on a stream
stream.fence(); // the transpose operation is asynchronous, so wait for it to finish

alps::transpose_yx(in, out, pencil, Kokkos::DefaultExecutionSpace()); // queue the transpose operation on the default stream
Kokkos::DefaultExecutionSpace().fence(); // wait for the transpose to finish
```

### Internal implementation
The function {cpp:func}`alps::transpose_xy` and {cpp:func}`alps::transpose_yx` calls {cpp:func}`alps::transpose::transpose` which looks up a transposition strategy in {cpp:class}`alps::transpose::TransposerPool` (similar to an `std::unordered_map`) using the transposed array dimensions and number of processes. Because some transposition strategy automatically tunes its internal parameters to achieve better performance, the tuned parameters need to be saved. {cpp:class}`alps::transpose::TransposerPool` provides a container to store the tuned strategies. When a transposition strategy is not found in the pool, a new strategy is created and tuned, and then inserted into the pool.

Each transposition strategy is an implementation of the base class {cpp:class}`alps::transpose::TransposerBase` defined in `src/transpose/transposer_base.h`. The following derived classes implement different transposition strategies:
- {cpp:class}`alps::transpose::TransposerSingle`: Single-process transpose (used when number of processes in the transposed dimension is 1)
- {cpp:class}`alps::transpose::TransposerMPIAll2All`: MPI all-to-all based transpose
- {cpp:class}`alps::transpose::TransposerMPIPoint2Point`: MPI point-to-point transpose
- {cpp:class}`alps::transpose::TransposerMPIPoint2PointSHM`: MPI shared-memory optimized transpose

## Halo (ghost) cell exchange

Multiple functions are provided to perform ghost cell exchanges, which all take an {cpp:type}`alps::HaloView` object and an {cpp:class}`alps::BlockPartition` as arguments. The updates are performed in-place. Currently, the exchanges are only supported for the z-axis. The functions provided are:
- Blocking exchange
    - {func}`alps::update_halo_upper_z`
    - {func}`alps::update_halo_lower_z`
    - {func}`alps::update_halo_z`
- Non-blocking exchange
    - {func}`alps::async_update_halo_upper_z`
    - {func}`alps::async_update_halo_lower_z`
    - {func}`alps::async_update_halo_z`

Most of the above functions only work for three-dimensional arrays. The only exception is {func}`alps::update_halo_z`, which works for rank-1 array and is used for exchanging the grid coordinates like `zz` and `zw`. The non-blocking functions use [non-blocking MPI communications](https://enccs.github.io/intermediate-mpi/non-blocking-communication-pt1/), which take or return a {class}`mpipp::irequest_pool` object. See the example below, which is adapted from {func}`alps::solver::div` in `div_curvilinear.cpp`
```cpp
  /* update the upper halos of velocities u, v, w */
  // the first call takes two arguments, returning an irequest_pool object
  auto requests = alps::async_update_halo_lower_z(mesh.grid.x_pencil(), w);
  // the next two calls take three arguments, appending new requests to the irequest_pool object
  alps::async_update_halo_z(requests, mesh.grid.x_pencil(), u);
  alps::async_update_halo_z(requests, mesh.grid.x_pencil(), v);

  /* perform local computation */
  auto stream = alps::get_next_stream();
  /* ... */
  /* the operations must be queued on an async stream rather than the default stream to allow overlap of communication and computation */

  requests.waitall(); // wait for the velocity halo update

  /* perform operations involving the halo cells */
```

### Internal implementation
See [`ghost_cell_exchange.h`](#file_decomp_ghost_cell_exchange.h) for the implementation of the ghost cell exchange functions.