Spectral operations
====================

The spectral operations are implemented in the [`alps::spectral`](#namespace_alps__spectral) namespace. The FFT library is abstracted in the [`alps::fft`](#namespace_alps__fft) namespace.

The {class}`alps::spectral::SpectralGrid` represents a computational grid that can be used for spectral operations. A {class}`alps::spectral::SpectralGrid` instance contains the {class}`alps::PencilPlan` that it relies on, the horizontal base wavenumbers, `pex` and `pey`, and various `alps::spectral::SpectralPlan` that are used for executing the spectral operations on different floating-point types and different backends.

- An alias {type}`alps::Grid` for {class}`alps::spectral::SpectralGrid` is provided for convenience.
```{eval-rst}
.. doxygentypedef:: alps::Grid
```

- Methods are provided to access the properties of the {class}`alps::spectral::SpectralGrid`:
```{eval-rst}
.. doxygenstruct:: alps::spectral::SpectralGrid
   :members:
```

- The list of operations that can be performed on a {class}`alps::spectral::SpectralGrid` are:
```{eval-rst}
.. doxygenfunction:: alps::spectral::ddx

.. doxygenfunction:: alps::spectral::d2dx2

.. doxygenfunction:: alps::spectral::ddy

.. doxygenfunction:: alps::spectral::ddy_and_add

.. doxygenfunction:: alps::spectral::d2dy2

.. doxygenfunction:: alps::spectral::d2dy2_and_add

.. doxygenfunction:: alps::spectral::cutoff_xy

.. doxygenfunction:: alps::spectral::dealias

.. doxygenfunction:: alps::spectral::fft_r2c_xy

.. doxygenfunction:: alps::spectral::fft_c2r_xy
```
Each operation takes a {class}`alps::spectral::SpectralGrid` object and a execution space as arguments. The execution space is used to determine the backend for the operation. For example, if the execution space is `Kokkos::DefaultExecutionSpace`, the operation is executed on the default execution space. If the execution space is `Kokkos::OpenMP`, the operation is executed on the host using OpenMP. The execution space also allows the [asynchronous execution of the operation](base/async.md), e.g. to overlap with MPI communication.

```{warning}
All operations listed above should be considered asynchronous with respect to the host. Therefore, proper synchronization is required to ensure the completion and avoid data races.

Internally, the program maintains temporary buffers to store intermediate data for some spectral operations, e.g. the buffer to store the transposed array in `alps::spectral::ddy` and `alps::spectral::dealias`. The buffers are shared. Therefore different spectral operations should __not__ be called concurrently unless you know internally they do not interfere.
```

Example usage:
```cpp
#include <spectral/spectral.h>

using namespace alps;
auto const grid = Grid(PencilPlan(comm, {nx, ny, nz}), pex, pey);

auto const& partition = grid.partition(); // get the x pencil partition
auto const& partition_y = grid.partition(Pencil::Y); // get the y pencil partition
auto const& comm = grid.comm(); // get the associated MPI communicator

spectral::ddx(out, in, grid, Kokkos::DefaultExecutionSpace());
Kokkos::DefaultExecutionSpace().fence(); // wait for the operation to finish

// Use get_r2c_xy_output_layout to get the layout of the output of fft_r2c_xy
auto d = MDView<double***, default_memory_pool>(
    Kokkos::view_alloc("d", Kokkos::WithoutInitializing),
    grid.get_r2c_xy_output_layout<double, Kokkos::DefaultExecutionSpace>())
spectral::fft_r2c_xy(d, in, grid, Kokkos::DefaultExecutionSpace());
```

## Internal implementation

The spectral operations are delegated to the methods implemented in {class}`alps::SpectralPlan`.

```{eval-rst}
.. doxygenclass:: alps::spectral::SpectralPlan
```
The template parameter `T` is the floating-point type, e.g. `double` or `float`. The template parameter `ExecutionSpace` is the execution space, e.g. `Kokkos::Cuda` or `Kokkos::OpenMP`. The template parameter `Backend` is the FFT backend, e.g. two backends, {class}`alps::fft::CUFFT` and {class}`alps::fft::VKFFT`, can be used with `Kokkos::Cuda`.

By default, `VKFFT` is used for `Kokkos::Cuda` and `FFTW` is used for `Kokkos::OpenMP`. The backend for the device execution space can be specified when constructing the {class}`alps::spectral::SpectralGrid` object as
```cpp
auto const grid = Grid(PencilPlan(comm, {nx, ny, nz}), pex, pey, fft::CUFFT());
```