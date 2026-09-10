# Memory management

To minimize the overhead of allocating and deallocating temporary arrays on GPUs, we implement a memory pool to manage the memory. The memory pool is implemented in the {class}`alps::memory::DynamicSizePool` and is interfaced to Kokkos through a new templated memory space {class}`alps::PoolSpace`. This {class}`alps::PoolSpace` may be used as a replacement for the default Kokkos memory space.


{class}`alps::PoolSpace` is templated on a Kokkos memory space, which is used to allocate the memory. For example, the space `alps::PoolSpace<Kokkos::CudaSpace>` represents the memory pool built on `Kokkos::CudaSpace`. Additionally, the following aliases are provided to represent the memory pool built on the default device memory space and the host memory space:
```{eval-rst}
.. doxygentypedef:: alps::memory_pool

.. doxygentypedef:: alps::default_memory_pool

.. doxygentypedef:: alps::default_host_memory_pool
```
Example usage of the memory pool space:
```cpp
// vec_u is a Vector3Field<double***> object
auto const& u = vec_u.x;
// allocate a HaloView with the same extents and halo regions as um
HaloView<double***, default_memory_pool> tmp("tmp", u.layout(), u.begins());

// allocate a MDView with the same extents as u on the host
MDView<double***, default_host_memory_pool> u_h("u_h", u.layout());
// create a mirror view of u_h on the device
auto const& u_d = Kokkos::create_mirror_view(default_memory_pool(), u_h);
```

```{seealso}
Kokkos documentation for
- [concept of memory spaces](https://kokkos.org/kokkos-core-wiki/ProgrammingGuide/ProgrammingModel.html#memory-spaces)
- [memory space and Kokkos::View](https://kokkos.org/kokkos-core-wiki/ProgrammingGuide/View.html#managing-data-placement)
- [list of memory spaces](https://kokkos.org/kokkos-core-wiki/API/core/memory_spaces.html)
```