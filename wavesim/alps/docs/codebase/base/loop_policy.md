# Loop policy

The loop policy is a concept that specifies the loop region as well as how to execute the loop. Kokkos provides a number of loop policies, including `Kokkos::RangePolicy`, `Kokkos::TeamPolicy`, `Kokkos::MDRangePolicy`, etc. See the following documentation from Kokkos for details:
- [Execution policies](https://kokkos.org/kokkos-core-wiki/ProgrammingGuide/ProgrammingModel.html#execution-policies)
- [Hierarchical parallelism](https://kokkos.org/kokkos-core-wiki/ProgrammingGuide/HierarchicalParallelism.html)
- [List of execution policies](https://kokkos.org/kokkos-core-wiki/API/core/Execution-Policies.html)

In this codebase, we provide two aliases for the loop policies:
```{eval-rst}
.. doxygentypedef:: alps::LoopPolicy

.. doxygentypedef:: alps::GridPolicy
```
{type}`alps::LoopPolicy` is an alias of `Kokkos::RangePolicy` or `Kokkos::MDRangePolicy` depending on the rank. {type}`alps::GridPolicy` is an alias of `Kokkos::TeamPolicy`.
The two aliases specify `int` (32-bit integer on most platforms) as the loop index type instead of the default `int64_t`, which may potentially reduce the register pressure and improve performance on GPUs. {type}`alps::LoopPolicy` also fixes the iteration order for `Kokkos::MDRangePolicy` to `Kokkos::Iterate::Left` to be consistent with the column-major ordering of {type}`alps::MDView` and {type}`alps::HaloView` (c.f. [here](multidimensional_array.md)).


```{seealso}
[CUDA compatibility](https://kokkos.org/kokkos-core-wiki/ProgrammingGuide/Interoperability.html#cuda-interoperability) describes how `Kokkos::TeamPolicy` is mapped to CUDA block dimensions.
```