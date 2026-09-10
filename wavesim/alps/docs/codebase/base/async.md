# Asynchronous execution

Asynchronous execution of kernels on GPUs can be achieved by assigning different device execution spaces to different parallel bodies. In Kokkos, execept for the default execution space `Kokkos::DefaultExecutionSpace`, every new CUDA execution space instance correspond to a CUDA stream. Using streams, one can launch multiple kernels concurrently on the same GPU. Streams are also necessary to overlap computation and communication on GPUs.

```{note}
This section is more relevant to executing multiple kernels concurrently on GPUs. The execution of kernels on GPUs should always be considered asynchronous with respect to the host. See [page 24 of this tutorial](https://github.com/kokkos/kokkos-tutorials/edit/main/LectureSeries/KokkosTutorial_05_SIMDStreamsTasking.pdf) for discussion about the host-device synchronization in Kokkos.

MPI communications from GPUs may also use streams. Because the default stream may block work on other asynchronous streams, __not__ scheduling work on the default stream can help overlap computation and communication on GPUs.
```

```{seealso}
The following tutorial for CUDA provides a good introduction to streams. The concepts are also applicable to Kokkos.
- [CUDA C/C++ Streams and Concurrent from NVIDIA](https://developer.download.nvidia.com/CUDA/training/StreamsAndConcurrencyWebinar.pdf)
- [CUDA training from ENCSS](https://enccs.github.io/cuda/3.02_TaskParallelism/)
```
## Stream management
In our code, to avoid the overhead of creating and destroying streams, we provide a pool of streams, which is implemented in the {class}`alps::StreamPool`. At the initialization of the program, an instance of `alps::StreamPool<Kokkos::DefaultExecutionSpace>` (mapped to {class}`alps::StreamPool\<Kokkos::Cuda>` or {class}`alps::StreamPool\<Kokkos::OpenMP>` depending on the default backend) is created and can be accessed globally from {var}`alps::default_execution_streams`. The stream can be retrieved by calling {func}`alps::get_next_stream()`.

```{note}
{class}`alps::StreamPool\<Kokkos::Cuda>` maintains a circular buffer of CUDA streams. The streams returned by {func}`alps::get_next_stream` are in a round-robin fashion. The calls to {func}`alps::get_next_stream` are thread-safe.
```

```cpp
auto stream1 = alps::get_next_stream();
auto stream2 = alps::get_next_stream();

Kokkos::parallel_for("kernel1", Kokkos::RangePolicy(stream1, 0, 100), KOKKOS_LAMBDA(int i) {
  // kernel 1 is executed on stream1
});
Kokkos::deep_copy(stream2, v_h, v); // copy v to v_h asynchronously on stream2

Kokkos::parallel_for("kernel2", Kokkos::RangePolicy(stream2, 0, 100), KOKKOS_LAMBDA(int i) {
  // kernel 2 is executed on stream2 after the copy is finished
});

stream1.fence(); // wait for all the work on stream1 to finish
stream2.fence(); // wait for all the work on stream2 to finish
```