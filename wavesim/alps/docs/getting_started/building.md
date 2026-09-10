# Building

## Requirements

### Compiler versions

Compilers are not tested extensively. Recent versions of the following compilers should work:

- CPU only (OpenMP)
    - gcc >= 8
    - Clang >= 11

When targeting NVIDIA GPUs, one can use either *nvcc* or *Clang* to compile with CUDA support. When *Clang* is used to compile CUDA code directly, a CUDA installation is still needed for CUDA header files and libraries.
The minimum  required versions are:
- GPU + CPU (CUDA)
    - nvcc >= 11.4
    - Clang >= 13 (with CUDA >= 11.0)

  > __Note:__  *nvcc* requires a host compiler, such as *gcc* or *Clang*, to function. 
The compatibility between CUDA and host compiler can be found [here](https://gist.github.com/ax3l/9489132) although there may be exceptions, e.g. certain versions of *gcc* may not work with certain versions of *nvcc* (see [known issues](../known_issues.md)). 

  > __Note:__ *NVHPC* compilers, which may be loaded using `PrgEnv-nvidia` module on Cray machines, are not the same as *nvcc* compiler. They are formerly known as PGI compilers and are not tested with this project. However, the *NVHPC* installation comes with *nvcc* compiler and the CUDA libraries, which may be used to compile the project.

### Dependencies
- FFTW >= 3
  (both double precision and single precision libraries are needed)
- MPI (MPI 3.0 compatible)
- HDF5 with parallel feature

The project also has a few other dependencies, such as [spdlog](https://github.com/gabime/spdlog), which are downloaded and compiled automatically by [CPM.cmake](https://github.com/cpm-cmake/CPM.cmake).

### Build dependency
- CMake >= 3.23

For certain Cuda installations, e.g. from NVHPC, CMake >= 3.22 may be needed to detect cuSparse and cuFFT correctly. Therefore, a newer CMake is preferred.

## Setup and build

Below is a general guide to configure and build the project. More specific instructions on certain systems can be found [here](./linux.md).
In the following steps, `${PROJECT_SRC_DIR}` and `${PROJECT_BUILD_DIR}` refer to the source directory and build directory of the project, respectively. Replace them with actual paths when executing commands. The build directory cannot be the same as the source directory.

1. Preparing dependencies listed above by installing using your package manager or by loading modules.

2. (Optional) It is recommended to set an environment variable `CPM_SOURCE_CACHE`, which is used to cache the dependencies to be downloaded during the configure step (see [__here__](https://github.com/cpm-cmake/CPM.cmake#cpm_source_cache) for more details). For example:
   ```bash
   export CPM_SOURCE_CACHE=${HOME}/.cache/CPM
   ```
   This step can save time when configuring the project multiple times. To make it permanent, add the line above to your shell configuration file, e.g. `~/.bashrc` or `~/.zshrc`.

3. Download sources
    ```bash
    $ git clone https://github.com/path/to/repository ${PROJECT_SRC_DIR}
    ```

4. Configure  
   Configure the project using:
   ```bash
   cmake -S ${PROJECT_SRC_DIR} -B ${PROJECT_BUILD_DIR} <additional options>
   ```
   By default, the project is configured to build for CPU only using OpenMP.
   It is recommended to make the compiler optimize for a specific CPU architecture, add a CMake option `ALPS_ARCH_{cpu}` where `{cpu}` is the targeted CPU arch. One can use `-DALPS_ARCH_NATIVE=ON` to let the compiler choose the best optimization for the current machine, e.g.:
   ```bash
   cmake -S ${PROJECT_SRC_DIR} -B ${PROJECT_BUILD_DIR} -DALPS_ARCH_NATIVE=ON
   ```
   One can also specify the CPU architecture directly. A full list of architectures can be found [__here__](https://kokkos.github.io/kokkos-core-wiki/keywords.html#architecture-keywords). 

   >  __Note:__ All the `ALPS_ARCH_*` options are forwarded as `Kokkos_ARCH_*` to the Kokkos configuration.

   >  __Note:__ The supported architectures depend on the compilers used. Older compilers do not support the compiler flags for newer architectures. For example, the optimization for Zen3 architecture (enabled by `-DALPS_ARCH_ZEN3`) is only supported by gcc >= 11. 

   To configure CUDA support, add `ALPS_ENABLE_CUDA` and `ALPS_ARCH_{gpu}` where `{gpu}` is the GPU architecture to target to the `cmake` line. For example:
   ```
   cmake -S ${PROJECT_SRC_DIR} -B ${PROJECT_BUILD_DIR} -DALPS_ARCH_ZEN2=ON -DALPS_ENABLE_CUDA=ON -DALPS_ARCH_VOLTA70=ON
   ```
   Similarly, the architectures are listed [__here__](https://kokkos.github.io/kokkos-core-wiki/keywords.html#architecture-keywords). For NVIDIA GPUs, the architecture name consists of the GPU family name and the compute capability. For example, `ALPS_ARCH_VOLTA70` targets the V100 GPU. The compute capability can be found [__here__](https://developer.nvidia.com/cuda-gpus).

   A full list of configure options, e.g. for specifying where to look for dependencies and forcing the program to disable GPU-aware MPI, see [CMake options](#cmake-options).

   > __Note:__ If CMake fails to find some dependencies, such as HDF5, consider passing their install locations using the options listed at [CMake options](#cmake-options).


5. Build
   ```bash
   cmake --build ${PROJECT_BUILD_DIR} -j 12 [-t <target>]
   ```
   Adjust the number of threads after `-j` to your need. Optionally, one can specify a target to build after `-t` to build only executables needed and reduce build time, e.g. `-t alps_wind_wave` and `-t alps_channel`.

## CMake options

| Keyword                 | Description                                     | Default |
|-------------------------|-------------------------------------------------|---------|
| `ALPS_ENABLE_CUDA`      | Enable the CUDA backend                         | `OFF`   |
| `ALPS_ARCH_*`           | Enable Kokkos archs, forwarded to Kokkos_ARCH_* | `OFF`   |
| `ALPS_MPI_IS_GPU_AWARE` | Whether MPI is GPU aware                        | `ON`    |
| `ALPS_BUILD_EXAMPLES`   | Whether to build examples                       | `ON`    |
| `ALPS_BUILD_BENCHMARKS` | Whether to build benchmarks                     | `OFF`   |
| `ALPS_BUILD_TESTS`      | Whether to build tests                          | `OFF`   |
| `ALPS_BUILD_HDF5`       | Whether to download and build HDF5 from source  | `OFF`   |

The following options control `find_package` paths

| Keyword            | Description                                                                 |
|--------------------|-----------------------------------------------------------------------------|
| `HDF5_ROOT`        | Top-level installation directory of HDF5                                    |
| `CUDAToolkit_ROOT` | Cuda installation directory (may help find cuFFT and cuSparse libraries)    |

The following options are mainly for code developments:

| Keyword                    | Description                             | Default |
|----------------------------|-----------------------------------------|---------|
| `ALPS_ENABLE_CLANG_TIDY`   | Enable static analysis with Clang-Tidy. | `OFF`   |
| `ALPS_ENABLE_CPPCHECK`     | Enable static analysis with Cppcheck.   | `OFF`   |
| `ALPS_ENABLE_CLANG_FORMAT` | Setup a custom `clangformat` target     | `OFF`   |
