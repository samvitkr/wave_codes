This file provides guidance to AI coding agents like Claude Code (claude.ai/code), Cursor AI, Codex, Gemini CLI, GitHub Copilot, and other AI coding assistants when working with code in this repository.

# ALPS Agent Development Guidelines

ALPS is a boundary-fitted-grid Navier-Stokes solver for turbulent flows using Kokkos for portable parallelism and MPI for distributed execution.

## Build System

Local CMake build configurations and toolchain management with Pixi are primarily designed for development, refactoring, and local execution phases. When executing code on supercomputers or HPC clusters, agents should prioritize system-specific configurations and utilize vendor-provided packages for critical dependencies, such as MPI, to ensure compatibility with specialized hardware interconnects and runtime environments.

### Environment Setup
- Uses pixi for toolchains and MPI: `pixi install -a --frozen`
- Optional: Set `CPM_SOURCE_CACHE` environment variable to cache dependencies (e.g., `export CPM_SOURCE_CACHE=${HOME}/.cache/CPM`)

### Build environment provided by Pixi
- Primary toolchains provided by pixi environments: GCC/Clang for CPU, CUDA/Clang or CUDA/NVCC for GPU
- MPI required for parallel execution
- HDF5
- FFTW

### Build Commands
Use `CMakeUserPresets.json` for local builds. If missing, create it from the template in `scripts/CMakeUserPresets.json`.

All build commands should run in Pixi environment by prefixing commands with `pixi run` or `pixi run -e cuda` to ensure a consistent toolchain.

#### CPU Presets (GCC)
```bash
# Configure build
pixi run cmake --preset local-cpu

# Build the project
pixi run cmake --build build-local-cpu

# Run tests
pixi run ctest --test-dir build-local-cpu --output-on-failure

# Run single test
pixi run ctest --test-dir build-local-cpu -R "test_name" --output-on-failure
```

#### CUDA Presets (Clang)

##### double precision
```bash
# Configure build
pixi run -e cuda cmake --preset local-cuda

# Build the project
pixi run -e cuda cmake --build build-local-cuda

# Run tests
pixi run -e cuda ctest --test-dir build-local-cuda --output-on-failure

# Run single test
pixi run -e cuda ctest --test-dir build-local-cuda -R "test_name" --output-on-failure
```

##### single precision, preferred during development
```bash
# Configure build
pixi run -e cuda cmake --preset local-cuda-single

# Build the project
pixi run -e cuda cmake --build build-local-cuda-single

# Run tests
pixi run -e cuda ctest --test-dir build-local-cuda-single --output-on-failure

# Run single test
pixi run -e cuda ctest --test-dir build-local-cuda-single -R "test_name" --output-on-failure
```

#### Other Presets
```bash
# Build with NVCC
pixi run -e cuda cmake --preset local-cuda-nvcc-gcc
pixi run -e cuda cmake --build build-local-cuda-nvcc-gcc
```

#### Detecting CUDA Architecture
**CRITICAL FOR AGENTS**: When running ANY CUDA preset
1. **FIRST** detect the GPU with `nvidia-smi -L`
2. **THEN** add the appropriate `-DALPS_ARCH_<NAME>=ON` flag to cmake configure command to override the default arch

```bash
# Step 1: Detect GPU
nvidia-smi -L

# Step 2: Use the detected architecture
cmake --preset local-cuda-single -DALPS_ARCH_ADA89=ON  # Ada (RTX 40xx, RTX 1000 Ada)
cmake --preset local-cuda-single -DALPS_ARCH_AMPERE86=ON  # Ampere (RTX 30xx, A100)
cmake --preset local-cuda-single -DALPS_ARCH_HOPPER90=ON  # Hopper (H100)
```

## High-Level Architecture

The codebase follows a layered architecture:

1. **Parallelization Layer**: Kokkos (CUDA/OpenMP), MPI (via `mpipp` wrapper)
2. **Base Layer** (`src/common/`):
   - `async/`: `StreamPool` for asynchronous execution, `Event` for synchronization
   - `container/`: `MDView` (alias for `Kokkos::View`), `HaloView` (alias for `Kokkos::Experimental::OffsetView`), `Vector3Field` (SoA format)
   - `kokkos_abstraction/`: `PoolSpace` for optimized memory allocation (`default_memory_pool`)
3. **Transpose Layer** (`src/transpose/`):
   - `TransposerPool`: Manages reusable transposition plans
   - `TransposerBase`: Base class for transpose strategies
   - MPI transpose strategies (All2All, Point2Point, Point2PointSHM)
   - Operations: `transpose`
4. **Domain Decomposition Layer** (`src/decomp/`):
   - `PencilPlan`: Manages 2D pencil decomposition (X-pencil, Y-pencil)
   - `GhostCellExchange`: Halo updates (`async_update_halo_z`)
5. **Spectral Operations Layer** (`src/spectral/`):
   - `SpectralGrid`: Manages spectral plans and FFT backends (VkFFT, CUFFT, FFTW)
   - Operations: `ddx`, `ddy`, `d2dx2`, `d2dy2`, `dealias`, `fft_r2c_xy`, `fft_c2r_xy`
6. **Solver Layer** (`src/solvers/`):
   - `Mesh`: Base class for 3D grids; `CurvilinearMesh` (e.g., `BottomWaveMesh`, `TopWaveMesh`) for boundary-fitted grids
   - `FlowField`: Stores velocity (`u`), pressure (`pp`), and boundary conditions
   - `ChannelFlowSolver`: Solvers for channel flows (`AB2`, `AB2CN`)
   - `FlowOverWaveSolver` / `FreeSurfaceSolver`: Specialized solvers for curvilinear domains
   - `PressureCurvilinearEqn`: Iterative Poisson solver using tridiagonal systems
7. **Application Layer** (`src/apps/`): Driver code and CLI utilities

## Codebase Organization

- `src/`: Core source code
  - `apps/`: Main application drivers and CLI utilities
  - `common/`: Base infrastructure
    - `async/`: Asynchronous execution (`StreamPool`) and synchronization (`Event`)
    - `base/`: Logging, macros, and basic utilities
    - `container/`: Multi-dimensional views (`MDView`, `HaloView`) and vector fields (`Vector3Field`)
    - `device/`: Device traits and management
    - `kokkos_abstraction/`: Kokkos-specific wrappers, execution policies, and `PoolSpace`
    - `memory/`: Memory allocation resources and dynamic pools
    - `program_options/`: CLI and configuration file parsing
    - `runtime/`: Runtime management and asynchronous utilities
  - `decomp/`: Domain decomposition and halo exchange logic
  - `transpose/`: Transposition between pencils and MPI communication strategies
  - `io/`: HDF5 and file I/O operations
  - `solvers/`: Physical solvers and mesh definitions
    - `field/`: Flow field definitions and boundary condition types
    - `hos/`: Higher-Order Spectral (HOS) wave solvers
    - `mesh/`: 3D grid definitions (Cartesian and Curvilinear)
    - `ns/`: Navier-Stokes solvers
      - `channel/`: Solvers for canonical channel flows (AB2, AB2CN)
      - `curvilinear_common/`: Shared kernels for curvilinear grid solvers (pressure, viscous, etc.)
      - `flow_over_wave/`: Solver for flow over fixed wavy boundaries
      - `free_surface/`: Solver for free-surface flows
    - `operators/`: Discrete differential operators (grad, div) for various grids
    - `poisson/`: Tridiagonal and Poisson equation solvers
    - `source_terms/`: Physical source terms (Coriolis, Boussinesq, Rayleigh damping)
    - `turbulence_model/`: Subgrid-scale (SGS) models (Smagorinsky, AMD)
  - `spectral/`: FFT-based spectral operations and derivatives
- `test/`: Unit and integration tests using Catch2
- `examples/`: Example configurations and simulation setups
- `benchmark/`: Performance benchmarks for FFT and transpositions
- `scripts/`: Utility scripts for building, profiling, and submission
- `docs/`: Project documentation (Sphinx/Doxygen)

## Code Style Guidelines

### Coding standard
- C++17 standard
- Style: LLVM-based (see `.clang-format`)
- Indentation: 2 spaces
- Prefer east const (`Type const name`, `Type const& name`, `Type const* name`)

### Headers and Include Order
```cpp
#pragma once
// System headers
#include <memory>
// Third-party headers
#include <Kokkos_Core.hpp>
// Project headers (relative to ./src)
#include <common/base/macros.h>
```

### Naming Conventions
- **Classes**: PascalCase (e.g., `FlowField`)
- **Functions/Variables**: snake_case (e.g., `calc_uhat`, `local_mesh`)
- **Constants**: UPPER_SNAKE_CASE

### Type System
- **Real type**: Use `alps::Real` (configurable precision)
- **Views**:
  - `MDView<Real***>`: Multi-dimensional arrays
  - `HaloView<Real***>`: Fields with halo regions (OffsetView)
  - `Vector3Field<Real***>`: SoA 3D vector field consisting of `HaloView`
- **Layout**: `MDView` and `HaloView` use the layout `Kokkos::LayoutLeft` (column-major).

## Development Patterns

### Asynchronous Execution
- Use `alps::get_next_stream()` to retrieve streams from the `StreamPool`
- All spectral operations and transpositions are asynchronous
- Use `stream.fence()` or `Kokkos::DefaultExecutionSpace().fence()` for synchronization
- Use `Event` and `enqueue`/`wait_for` for cross-stream synchronization

### Memory Management
- Use `alps::default_memory_pool` for frequent allocations to minimize GPU overhead
- Prefer RAII and `std::unique_ptr` for object ownership

### GPU Programming
- Use `KOKKOS_LAMBDA` for device lambdas
- Use `alps::LoopPolicy` (Range/MDRange) and `alps::GridPolicy` (TeamPolicy) for parallel kernels

### Domain Decomposition
- Data is typically decomposed into pencils
- Transposition is required when switching between X-derivative (X-pencil) and Y-derivative (Y-pencil)

### Spectral Operations
- FFT backends: VKFFT (default for CUDA), CUFFT, FFTW
- Operations are asynchronous and may share temporary buffers
- Dealiasing (2/3 rule) is typically required after non-linear terms

## Testing
- Framework: Catch2
- Test files: `test/**/*_test.cpp`
- Use `TEMPLATE_TEST_CASE_SIG` for testing across different types and execution spaces

## Development Workflow

### After Making Changes
1. Format code: `pixi run clang-format -i file.cpp`
2. Format CMake: `pixi run cmake-fmt`
