# System-specific setup

### All

It is recommended to set an environment variable `CPM_SOURCE_CACHE` and add to your shell configuration file (e.g. `~/.bashrc` or `~/.zshrc`). For example:
```bash
export CPM_SOURCE_CACHE=${HOME}/.cache/CPM
```

### MSI

Dependencies are installed in the shared project directory at `/home/shenl/shared/opt/spack`. To load the necessary modules, run the following commands:
```bash
module use $MSIPROJECT/shared/opt/spack/share/modules/linux-rocky8-x86_64_v3
module purge
source $MSIPROJECT/shared/opt/spack/environments/alps-deps/loads
```
  > __Note:__ The installed *OpenMPI*, provided via Spack, disables `mpirun` and `mpiexec` commands due to their slow startup times. Use `srun` instead. If you must use `mpirun`, you can call `orterun`, which is identical to `mpirun` in OpenMPI.

Configure:
MSI has multiple GPU types. Use different presets for different GPUs.
```bash
# For A100
cmake -S ${PROJECT_SRC_DIR} -B ${PROJECT_BUILD_DIR} --preset=msi-a100
# For V100
cmake -S ${PROJECT_SRC_DIR} -B ${PROJECT_BUILD_DIR} --preset=msi-v100
# For A40, single precision (double precision on A40 is not very performant)
cmake -S ${PROJECT_SRC_DIR} -B ${PROJECT_BUILD_DIR} --preset=msi-a40 -DALPS_DEFAULT_SINGLE_PRECISION=ON
```

Run:
Use the submit script in the `tools/msi` directory as a template. Keep in mind that the `mpirun` or `mpiexec` in the compiled OpenMPI are disabled because `srun` is the preferred method to invoke the application (see [this page](https://github.com/spack/spack/pull/10340) for more information). If you absolutely needs to use `mpirun`, call `orterun`.

### Nautilus
Load modules:
```bash
module load nvidia/hpc-x/v2.14-gcc/hpcx-mt-ompi
module load scl/gcc-toolset-9
module load cuda/cuda-11.6
```

Install dependencies:
  - HDF5: needs to be built from source

Configure:
```bash
cmake -S ${PROJECT_SRC_DIR} -B ${PROJECT_BUILD_DIR} --preset=nautilus-a100
```

### Narwhal
Load modules:
```bash
module swap PrgEnv-cray PrgEnv-gnu
module load cpe-cuda
module swap cray-mpich cray-mpich-ucx
module swap craype-network-ofi craype-network-ucx
module load craype-accel-nvidia70
module unload cuda/10.1
module swap gcc gcc/11.2.0
module load cray-fftw
module load cray-hdf5-parallel
export CUDA_HOME=/opt/nvidia/hpc_sdk/Linux_x86_64/22.11/cuda/11.8
export LD_LIBRARY_PATH="${CUDA_HOME}/lib64:${LD_LIBRARY_PATH}"
```

Configure:
```bash
cmake -S ${PROJECT_SRC_DIR} -B ${PROJECT_BUILD_DIR} --preset=narwhal-v100
```

Run:
It is expected that the commands in "Load modules" are executed (or included in the job script) before running the application.
```bash
export OMP_NUM_THREADS=16
export OMP_PROC_BIND=close
export OMP_PLACES=cores
mpiexec -n 4 -ppn 2 -d ${OMP_NUM_THREADS} -cpu-bind depth -env MPICH_GPU_SUPPORT_ENABLED=1 -env UCX_TLS=sm,ud,cuda,self ./alps_fst
```
The `mpiexec` command runs the application with 4 MPI ranks, 2 ranks per node (because Narwhal has 2 GPUs per node), and 16 OpenMP threads per rank. The `OMP_PROC_BIND` and `OMP_PLACES` environment variables are set to `close` and `cores`, respectively, to bind OpenMP threads to cores. The `MPICH_GPU_SUPPORT_ENABLED` and `UCX_TLS` environment variables are required to enable GPU support in the MPI library and UCX transport layer, respectively. They are not needed if the application does not use GPUs.
