#!/bin/bash -l
#SBATCH --job-name=re180  # Job name
#SBATCH --nodes=1              # Total number of nodes
#SBATCH --ntasks-per-node=128    # MPI ranks per node
##SBATCH --gres=gpu:v100:2       # Type (must also math partition) and number of GPUs per node
##SBATCH --mem=80g              # This slightly more than the sum of GPU memory per node
#SBATCH -p shen             # Partition
#SBATCH --time=24:00:00         # Run time (d-hh:mm:ss)

SRC_DIR=${SLURM_SUBMIT_DIR}/alps
BUILD_DIR=/tmp/cc/build_cpu
EXE=alps_wind_wave

# Load modules
module use $MSIPROJECT/shared/opt/spack/share/modules/linux-rocky8-x86_64_v3
module purge
source $MSIPROJECT/shared/opt/spack/environments/alps-deps/loads
# export CPM_SOURCE_CACHE=${HOME}/.cache/cpm # Set up the path to dependency cache

# Build the executable
cmake -S ${SRC_DIR} -B ${BUILD_DIR} -DALPS_ARCH_ZEN2=ON  -GNinja -DALPS_DEFAULT_SINGLE_PRECISION=ON
cmake --build ${BUILD_DIR} -t ${EXE} -j 16
rm ${SLURM_SUBMIT_DIR}/run/${EXE}
cp ${BUILD_DIR}/${EXE} ${SLURM_SUBMIT_DIR}/run/${EXE}

# Run
cd ${SLURM_SUBMIT_DIR}/run
srun -n 32 ${SRC_DIR}/scripts/hpcbind --distribute=32 --output-mode=all --openmp-places=cores --lstopo -- ./${EXE} -c param.toml
ret=$? # capture exit status of the previous command

if ((ret)); then # non-zero means failure
  exit "$ret"
fi

# If no error, resubmit
# sbatch submit_cpu.sh
# sbatch submit.sh
