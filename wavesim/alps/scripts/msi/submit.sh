#!/bin/bash -l
#SBATCH --job-name=examplejob  # Job name
#SBATCH --nodes=1              # Total number of nodes 
#SBATCH --cpus-per-task=8      # Cores per task (>1 if multi-threaded tasks)
#SBATCH --ntasks-per-node=2    # MPI ranks per node
#SBATCH --gres=gpu:a100:2      # Type (must also math partition) and number of GPUs per node
#SBATCH --mem=100g             # This slightly more than the sum of GPU memory per node
#SBATCH -p a100-4              # Partition
#SBATCH --time=8:00:00         # Run time (d-hh:mm:ss)


# Load modules
module use /home/shenl/shared/opt/spack/share/modules/linux-rocky8-x86_64_v3
module purge
source /home/shenl/shared/opt/spack/environments/alps-deps/loads

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
## TEST! Sometimes binding may hurt performance if not exclusively using a node.
# export OMP_PLACES=cores
# export OMP_PROC_BIND=close

cd ${SLURM_SUBMIT_DIR}
srun --mpi=pmi2 -n $SLURM_NTASKS ./alps_channel