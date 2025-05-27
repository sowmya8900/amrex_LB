#!/bin/bash
#SBATCH -N 6                # Number of nodes
#SBATCH -C cpu             # Use CPU nodes
#SBATCH -q debug           # Use debug queue
#SBATCH -t 00:10:00        # Run time
# Load required modules
module load craype-x86-milan
module load PrgEnv-gnu
module load cray-mpich

# Run the example
srun -n 6 ./main3d.gnu.x86-milan.MPI.OMP.ex test_type=example 