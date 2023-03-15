#!/bin/bash
#SBATCH -n 27 -N 1 --exclusive --threads-per-core=1 -p backus  --error=slurm-%j.err --output=slurm-%j.out
export OMP_NUM_THREADS=1
ulimit -s unlimited
module load intel/18.0.1.163 openmpi/4.1.3-intel
mpirun -n 27 ../exec/main
