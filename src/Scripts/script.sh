#!/bin/bash
#SBATCH -n 2 -N 1 --exclusive --threads-per-core=1 -p mixedp --error=slurm-%j.err --output=slurm-%j.out
export OMP_NUM__THREADS=1
module load intel/18.0.1.163 openmpi/4.0.4_intel mpi/impi-intel2018 cmake anaconda2
mpi -n 2 /home/reemh/eclipse-workspace_ronw/Backus/src/exec/main
