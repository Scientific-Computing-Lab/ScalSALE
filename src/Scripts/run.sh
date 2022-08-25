#!/bin/bash
#module purge
#module load cmake/X.XX.2 eclipse/2018 anaconda/2.5.0 gcc/9.1.0 intel/2017 mpi/openmpi-1.6.4-gcc-9.1.0 pFUnit/3.2.9-intel-openmpi hdf5/1.8.9-openmpi-1.6.4-gcc-9.1.0 silo/4.8-openmpi-1.6.4-gcc-9.1.0
#module load cmake/X.XX.2 eclipse/2018 gcc/9.1.0 intel/2017 json-fortran/intel-2017 mpi/openmpi-1.10.4-intel-2017 pFUnit/3.2.9-openmpi-1.10.4-intel-2017 hdf5/1.8.9-openmpi-1.10.4-intel-2017 silo/4.8-openmpi-1.10.4-intel-2017 scr/1.2.0-openmpi-1.10.4-intel-2017
#module load anaconda/2.5.0
#source activate backus-openmpi-1.10.4-intel-2017
#module load intel/18.0.1.163 openmpi/4.0.4_intel cmake anaconda2

#SCR_LIB_FLAGS="-lscrf -L${SCR_PATH}/scr/lib64 -lscr"
#SCR_INCLUDE_FLAGS="-I${SCR_PATH}/scr/include -I/usr/include -I."
rm -rf Silo_Diagnostics/*
rm *.so 2>/dev/null

export SCR_CONF_FILE=`pwd`/../CR/scr_conf.conf
export SCR_RUNS=4
if test "$1" = "cr"
then
	scr_mpirun -n 1 -mca btl self,sm,openib python ../Main/main.py
else
#	mpirun -n 1 -mca btl self,sm,openib ../exec/main
    #mpirun -n 1 -mca btl self,sm,openib ../exec/main
    #mpirun -n 8 ../exec/main
    mpirun -n $1 ../exec/main
fi
