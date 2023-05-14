#! /bin/bash

#module load cmake/X.XX.2 eclipse/2018 gcc/9.1.0 intel/2017 mpi/openmpi-1.6.4-gcc-9.1.0 pFUnit/3.2.9-intel-openmpi hdf5/1.8.9-openmpi-1.6.4-gcc-9.1.0 silo/4.8-openmpi-1.6.4-gcc-9.1.0
#module load cmake/X.XX.2 gcc/9.1.0 eclipse/2018 intel/2017 json-fortran/intel-2017 mpi/openmpi-1.10.4-intel-2017 pFUnit/3.2.9-openmpi-1.10.4-intel-2017 hdf5/1.8.9-openmpi-1.10.4-intel-2017 silo/4.8-openmpi-1.10.4-intel-2017 scr/1.2.0-openmpi-1.10.4-intel-2017 
#module purge
#module load intel/18.0.1.163 openmpi/4.0.4_intel cmake 
#module load cmake/X.XX.2 eclipse/2018 anaconda/2.5.0 gcc/9.1.0 intel/2017 pFUnit/3.2.9-intel mpi/openmpi-1.6.4-gcc-9.1.0 hdf5/1.8.9-openmpi-1.6.4-gcc-9.1.0 silo/4.8-openmpi-1.6.4-gcc-9.1.0
#module load anaconda2
#source activate backus-openmpi-1.10.4-intel-2017
#SCR_LIB_FLAGS="-lscrf -L${SCR_PATH}/scr/lib64 -lscr"
#SCR_INCLUDE_FLAGS="-I${SCR_PATH}/scr/include -I/usr/include -I."
rm -rf ../../build
mkdir ../../build
cd ../../build

if test "$1" = "GNU"
then
  echo "@@@@ GNU @@@@"
  FC=mpif90 cmake ../src/
else
  echo "@@@@ INTEL @@@@"
  FC=mpif90 cmake ../src 
fi
#if test "$1" = tests
#then 
#  echo "@@@@ TEST MODE @@@@"
#  FC=mpif90 DEFS="-DTEST" cmake ../src/
#elif test "$1" = debug
#then
#  echo "@@@@ DEBUG MODE @@@@"
#  FC=if90 DEFS="-DDEBUG" cmake ../src/
#else
 # FC=mpiifort cmake ../src/
#fi
make
cd ../src/Scripts
rm core* 2>/dev/null
rm *.so 2>/dev/null
