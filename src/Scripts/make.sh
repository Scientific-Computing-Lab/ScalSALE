#! /bin/bash
#module purge
#module load intel/18.0.1.163 openmpi/4.0.4_intel cmake
#module load anaconda2

cd ../../build
#FC=mpiifort #cmake ../src/
#cmake ../src/
#if test "$1" = "GNU"
#then
#  echo "@@@@ GNU @@@@"
#  FC=mpifort cmake ../src/
#else
#  echo "@@@@ INTEL @@@@"
#  FC=mpiifort cmake ../src
#fi
make
cd ../src/Scripts
