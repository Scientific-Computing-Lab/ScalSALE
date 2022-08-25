#! /bin/bash
module purge
module load gcc/4.9.1 intel/2018 cmake/3.5.2 pFUnit/3.2.9-intel gnuplot/5.0.1 eclipse/2018 
module load anaconda/2.5.0
#module load mlog/1.0.1-gcc-4.9.1 mtime/1.0.0-gcc-4.9.1 hdf5/1.8.9-openmpi-1.6.4-gcc-4.9.1 memtracker/1.0.0-gcc-4.9.1
rm core*
rm *.so
f2py -c -I../../build/modules -m datafile_fortran ../../build/General/CMakeFiles/General.dir/*.o ../../build/Mesh/CMakeFiles/Mesh.dir/*.o ../../build/Main/CMakeFiles/Main.dir/*.o ../../build/Time_step/CMakeFiles/Time_step.dir/*.o ../../build/Quantities/Cell/CMakeFiles/Q_Cell.dir/*.o ../../build/Quantities/Vertex/CMakeFiles/Q_Vertex.dir/*.o ../../build/Quantities/CMakeFiles/Quantities.dir/*.o ../../build/Material/Equation_of_state/CMakeFiles/EOS.dir/*.o ../../build/Material/CMakeFiles/Material.dir/*.o ../../build/Boundary_conditions/Cell/CMakeFiles/BC_Cell.dir/*.o ../../build/Boundary_conditions/Vertex/CMakeFiles/BC_Vertex.dir/*.o ../../build/archive/lib*
python ../Main/main.py