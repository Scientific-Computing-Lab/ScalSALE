#! /bin/bash
#module purge
module load cmake/3.5.2 intel/2017 pFUnit/3.2.9-intel gnuplot/5.0.1 eclipse/2018 
module load anaconda/2.5.0
#module load mlog/1.0.1-gcc-4.9.1 mtime/1.0.0-gcc-4.9.1 hdf5/1.8.9-openmpi-1.6.4-gcc-4.9.1 memtracker/1.0.0-gcc-4.9.1
rm core* 2>/dev/null
rm *.so 2>/dev/null
f2py -c -I.build/modules -m datafile_fortran build/General/CMakeFiles/General.dir/*.o build/Mesh/CMakeFiles/Mesh.dir/*.o build/Main/CMakeFiles/Main.dir/*.o build/Quantities/Cell/CMakeFiles/Q_Cell.dir/*.o build/Quantities/Vertex/CMakeFiles/Q_Vertex.dir/*.o build/Quantities/CMakeFiles/Quantities.dir/*.o build/Material/Equation_of_state/CMakeFiles/EOS.dir/*.o  build/Rezone_and_Advect/CMakeFiles/Rezone_and_Advect.dir/*.o  build/archive/lib*  src/Main/datafile_interface.f90 2>/dev/null
mv *.so src/Main
python src/Main/main.py
