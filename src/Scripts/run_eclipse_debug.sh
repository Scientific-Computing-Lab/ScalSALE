#! /bin/bash
module load anaconda/2.5.0
rm src/Input_file/core*
rm src/Input_file/*.so
f2py --debug-capi -c -I./build/modules -m datafile_fortran build/General/CMakeFiles/General.dir/*.o build/Mesh/CMakeFiles/Mesh.dir/*.o build/Problem/CMakeFiles/Problem.dir/*.o build/Quantities/Cell/CMakeFiles/Q_Cell.dir/*.o build/Equation_of_state/CMakeFiles/Equation_of_state.dir/*.o build/Quantities/CMakeFiles/Quantities.dir/*.o build/archive/lib* src/Input_file/datafile_interface.f90
mv *.so src/Input_file
gdb python 
