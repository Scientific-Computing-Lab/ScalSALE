cmake_minimum_required (VERSION 3.5.2)

project(BACKUS)
enable_language (Fortran)

#set(CMAKE_Fortran_COMPILER ifort)
#set(FC "ifort")
#set(CMAKE_Fortran_FLAGS "-fPIC -O0 -cpp -g -fpp -fpe0 -traceback -fp-model source -fimf-arch-consistency=true -fp-stack-check -debug -check bounds -check format -check output_conversion -check pointers  -check uninit -qopenmp")
#set(CMAKE_Fortran_FLAGS "-fPIC -O2 -cpp -g -fpp -fpe0 -traceback -fp-model source -fimf-arch-consistency=true -fp-stack-check -debug -check bounds -check format -check output_conversion -check pointers  -check uninit -qopenmp")
#set(CMAKE_Fortran_FLAGS "-fPIC -cpp -fpp -O2 -qopenmp")
add_definitions($ENV{DEFS})

if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
set(ENV{JSONFORTRAN} "/home/reemh/Downloads/json-fortran_gnu")
set(CMAKE_Fortran_FLAGS "-fPIC -cpp -fopenmp -O0 -ffree-line-length-512")
else()
set(ENV{JSONFORTRAN} "/home/reemh/Downloads/json-fortran")
set(CMAKE_Fortran_FLAGS "-fPIC -cpp -fpp -O2 -g -qopenmp ")
endif()

#set(CMAKE_Fortran_FLAGS "-g -qopenmp -fpp -fPIC -traceback -fpe0")
#set(CMAKE_Fortran_FLAGS "-fPIC -O0 -g -fpp -fpe0 -traceback -qopenmp -cpp -fp-model precise -fp-model precise -fp-model source -fimf-arch-consistency=true -fp-stack-check -debug")
add_definitions($ENV{DEFS})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/archive)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/exec)
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/Scripts)
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/modules)
set(SRC_FILES 
              ${PROJECT_SOURCE_DIR}/Parallel
              ${PROJECT_SOURCE_DIR}/General
              ${PROJECT_SOURCE_DIR}/Input
	     	  ${PROJECT_SOURCE_DIR}/Boundary_conditions
  	     	  ${PROJECT_SOURCE_DIR}/Quantities
	     	  ${PROJECT_SOURCE_DIR}/Boundary_conditions/Vertex
              ${PROJECT_SOURCE_DIR}/Boundary_conditions/Cell
	     	  ${PROJECT_SOURCE_DIR}/Quantities/Cell 
	     	  ${PROJECT_SOURCE_DIR}/Quantities/Vertex
	          ${PROJECT_SOURCE_DIR}/Mesh     	  	     	  
	     	  ${PROJECT_SOURCE_DIR}/Material/Equation_of_state 
	     	  ${PROJECT_SOURCE_DIR}/Material 
	     	  ${PROJECT_SOURCE_DIR}/Rezone_and_Advect
  			  ${PROJECT_SOURCE_DIR}/Time_step
  			  ${PROJECT_SOURCE_DIR}/Diagnostics
#  			  ${PROJECT_SOURCE_DIR}/CR
			  ${PROJECT_SOURCE_DIR}/Main
			  )

foreach(dir ${SRC_FILES})
	add_subdirectory(${dir})
endforeach(dir)
