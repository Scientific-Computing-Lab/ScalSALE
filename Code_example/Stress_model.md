# Stress model addition to ScalSALE
Physical attributes such as stress play an important role in describing the motion of solids in hydrodynamic codes. As a case study, the following workflow describes how the stress module can be easily added to ScalSALE. Thus, expanding the code and demonstrating the idea of bridging the gap between benchmarks and the physical application.

## Prerequisits

This code was tested with:
1. intel/2017 or intel/2018 and OneAPI 2021
2. OpenMPI 1.10.4, 4.0.4, 4.1.3 - any OpenMPI that supports MPI3+ standard.
3. json-fortran https://github.com/jacobwilliams/json-fortran compiled with the same intel.
4. cmake 3.15 or higher

## Folders Documentation

The src folder contains the source code files of ScalSALE, here is a short documentation of its subfolders.

|Folder                    |Documentation                                                                   |
|:---:                     |:---                                                                            |
|**Boundary_conditions**   |Contains the classes that implement the Boundary Conditions                     |
|**CR**                    |Contains the classes that implement the Checkpoint Restart                      |
|**Datafiles**             |Contains all the input datafiles for ScalSALE                                     |
|**General**               |Contains General modules and code files for ScalSALE                              |
|**Input**                 |Contains the classes that parse the input datafile                              |
|**Main**                  |Contains the main code files of ScalSALE                                          |
|**Material**              |Contains the classes that belong to the materials                               |
|**Mesh**                  |Contains the mesh implementation classes                                        |
|**Parallel**              |Contains the Parallelization implementation classes                             |
|**Quantities**            |Contains all the Physical quantities classes in ScalSALE                          |
|**Rezone_and_Advect**     |Contains the implementation of the rezone and advection classes                 |
|**Scripts**               |Contains the Scripts code files                                                 |
|**Time_step**             |Contains the hydrodynamic time step implementation                              |
|**exec**                  |Contains the executable file, created after compilation                         |


