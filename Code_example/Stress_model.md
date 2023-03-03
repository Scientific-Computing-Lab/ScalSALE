# Stress model addition to ScalSALE
Physical attributes such as stress play an important role in describing the motion of solids in hydrodynamic codes. In non-compressible, visco-elastic, or visco-inelastic fluids, forces between the different layers of the fluid result in stress forces that oppose the fluid flow. The stress is strain and velocity dependent, influencing the fluid's momentum distribution and other hydrodynamic properties. Hence, the stress-strain-velocity relation, which is material-dependent, is crucial for describing the motion of visco-fluids correctly.

 As a case study, the following workflow describes how the stress module can be easily added to ScalSALE. Thus, expanding the code and demonstrating the idea of bridging the gap between benchmarks and the physical application.

In order to add new physical model to ScalSALE a class of the new physical model class need to be created. This class will contain all the relevant quantities and calculations for the new model. The stress model class can be seen below:
```fortran
type :: stress_model
  type(shear_modulus_t) , pointer :: shear
  type(stress_yield_t)  , pointer :: yield
  type(stress_tensor_t) , pointer :: stress_tensor
  contains
    procedure :: Calculate_stress
    procedure :: Update_acceleration
end type stress_model
```

*This class described the new physical stress module. It contains two main functions, a function that calculates the stress in the current time step and a function that update the acceleration values using the updated stress tensor.*
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


