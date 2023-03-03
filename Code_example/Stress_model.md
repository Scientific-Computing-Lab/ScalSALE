# Stress model addition to ScalSALE
Physical attributes such as stress play an important role in describing the motion of solids in hydrodynamic codes. In non-compressible, visco-elastic, or visco-inelastic fluids, forces between the different layers of the fluid result in stress forces that oppose the fluid flow. The stress is strain and velocity dependent, influencing the fluid's momentum distribution and other hydrodynamic properties. Hence, the stress-strain-velocity relation, which is material-dependent, is crucial for describing the motion of visco-fluids correctly.

 As a case study, the following workflow describes how the stress module can be easily added to ScalSALE. Thus, expanding the code and demonstrating the idea of bridging the gap between benchmarks and the physical application.

Generally, in numerical codes, the physical stress model contains the stress tensor calculation. Consequently, the acceleration is updated accordingly. Therefore, a trivial connection to the hydrodynamics calculation is apparent due to the acceleration updates. However, the stress tensor calculation is independent to the hydrodynamic model. Hence, the stress model in ScalSALE is implemented in a separate class named stress_model (as seen in the code bellow).
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

```fortran
module material_module
  type, extends(material_base_t) :: material_t
    integer, dimension(:) :: stress_module_type 
    ! mat dependent module type (steinberg, perfect, etc)
    contains
      procedure :: Apply_stress_module 
      !calculations of the stress yield and shear
  end type material_t
contains
  subroutine Apply_stress_module(this, shear, yield)
    implicit none
    class(material_t), intent(in out) :: this
    type(shear_modulus_t), pointer, intent(in out) :: shear 
    ! the shear modulus for every material
    type(stress_yield_t), pointer, intent(in out) :: yield 
    ! the yield for every material
    integer :: i, j, k, m
    do m=1, nmats ! loop over all materials
      do k=1, nz ! in the loops we calculate shear and yield for the steinberg model
        do j=1, ny
          do i=1, nx
            ! yield and shear calculations using the Steinberg stress hardening model
            eta = density0(i,j,k)/density(i,j,k)
            shear(m,i,j,k) = shear0+GP_steinberg(m)* &
              pressure(i,j,k)*eta**(1d0/3)+GT_steinberg(m)* &
             (temperature(i,j,k)-init_temperature(i,j,k))
            yield(m,i,j,k) = min(Y0(m)*(1+BETA_steinberg(m)* &
              plastic_strain(i,j,k)) ** N_steinberg(m) &
              ,YIELD_max(m))*(shear(m,i,j,k)/shear0(m))
          end do
        end do
      end do
    end do
    !the MPI synchronization 
    call shear%Exchange_virtual_space_blocking() 
    call yield%Exchange_virtual_space_blocking()
  end subroutine Apply_stress_module
```

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


