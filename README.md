# ScalSALE: A Scalable Framework for Scientific Code in Modern Fortran and MPI-3
This repository contains the code that originates for the idea presented [https://arxiv.org/abs/1910.06415].

For the datafiles for scaling (lagrange/euler and weak/strong) please refer to https://drive.google.com/file/d/1NXgf8GsVuGG-3ZO87vrbQgpyJ8QRM5Sr/view?usp=sharing
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


## Compile & Install
This project is built via cmake. In order to create the makefiles and the executable you need to run the following commands from the ScalSALE main directory:
```
mkdir build
cd build
cmake ../src
make
```
Or alternatively, enter the directory `src/Scripts` using `cd` command and execute the following bash script:
```
./clean.sh
```
The script creates a folder named build that contains the .o and .a files. In addition, the execution file is located at the directory:
```
src/exec/main
```
Once the makefiles are created, you can simple compile the project with `./make.sh` bash script, located in `src/Scripts` directory.\
\
**Please load the modules before compiling**


## Execution
### Execute Locally
To run the code simply execute the bash script `./run.sh` located in `src/Scripts` directory.\
The script executes the following mpirun command:
```
mpirun -n np ../exec/main
```
Where `np` is the number of mpi processes that are being used.\
In order to run ScalSALE the number of processes in the executation command needs to be similar to the number of processes in the datafile.\
So if you change the number of processes in the datafile **you need to update the number of processes in `run.sh` accordingly**.

#### Execute on a Cluster
To run the code using slurm execute the python script `slurm_run.py` located in `src/Scripts` directory.\
The script checks the number of processes in the datafile and send the relevant job to slurm. The scripts runs the datafile located in `src/Datafiles/datafile.json`\
Please note that you need to change the partition in the job request script:
```
-p mixedp
```
In our case the partition is mixedp, you need to change it to the partition you are using on your cluster.\
In addition, you need to change the module load on the job request to the modules you are loading on your claster (if needed). The relevant part in the job request script is:
```
module load intel/18.0.1.163 openmpi/4.0.4_intel mpi/impi-intel2018 cmake anaconda2
```

## Datafile
The datafile is written as a json file (parsed via the json-fortran module explained above -- it is possible to write a new parser as long as it knows how to read a json file). The main objective of the datafile is to be parsed by the src/Input/datafile_object.f90 module which reads and parses the datafile in order to properly define the execution. 

For the general case, expanding new options, or exploring existing ones can be found in the files datafile_object.f90, replace_words.f90 and defaults.f90 in the folder src/Input. Expanding the datafile can be simply done by following the same coding pattern - adding a default value if needed and allowing words instead of integers and finally parsing them and adding them to the datafile_object.f90 object (that will be passed to src/Main/problem.f90 which builds the classes and physical modules). 

The following segment will describe the datafile in pieces.

1. Materials - Each material is defined by a name, its EoS (currently only ideal gas is allowed) and other properties such as density, initial energy etc:
```json
"Your_Material": {
        "A": 4.0, # atomic number
        "index": 1,  # Index that will appear in the order of materials
        "Z_2": 4.0,  
        "gamma_gas": 1.667, 
        "rho_0": 1.0, # initial density
        "Z": 2.0, 
        "sie_0": 0.0, # if set to 0, the energy will be calculated from the EoS (given a temperature)
        "eos_type": "ideal" # type of EoS, can see in src/Input/replace_words.f90 for more options
    } 
```
2. The layers of the materials - to define the mesh, you first need to define the layers and materials that fill that mesh. It also defines number of cells in each layer. Of course, for the 3d case you should add "number_layers_k"! The following json creates 4 layers, 20x20 cells in each layer - (in case it is an xy mesh, imagine it is a box split into four right in the middle and from left to right then bottom to up it goes material1, material2, material2 and material1.
```json
 "layers_materials": {
    	"number_layers_j": 2,
        "number_cells_j": [
            20,
            20
        ], 
        "number_layers_i": 2, 
        "number_cells_i": [
            20, 
            20
        ], 
        "materials": [
            "Your_material1", 
            "Your_material2",
            "Your_material2", 
            "Your_material1"
        ]
    }, 
```

3. Boundary condition & Mesh type - Every phyiscal problems requires boundary conditions and defenition of which mesh_type. Currently, it supports 3 different meshes: 2d xy, 3d xyz and 3d pyramid. To add a new mesh, simply expand the src/Mesh/mesh_2d or src/Mesh/mesh_3d to build a new way for the X,Y,Z coordinates. the values are integers, which are explained in src/Main/problem. 
```json
"cell_set": {
        "mesh_type": "x_y", 
        "boundary_conditions": [
            2, 
            2, 
            2, 
            2
        ]
    }, 
```

4. Contours - Between each two layers (and the different axis) a contour must be defined for the maximal value, and minimal value of the coordinates. We allow multiple ways to define the contours for example an elipse or a simple straight line. Usually, the j and the k contours are defined as the angles. For a more general case of building a mesh, please edit the code such that the contours_j and contours_k can support such a case. Please note that the symmetry axis is I, therefore usually you will want to define a simple box by giving I the boundaries of the box and the J as a "dummy angle argument" such as follows: this will create a box, specifically 4 boxes with 4 different layers: the first layer is bounded from (0,0) to (0.5,0) , the second layer (0.5,0) to (1,0), the third layer from (0.5,0) to (0.5,0.5) and the last layer from (0.5,0.5) to (1,1). 

```json
"contours": {
        "contours_j": [
            {
                "units": "regular", 
                "theta0": 0.0
            }, 
            {
                "units": "regular", 
                "theta0": 0.5
            }, 
            {
                "units": "regular", 
                "theta0": 1.0
            }
        ], 
        "contours_i": [
            {
                "y1": -1, 
                "x2": 0, 
                "x1": 0, 
                "y2": 2, 
                "contour_type": "line"
            }, 
            {
                "y1": -1, 
                "x2": 0.5, 
                "x1": 0.5, 
                "y2": 2, 
                "contour_type": "line"
            }, 
            {
                "y1": -1, 
                "x2": 1, 
                "x1": 1, 
                "y2": 2, 
                "contour_type": "line"
            }
        ]
    }, 
```

5. Zones - After defining the contour, sometimes the user will want to give each cell a different size, thie is defined in the zone segement. Whereas, constant defines each cell to be in the same size (dr and d_theta are for other types of zones, such a geometry series etc...)
```json
 "zone": {
        "zone_i": [
            {
                "type": "constant", 
                "dr": 0.0
            }, 
            {
                "type": "constant", 
                "dr": 0.0
            }
        ], 
        "zone_j": [
            {
                "d_theta": 0.0, 
                "type": "constant"
            },
            {
                "d_theta": 0.0, 
                "type": "constant"
            }
        ]
    }, 
```

6. Simulation parameters - to change the parameters of the simulation such as final time, and other predefined constants, alter this segment:
```json
   "simulation_parameters": {
        "time_final": 0.3, 
	"cyl": 0,
        "init_temperature": 3.22e-08, 
        "dt": 1e-14, 
        "dt_max": 0.1
    },
```

7. Rezone type - currently only Lagrange and Euler rezoneing are allowed. 0 for lagrange 1 for euler 
```json
"rezone_advect": {
        "rezone_type": 0
    }
```

For more examples please refer to the different datafiles in src/Datafiles. 

## Sedov-Taylor

#### Sedov-Taylor run example
For example this is the initial state of the problem (3D, 15^3 cells)

<p align="center">
<img src="https://github.com/ronwagner311/ExaSALE/blob/main/Images/visit0005.png" width="300" height="300">
</p>


Lagrangian mesh: &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; ALE:

<img align="left" src="https://github.com/ronwagner311/ExaSALE/blob/main/Images/visit0007.png" width="300" height="300">
<img align="right" src="https://github.com/ronwagner311/ExaSALE/blob/main/Images/visit0006.png" width="300" height="300">
