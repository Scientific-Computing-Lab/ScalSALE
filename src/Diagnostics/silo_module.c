/* ############################################################################# */
/* Headers and prototypes */

/* C headers */

#include <silo.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <sys/syscall.h>
#include <errno.h>
#include <sys/types.h>
#include <time.h>
#include <errno.h>


/* Functions prototypes */
int dimensions;

int c_silo_set_file_num(int file_num);
int c_silo_set_files_path(char *path);

int c_silo_init_2d(double T, int NCYC, int file_sn, int MEMAD, 
		int NXP, int NYP, 
		double *XS, double *YS,
		char *SILO_DIAGNOSTIC_LIBRARY);
int c_silo_init_3d(double T, int NCYC, int file_sn, int MEMAD, 
					int NXP, int NYP, int NZP,
					double *XS, double *YS, double *ZS,
					char *SILO_DIAGNOSTIC_LIBRARY);

int silo_initializing_values(double T, int NCYC, int file_sn, int MEMAD, 
		char *SILO_DIAGNOSTIC_LIBRARY); //// Modify by c_silo_init
int silo_set_paths();
//int silo_PMPIO_connectivity();
int silo_mesh_setting_2d(double *XS, double *YS, int NXP, int NYP);
int silo_mesh_setting_3d(double *XS, double *YS, double *ZS, int NXP, int NYP, int NZP);

int c_silo_write_node_data_int(char *var_name, int *values);
int c_silo_write_node_data_float(char *var_name, float *values);
int c_silo_write_node_data_double(char *var_name, double *values);

int c_silo_write_zone_data_int(char *var_name, int *values);
int c_silo_write_zone_data_float(char *var_name, float *values);
int c_silo_write_zone_data_double(char *var_name, double *values);

int c_silo_write_node_data_double_vector(char *var_name, double *values, double *values2, double *values3);
//int c_silo_write_zone_data_double_vector(char *var_name, double *values);

int c_silo_finalize();
int silo_write_multi_mesh_labels();
int silo_write_multi_vars_labels();

int write_to_bugfile(char* error_message);


/* ############################################################################# */
/* Global data struct */

typedef struct _c_silo_info {

	// Silo mesh
	int memad;                      // Dimentions of the simulation [1,2,3]
	int nv;                         // Number of vertexes [1..]
	int file_sn;                    // File serial number [1, 2, ..]
	int *dims_sizes;                // dimentionas sizes [235423, 4234, 423423]

	// Cycle and time
	int cycle;                      // Cycle of the run [1, 345, 678..]
	double time;                    // Time of the run [1e-06, ..]
	DBoptlist *silo_opt_list;       // Contains the time and the cycle for the plot

	// Coordinates
	char **coordnames;              // Coordinates names
	double **coords;                // Coordinates array

	// Formats
	char *silo_file_ext;            // Files name extantion [hdf5, pdb]
	int silo_driver;                // Files driver [hdf5, pdb]

	// Silo files
	char *current_dir_name;         // Current Silo directory name 
	char *directory_path_name;      // Current Silo directory path
	char *silo_file_name;           // Silo file name ["silo_diagnostics_for_leeor3d"]
	char *silo_domain_name;         // Silo domain (=rank) name ["domain000", "domain001",..]
	DBfile *silo_file;              // Pointer to the file [*]

	// Values variables
	int silo_var_count;             // Amount of data vectors which the Silo file will get [1, 435, 554..]
	char **silo_var_names;          // Names of the data vectors ["pressure", ..]

	// Memory sum
	double rank_total_memory_write; // Total amount of data (in bytes) which been written to the Silo files [34543543, 43563454..]

	// Log files
	FILE *rank_writing_data;        // Log file for the ranks [.txt]
	FILE *ranks_meta_data;          // Log file for the data [.txt]

	// Timers
	time_t t_start, t_end;                  // start-end time
	clock_t clock_start, clock_end;         // start-end clock time

	// Error
	char *error_message; // Error message 

} c_silo_info_t;

// Global struct object of type c_silo_info_t
c_silo_info_t _gsi;


/* ############################################################################# */
/* Seperated initialization methods (for further use) */

// Setting the files path for silo diagnostics
// (in seperate function for need to set directly from Leeor Fortran code)
int c_silo_set_files_path(char *path) {
	_gsi.current_dir_name = strdup(path);
	return 0;
}

/* ############################################################################# */
/* Silo initialization */

/*-----------------------------------------------------------------------------
 * Values c_silo_init_2d() gets from Leeor2d 2d.f90:
 *
 * 1)  double T                      : real time of the current plot (0, 1e-6, 2.0006e-6, ...)
 * 2)  int NCYC                      : real cycle of the current plot (0, 456, 674, 2342, ...)
 * 3)  int file_s                    : Serial number of the current plot (0,1,2.., ...)  
 * 4)  int MEMAD                     : real dimention of the current run of the leeor2d (1,2,3)
 * 5)  int NXP                       : Dimention across axis X
 * 6)  int NYP                       : Dimention across axis Y
 * 7)  double XS*                    : Coordinates array with X points (size = NV) 
 * 8)  double YS*                    : Coordinates array with Y points (size = NV)
 * 9) char *SILO_DIAGNOSTIC_LIBRARY : Name of the Silo diagnostic library ("HDF5", "PDB")
 *
 *-----------------------------------------------------------------------------
 */

int c_silo_init_2d(double T, int NCYC, int file_sn, int MEMAD,
		int NXP, int NYP,
		double *XS, double *YS,   
		char *SILO_DIAGNOSTIC_LIBRARY)
{

	if (silo_initializing_values(T, NCYC, file_sn, MEMAD, SILO_DIAGNOSTIC_LIBRARY) != 0) return 1;
	if (silo_set_paths() != 0) return 1;
	if (silo_mesh_setting_2d(XS, YS, NXP, NYP) != 0) return 1;                       
        dimensions = 2;

	return 0;

}

/*-----------------------------------------------------------------------------
 * Values c_silo_init_3d() gets from Leeor3d 3d.f90:
 *
 * 1)  double T                      : real time of the current plot (0, 1e-6, 2.0006e-6, ...)
 * 2)  int NCYC                      : real cycle of the current plot (0, 456, 674, 2342, ...)
 * 3)  int file_s                    : Serial number of the current plot (0,1,2.., ...)  
 * 4)  int MEMAD                     : real dimention of the current run of the leeor3d (1,2,3)
 * 5)  int NXP                       : Dimention across axis X
 * 6)  int NYP                       : Dimention across axis Y
 * 7)  int NZP                       : Dimention across axis Z
 * 8)  double XS*                    : Coordinates array with X points (size = NV) 
 * 9)  double YS*                    : Coordinates array with Y points (size = NV)
 * 10) double YS*                    : Coordinates array with Z points (size = NV)
 * 11) char *SILO_DIAGNOSTIC_LIBRARY : Name of the Silo diagnostic library ("HDF5", "PDB")
 *
 *-----------------------------------------------------------------------------
 */

int c_silo_init_3d(double T, int NCYC, int file_sn, int MEMAD,
				int NXP, int NYP, int NZP,
				double *XS, double *YS, double *ZS,  
				char *SILO_DIAGNOSTIC_LIBRARY)
{
	
	if (silo_initializing_values(T, NCYC, file_sn, MEMAD, SILO_DIAGNOSTIC_LIBRARY) != 0) return 1;
	if (silo_set_paths() != 0) return 1;
	//if (silo_PMPIO_connectivity() != 0) return 1;
	if (silo_mesh_setting_3d(XS, YS, ZS, NXP, NYP, NZP) != 0) return 1;                       
        dimensions = 3;

	return 0;
	
}

/*
 * silo_initializing_values set the initial values from the LEEOR2D to the global struct gsi.
 */
int silo_initializing_values(double T, int NCYC, int file_sn, int MEMAD, 
		char *SILO_DIAGNOSTIC_LIBRARY) {

	/* Values form from Fortran LEEOR2D */

     /* Values form from Fortran LEEOR3D */
	/* Standard MPI initialization stuff */
//	if ( dimensions == 3)
//	{
//		MPI_Comm_size(MPI_COMM_WORLD, &_gsi.mpi_size);
//		MPI_Comm_rank(MPI_COMM_WORLD, &_gsi.mpi_rank);
//	}
	
	/* Setting the Error message string */
	if (!(_gsi.error_message = (char*) malloc(1024 * sizeof(char)))) {
		write_to_bugfile("Memory allocation failure at the allocation of the error_message variable");
		return 1;
	}

	// Taking start times
	_gsi.t_start = time(NULL);
	_gsi.clock_start = clock();

	// Setting default values of global structure
	_gsi.silo_file_ext = strdup(SILO_DIAGNOSTIC_LIBRARY);

	// Setting the silo diagnostic 
	if (!strcmp(SILO_DIAGNOSTIC_LIBRARY, "pdb"))
		_gsi.silo_driver = DB_PDB;
	else
		_gsi.silo_driver = DB_HDF5;

	// Setting values from function arguments
	_gsi.memad = MEMAD;
	_gsi.file_sn = file_sn;
	_gsi.cycle = NCYC;
	_gsi.time = T;
	_gsi.silo_var_count = 0;

	/* Allocate the variables which are defined by the memad */
	if (!(_gsi.coordnames = (char **) malloc(_gsi.memad * sizeof(char*)))) {
		sprintf(_gsi.error_message, "Memory allocation failure. [%s]\n", strerror(errno));
		write_to_bugfile(_gsi.error_message);
		free(_gsi.error_message);
		return 1;
	}

	if (!(_gsi.coords = (double**) malloc(_gsi.memad * sizeof(double*)))) {
		sprintf(_gsi.error_message, "Memory allocation failure. [%s]\n", strerror(errno));
		write_to_bugfile(_gsi.error_message);
		free(_gsi.error_message);
		return 1;
	}

	return 0;

}

/*
 * silo_set_paths set the paths and creates the directories needed for the diagnostics.
 */
int silo_set_paths(){

	/* Setting the current directory for the cycle silo files */
	if (!(_gsi.directory_path_name = (char*) malloc(1024 * sizeof(char)))) {
		sprintf(_gsi.error_message, "Memory allocation failure. [%s]\n", strerror(errno));
		write_to_bugfile(_gsi.error_message);
		free(_gsi.error_message);
		return 1;
	}


	/* Set the files path for silo diagnostics */
	c_silo_set_files_path(getcwd(NULL, 0));

	sprintf(_gsi.directory_path_name, "%s/Silo_Diagnostics",
			_gsi.current_dir_name);	

	if (fopen(_gsi.directory_path_name, "r") != NULL){
		// Checks if the directory already exsist
	}
	else if ((mkdir(_gsi.directory_path_name, 0777)) == -1) {
		sprintf(_gsi.error_message, "Directory creation failure. error: %s [%s]\n", strerror(errno), _gsi.directory_path_name);
		write_to_bugfile(_gsi.error_message);
		free(_gsi.error_message);
		return 1;
	}

	sprintf(_gsi.directory_path_name, "%s/Silo_Diagnostics/SN.%d",
			_gsi.current_dir_name, _gsi.file_sn);
	if (fopen(_gsi.directory_path_name, "r") != NULL){
		// Checks if the directory already exsist
	}
	else if ((mkdir(_gsi.directory_path_name, 0777)) == -1) {
		sprintf(_gsi.error_message, "Directory creation failure. error: %s [%s]\n", strerror(errno), _gsi.directory_path_name);
		write_to_bugfile(_gsi.error_message);
		free(_gsi.error_message);
		return 1;
	}

	chdir(_gsi.directory_path_name);

	return 0;

}


/*
 * silo_mesh_setting sets the data needed to the creation of the Silo mesh and creates it.
 */
int silo_mesh_setting_2d(double *XS, double *YS, int NXP, int NYP){
	return silo_mesh_setting_3d(XS, YS, NULL, NXP, NYP, 0);
}
int silo_mesh_setting_3d(double *XS, double *YS, double *ZS, int NXP, int NYP, int NZP){

	/* Construct names for the silo files and the directories in them */
	if (!(_gsi.silo_file_name = (char*) malloc(1024 * sizeof(char)))) {
		sprintf(_gsi.error_message, "Memory allocation failure. [%s]\n", strerror(errno));
		write_to_bugfile(_gsi.error_message);
		free(_gsi.error_message);
		return 1;
	}

	if (!(_gsi.silo_var_names = (char**) malloc(1024 * sizeof(char*)))) {
		sprintf(_gsi.error_message, "Memory allocation failure. [%s]\n", strerror(errno));
		write_to_bugfile(_gsi.error_message);
		free(_gsi.error_message);
		return 1;
	}

	sprintf(_gsi.silo_file_name, "silo_%03d.%s", _gsi.file_sn, _gsi.silo_file_ext);

	/* Open the Silo file */
	if (!(_gsi.silo_file = DBCreate(_gsi.silo_file_name, DB_CLOBBER, DB_LOCAL, NULL, _gsi.silo_driver))) {
		sprintf(_gsi.error_message, "Mesh allocation failure. [%s]\n", strerror(errno));
		write_to_bugfile(_gsi.error_message);
		free(_gsi.error_message);
		return 1;
	}

	/* Node coordinates */

	switch (_gsi.memad) {
	case 1:
		_gsi.coordnames[0] = strdup("Xcoords");
		_gsi.coords[0] = XS;
		break;
	case 2:
		_gsi.coordnames[0] = strdup("Xcoords");
		_gsi.coordnames[1] = strdup("Ycoords");
		_gsi.coords[0] = XS;
		_gsi.coords[1] = YS;
		break;
	case 3:
		_gsi.coordnames[0] = strdup("Xcoords");
		_gsi.coordnames[1] = strdup("Ycoords");
		_gsi.coordnames[2] = strdup("Zcoords");
		_gsi.coords[0] = XS;
		_gsi.coords[1] = YS;
		_gsi.coords[2] = ZS;
		break;
	default:
		sprintf(_gsi.error_message, "Error wrong memad");
		write_to_bugfile(_gsi.error_message);
		free(_gsi.error_message);

		exit(1);
	}

	/* Connectivity */ 

	if(!(_gsi.dims_sizes = (int*) malloc(_gsi.memad * sizeof(int)))){
		sprintf(_gsi.error_message, "Memory allocation failure. [%s]\n", strerror(errno));
		write_to_bugfile(_gsi.error_message);
		free(_gsi.error_message);
		return 1;
	}

	switch (_gsi.memad) {
	case 1: 
		_gsi.dims_sizes[0] = NXP;
		_gsi.nv = _gsi.dims_sizes[0];
		break;
	case 2: 
		_gsi.dims_sizes[0] = NXP;
		_gsi.dims_sizes[1] = NYP;
		_gsi.nv = _gsi.dims_sizes[0] * _gsi.dims_sizes[1];
		break;
	case 3:
		_gsi.dims_sizes[0] = NXP;
		_gsi.dims_sizes[1] = NYP;
		_gsi.dims_sizes[2] = NZP;
		_gsi.nv = _gsi.dims_sizes[0] * _gsi.dims_sizes[1] * _gsi.dims_sizes[2];
		break;
	default :
		sprintf(_gsi.error_message, "Error wrong memad");
		write_to_bugfile(_gsi.error_message);
		free(_gsi.error_message);
		exit(1);
	}

	_gsi.silo_opt_list = DBMakeOptlist(2);

	DBAddOption(_gsi.silo_opt_list, DBOPT_DTIME, &_gsi.time);
	DBAddOption(_gsi.silo_opt_list, DBOPT_CYCLE, &_gsi.cycle);

	/* Write the mesh geometry into the silo file */
	if (DBPutQuadmesh(_gsi.silo_file, "mesh", _gsi.coordnames, _gsi.coords, _gsi.dims_sizes, 
			_gsi.memad, DB_DOUBLE, DB_NONCOLLINEAR, _gsi.silo_opt_list)) {

		sprintf(_gsi.error_message, "DBPutQuadmesh at silo_mesh_settings. [%s]\n", strerror(errno));
		write_to_bugfile(_gsi.error_message);
		free(_gsi.error_message);
		return 1;
	}

	_gsi.rank_total_memory_write += _gsi.nv * sizeof(double);


	int i;
	for (i = 0; i < _gsi.memad; i++) {
		free(_gsi.coordnames[i]);
	}

	free(_gsi.coords);
	free(_gsi.coordnames);

	return 0;

}

/* ############################################################################# */

/* 
 *  c_silo_write_node_data_int writes an int values of the diagnostic to on the node-center mesh
 */
int c_silo_write_node_data_int(char *var_name, int *values) {

	if ((DBPutQuadvar1(_gsi.silo_file, var_name, "mesh", values,
			_gsi.dims_sizes, _gsi.memad, NULL, 0, DB_INT, DB_NODECENT, _gsi.silo_opt_list))) {

		sprintf(_gsi.error_message, "DBPutQuadvar1 at c_silo_write_node_data_int failure. [%s]\n", strerror(errno));
		write_to_bugfile(_gsi.error_message);
		free(_gsi.error_message);
		return 1;
	}
	_gsi.silo_var_names[_gsi.silo_var_count++] = strdup(var_name);

	_gsi.rank_total_memory_write += _gsi.nv * sizeof(int);

	return 0;

}

/* 
 *  c_silo_write_node_data_float writes an float values of the diagnostic to on the node-center mesh
 */
int c_silo_write_node_data_float(char *var_name, float *values) {

	if ((DBPutQuadvar1(_gsi.silo_file, var_name, "mesh", values, 
			_gsi.dims_sizes, _gsi.memad, NULL, 0, DB_FLOAT, DB_NODECENT, _gsi.silo_opt_list))) {

		sprintf(_gsi.error_message, "DBPutQuadvar1 at c_silo_write_node_data_float failure. [%s]\n", strerror(errno));
		write_to_bugfile(_gsi.error_message);
		free(_gsi.error_message);
		return 1;
	}
	_gsi.silo_var_names[_gsi.silo_var_count++] = strdup(var_name);

	_gsi.rank_total_memory_write += _gsi.nv * sizeof(float);

	return 0;

}

/* 
 *  c_silo_write_node_data_double writes an double values of the diagnostic to on the node-center mesh
 */
int c_silo_write_node_data_double(char *var_name, double *values) {

	if ((DBPutQuadvar1(_gsi.silo_file, var_name, "mesh", values,
			_gsi.dims_sizes, _gsi.memad, NULL, 0, DB_DOUBLE, DB_NODECENT, _gsi.silo_opt_list))) {

		sprintf(_gsi.error_message, "DBPutQuadvar1 at c_silo_write_node_data_double failure. [%s]\n", strerror(errno));
		write_to_bugfile(_gsi.error_message);
		free(_gsi.error_message);
		return 1;
	}
	_gsi.silo_var_names[_gsi.silo_var_count++] = strdup(var_name);

	_gsi.rank_total_memory_write += _gsi.nv * sizeof(double);

	return 0;

}


/* 
 *  c_silo_write_node_data_double_vector writes an double values of the diagnostic to on the node-center mesh
 */
int c_silo_write_node_data_double_vector(char *var_name, double *values, double *values2, double *values3) {

	int i,j;
	int counter =0;
	char **varnames;
	double **comp;
	
	/* Allocate the variables which are defined by the memad */
	if (!(varnames = (char **) malloc(_gsi.memad * sizeof(char*)))) {
		sprintf(_gsi.error_message, "Memory allocation failure. [%s]\n", strerror(errno));
		write_to_bugfile(_gsi.error_message);
		free(_gsi.error_message);
		return 1;
	}
	
	/* Allocate the variables which are defined by the memad */
	if (!(comp = (double **) malloc(_gsi.memad * sizeof(double*)))) {
		sprintf(_gsi.error_message, "Memory allocation failure. [%s]\n", strerror(errno));
		write_to_bugfile(_gsi.error_message);
		free(_gsi.error_message);
		return 1;
	}
	
	switch (_gsi.memad) {
	case 1:
		
		varnames[0] = strdup("Ucoords");
		
		comp[0] = (double *) malloc (sizeof(double) * _gsi.nv);
		
		for(i=0; i < _gsi.nv; i++){
			comp[0][i] = (double) values[i];
		}
		
		break;
		
	case 2:
		
		varnames[0] = strdup("Ucoords");
		varnames[1] = strdup("Vcoords");
		
		comp[0] = (double *) malloc (sizeof(double) * _gsi.nv);
		comp[1] = (double *) malloc (sizeof(double) * _gsi.nv);
		
		for(i=0; i < _gsi.nv; i++){
			comp[0][i] = (double) values[i];
			comp[1][i] = (double) values2[i];
		}
		
		break;
		
	case 3:
		
		varnames[0] = strdup("Ucoords");
		varnames[1] = strdup("Vcoords");
		varnames[2] = strdup("Wcoords");
		
		comp[0] = (double *) malloc (sizeof(double) * _gsi.nv);
		comp[1] = (double *) malloc (sizeof(double) * _gsi.nv);
		comp[2] = (double *) malloc (sizeof(double) * _gsi.nv);
		
		for(i=0; i < _gsi.nv; i++){
			comp[0][i] = (double) values[i];
			comp[1][i] = (double) values2[i];
			comp[1][i] = (double) values3[i];
		}	
		
		break;
		
	default:
		sprintf(_gsi.error_message, "Error wrong memad");
		write_to_bugfile(_gsi.error_message);
		free(_gsi.error_message);
		exit(1);
	}

	if ((DBPutQuadvar(_gsi.silo_file, var_name, "mesh", _gsi.memad, varnames, comp, 
			_gsi.dims_sizes, _gsi.memad, NULL, 0, DB_DOUBLE, DB_NODECENT, _gsi.silo_opt_list))) {

		sprintf(_gsi.error_message, "DBPutQuadvar at c_silo_write_node_data_double_vector failure. [%s]\n", strerror(errno));
		write_to_bugfile(_gsi.error_message);
		free(_gsi.error_message);
		return 1;
	}
	
	_gsi.silo_var_names[_gsi.silo_var_count++] = strdup(var_name);

	_gsi.rank_total_memory_write += _gsi.nv * sizeof(double);
	
	for (i = 0; i < _gsi.memad; i++) {
		free(comp[i]);
		free(varnames[i]);
	}
	free(comp);

	return 0;

}

/* 
 *  c_silo_write_zone_data_int writes an int values of the diagnostic to on the zone-center mesh
 */
int c_silo_write_zone_data_int(char *var_name, int *values) {

	int *dims_sizes1, i;
	dims_sizes1 = malloc(dimensions * sizeof(int));

	for ( i=0 ; i<dimensions ; ++i)
		dims_sizes1[i] = _gsi.dims_sizes[i] -1;

	if ((DBPutQuadvar1(_gsi.silo_file, var_name, "mesh", values,
			dims_sizes1, _gsi.memad, NULL, 0, DB_INT, DB_ZONECENT, _gsi.silo_opt_list))) {

		sprintf(_gsi.error_message, "DBPutQuadvar1 at c_silo_write_zone_data_int failure. [%s]\n", strerror(errno));
		write_to_bugfile(_gsi.error_message);
		free(_gsi.error_message);
		return 1;
	}
	_gsi.silo_var_names[_gsi.silo_var_count++] = strdup(var_name);

	_gsi.rank_total_memory_write += _gsi.nv * sizeof(int);

	return 0;

}


/* 
 *  c_silo_write_zone_data_float writes an float values of the diagnostic to on the zone-center mesh
 */
int c_silo_write_zone_data_float(char *var_name, float *values) {

	int *dims_sizes1, i;
	dims_sizes1 = malloc(dimensions * sizeof(int));

	for ( i=0 ; i<dimensions ; ++i)
		dims_sizes1[i] = _gsi.dims_sizes[i] -1;

	if ((DBPutQuadvar1(_gsi.silo_file, var_name, "mesh", values,
			dims_sizes1, _gsi.memad, NULL, 0, DB_FLOAT, DB_ZONECENT, _gsi.silo_opt_list))) {

		sprintf(_gsi.error_message, "DBPutQuadvar1 at c_silo_write_zone_data_float failure. [%s]\n", strerror(errno));
		write_to_bugfile(_gsi.error_message);
		free(_gsi.error_message);
		return 1;
	}
	_gsi.silo_var_names[_gsi.silo_var_count++] = strdup(var_name);

	_gsi.rank_total_memory_write += _gsi.nv * sizeof(float);

	return 0;

}

/* 
 *  c_silo_write_zone_data_double writes an double values of the diagnostic to on the zone-center mesh
 */
int c_silo_write_zone_data_double(char *var_name, double *values) {

	int *dims_sizes1, i;
	dims_sizes1 = malloc(dimensions * sizeof(int));

	for (  i=0 ; i<dimensions ; ++i)
		dims_sizes1[i] = _gsi.dims_sizes[i] -1;

	if ((DBPutQuadvar1(_gsi.silo_file, var_name, "mesh", values,
			dims_sizes1, _gsi.memad, NULL, 0, DB_DOUBLE, DB_ZONECENT, _gsi.silo_opt_list))) {

		sprintf(_gsi.error_message, "DBPutQuadvar1 at c_silo_write_zone_data_double failure. [%s]\n", strerror(errno));
		write_to_bugfile(_gsi.error_message);
		free(_gsi.error_message);
		return 1;
	}
	_gsi.silo_var_names[_gsi.silo_var_count++] = strdup(var_name);

	_gsi.rank_total_memory_write += _gsi.nv * sizeof(double);

	return 0;

}


///* 
// *  c_silo_write_zone_data_double writes an double values of the diagnostic to on the zone-center mesh
// */
//int c_silo_write_zone_data_double_vector(char *var_name, double *values) {
//
//	int dims_sizes1[2];	
//	char **varnames;
//	
//	if (!(varnames = (char **) malloc(_gsi.memad * sizeof(char*)))) {
//		sprintf(_gsi.error_message, "Memory allocation failure. [%s]\n", strerror(errno));
//		write_to_bugfile(_gsi.error_message);
//		free(_gsi.error_message);
//		return 1;
//	}
//
//	varnames[0] = strdup("var1");
//	varnames[1] = strdup("var2");	
//
//	dims_sizes1[0] = _gsi.dims_sizes[0] -1;
//	dims_sizes1[1] = _gsi.dims_sizes[1] -1;
//
//	if ((DBPutQuadvar(_gsi.silo_file, var_name, "mesh", 2, varnames, values, 
//			dims_sizes1, _gsi.memad, NULL, 0, DB_DOUBLE, DB_ZONECENT, _gsi.silo_opt_list))){
//
//		sprintf(_gsi.error_message, "DBPutQuadvar at c_silo_write_zone_data_double_vector failure. [%s]\n", strerror(errno));
//		write_to_bugfile(_gsi.error_message);
//		free(_gsi.error_message);
//		return 1;
//	}
//	_gsi.silo_var_names[_gsi.silo_var_count++] = strdup(var_name);
//
//	_gsi.rank_total_memory_write += _gsi.nv * sizeof(double);
//
//	return 0;
//
//}


/* ############################################################################# */

/*
 * c_silo_finalize close the files and creates the movie links
 */
int c_silo_finalize() {

	char *movie_file_name;
	char *temp_current_directory;

	int sum_nc, sum_nv;
	float max_CPU_time, max_Clock_time, temp;
	double sum_memory;

	// Creating soft links from the movie file to the root silo file for each frame

	int ret;

	if (!(movie_file_name = (char*) malloc(1024 * sizeof(char)))) {
		sprintf(_gsi.error_message, "Memory allocation failure. [%s]\n", strerror(errno));
		write_to_bugfile(_gsi.error_message);
		free(_gsi.error_message);
		return 1;
	}
	if (!(temp_current_directory = (char*) malloc(1024 * sizeof(char)))) {
		sprintf(_gsi.error_message, "Memory allocation failure. [%s]\n", strerror(errno));
		write_to_bugfile(_gsi.error_message);
		free(_gsi.error_message);
		return 1;
	}


	sprintf(movie_file_name, "silo_%03d.%s", _gsi.file_sn, _gsi.silo_file_ext);
	sprintf(_gsi.directory_path_name, "%s/Silo_Diagnostics/SN.%d/%s",
			_gsi.current_dir_name, _gsi.file_sn, _gsi.silo_file_name);
	strcpy(temp_current_directory, _gsi.current_dir_name);
	strcat(temp_current_directory, "/Silo_Diagnostics/Silo_Movie_Files");

	if ((mkdir(temp_current_directory, 0777)) == -1) {
		if (fopen(temp_current_directory, "r") == NULL){
			sprintf(_gsi.error_message, "Directory creation failure. error: %s [%s]\n", strerror(errno), _gsi.directory_path_name);
			write_to_bugfile(_gsi.error_message);
			free(_gsi.error_message);
			return 1;
		}
	}

	chdir(temp_current_directory);
	ret = symlink(_gsi.directory_path_name, movie_file_name);
	free(movie_file_name);
	free(temp_current_directory);

	sprintf(_gsi.directory_path_name, "%s/Silo_Diagnostics/SN.%d",
			_gsi.current_dir_name, _gsi.file_sn);
	chdir(_gsi.directory_path_name);
	free(_gsi.silo_file_name);
	free(_gsi.silo_file_ext);

	DBFreeOptlist(_gsi.silo_opt_list);
	DBClose(_gsi.silo_file);

	/* End timing */
	_gsi.t_end = time(NULL);
	_gsi.clock_end = clock();

	max_CPU_time = (float) (_gsi.t_end - _gsi.t_start);
	max_Clock_time = (float) (_gsi.clock_end - _gsi.clock_start);

	/* log file */
	sprintf(_gsi.directory_path_name, "%s/Silo_Diagnostics/SN.%d",
			_gsi.current_dir_name, _gsi.file_sn);
	chdir(_gsi.directory_path_name);

	/* Log files of the writing data and the metadata */
	_gsi.rank_writing_data = fopen("rank_writing_data.txt", "a");

	/* Writing data about the rank to "ranks_meta_data.txt" */

	fprintf(_gsi.rank_writing_data, "Cycle:%03d, ", _gsi.cycle);
	fprintf(_gsi.rank_writing_data, "Vertex:%03d, ", _gsi.nv);
	fprintf(_gsi.rank_writing_data, "Storage:%f (MB), ",
			(_gsi.rank_total_memory_write / (1024 * 1024)));
	fprintf(_gsi.rank_writing_data, "CPU_Time:%f, ", (max_CPU_time));
	fprintf(_gsi.rank_writing_data, "Clock_Time:%f (seconds)\n",
			(max_Clock_time / CLOCKS_PER_SEC));
	fclose(_gsi.rank_writing_data);

	/* Deallocation */
	free(_gsi.silo_var_names);
	free(_gsi.directory_path_name);

	/* Return to root directory */
	chdir(_gsi.current_dir_name);
	free(_gsi.current_dir_name);

	return 0;
}

/*
 * write_to_bugfile() writes a error message to fortran bug file
 * In leeor2d there is no bugfile, therefor we use the out file as the output file for the write_to_bugfile method
 */
int write_to_bugfile(char* error_message){

	FILE *bugfile;
	char file_name[20];

	bugfile = fopen("out", "a");

	fprintf(bugfile, "%s", error_message);
	fclose(bugfile);

	return 0;
}







