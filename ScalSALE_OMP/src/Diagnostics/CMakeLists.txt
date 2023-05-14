include_directories($ENV{SILO_PATH}/include)
#include_directories($ENV{HDF5_PATH}/include)
link_directories($ENV{SILO_PATH}/lib)
#link_directories($ENV{HDF5_PATH}/lib)

#add_library(Silo SHARED silo_module.c silo_fortran_to_c_interface.f90)
#target_link_libraries(Silo silo)

add_library(Diagnostics SHARED textual_diagnostic.f90 )#silo_diagnostic.f90 )#textual_diagnostic_hdf5.f90 silo_diagnostic.f90 plot_diagnostic.f90)
#target_link_libraries(Diagnostics Silo)
#target_link_libraries(Diagnostics hdf5 )

#diagnostics_group.f90 plot.f90 radius_1d.f90 shock_location.f90)
