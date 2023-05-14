
MODULE SILO_FORTRAN_TO_C_INTERFACE

    use, intrinsic :: iso_c_binding

      
    IMPLICIT NONE

    INTERFACE

        INTEGER(C_INT) function c_silo_init_2d(T, NCYC, file_sn, MEMAD,&
            NXP, NYP, XS, YS, SILO_DIAGNOSTIC_LIBRARY) BIND(C, NAME="c_silo_init_2d")
            USE, INTRINSIC :: iso_c_binding
            INTEGER(C_INT), VALUE :: NCYC, file_sn, MEMAD, NXP, NYP
            real(C_DOUBLE), VALUE :: T
            real(C_DOUBLE) :: XS(NXP, NYP), YS(NXP, NYP)
            character(kind=C_CHAR), intent(in) :: SILO_DIAGNOSTIC_LIBRARY(*)
        END FUNCTION c_silo_init_2d
         
        INTEGER(C_INT) function c_silo_init_3d(T, NCYC, file_sn, MEMAD,&
            NXP, NYP, NZP, XS, YS, ZS, SILO_DIAGNOSTIC_LIBRARY) BIND(C, NAME="c_silo_init_3d")
            USE, INTRINSIC :: iso_c_binding
            INTEGER(C_INT), VALUE :: NCYC, file_sn, MEMAD, NXP, NYP, NZP
            real(C_DOUBLE), VALUE :: T
            real(C_DOUBLE) :: XS(NXP, NYP, NZP), YS(NXP, NYP, NZP), ZS(NXP, NYP, NZP)
            character(kind=C_CHAR), intent(in) :: SILO_DIAGNOSTIC_LIBRARY(*)
        END FUNCTION c_silo_init_3d

        INTEGER(C_INT) function c_silo_write_node_data_int(var_name, var_values) &
            BIND(C, NAME="c_silo_write_node_data_int")
            USE, INTRINSIC :: iso_c_binding
            implicit none
            character(kind=C_CHAR), intent(in) :: var_name(*)
            INTEGER(C_INT) :: var_values(*)
        END FUNCTION c_silo_write_node_data_int

        INTEGER(C_INT) function c_silo_write_node_data_float(var_name, var_values) &
            BIND(C, NAME="c_silo_write_node_data_float")
            USE, INTRINSIC :: iso_c_binding
            implicit none
            character(kind=C_CHAR), intent(in) :: var_name(*)
            real(C_FLOAT) :: var_values(*)
        END FUNCTION c_silo_write_node_data_float

        INTEGER(C_INT) function c_silo_write_node_data_double(var_name, var_values) &
            BIND(C, NAME="c_silo_write_node_data_double")
            USE, INTRINSIC :: iso_c_binding
            implicit none
            character(kind=C_CHAR), intent(in) :: var_name(*)
            real(C_DOUBLE) :: var_values(*)
        END FUNCTION c_silo_write_node_data_double
        
        INTEGER(C_INT) function c_silo_write_node_data_double_vector(var_name, var_values, var_values2, var_values3) &
            BIND(C, NAME="c_silo_write_node_data_double_vector")
            USE, INTRINSIC :: iso_c_binding
            implicit none
            character(kind=C_CHAR), intent(in) :: var_name(*)
            real(C_DOUBLE) :: var_values(*), var_values2(*), var_values3(*)
        END FUNCTION c_silo_write_node_data_double_vector        

        INTEGER(C_INT) function c_silo_write_zone_data_int(var_name, var_values) &
            BIND(C, NAME="c_silo_write_zone_data_int")
            USE, INTRINSIC :: iso_c_binding
            implicit none
            character(kind=C_CHAR), intent(in) :: var_name(*)
            INTEGER(C_INT) :: var_values(*)
        END FUNCTION c_silo_write_zone_data_int

        INTEGER(C_INT) function c_silo_write_zone_data_float(var_name, var_values) &
            BIND(C, NAME="c_silo_write_zone_data_float")
            USE, INTRINSIC :: iso_c_binding
            implicit none
            character(kind=C_CHAR), intent(in) :: var_name(*)
            real(C_FLOAT) :: var_values(*)
        END FUNCTION c_silo_write_zone_data_float

        INTEGER(C_INT) function c_silo_write_zone_data_double(var_name, var_values) &
            BIND(C, NAME="c_silo_write_zone_data_double")
            USE, INTRINSIC :: iso_c_binding
            implicit none
            character(kind=C_CHAR), intent(in) :: var_name(*)
            real(C_DOUBLE) :: var_values(*)
        END FUNCTION c_silo_write_zone_data_double
        
        INTEGER(C_INT) function c_silo_write_zone_data_double_vector(var_name, var_values) &
            BIND(C, NAME="c_silo_write_zone_data_double_vector")
            USE, INTRINSIC :: iso_c_binding
            implicit none
            character(kind=C_CHAR), intent(in) :: var_name(*)
            real(C_DOUBLE) :: var_values(*)
        END FUNCTION c_silo_write_zone_data_double_vector        

        INTEGER(C_INT) function c_silo_finalize() &
            BIND(C, NAME="c_silo_finalize")
            USE, INTRINSIC :: iso_c_binding
        END FUNCTION c_silo_finalize
        
        
        INTEGER(C_INT) function c_silo_set_file_num_vof(file_num) &
            BIND(C, NAME="c_silo_set_file_num_vof")
            USE, INTRINSIC :: iso_c_binding
            implicit none
            INTEGER(C_INT), VALUE :: file_num
        END FUNCTION c_silo_set_file_num_vof
        
        INTEGER(C_INT) function c_silo_init_vof(T, NCYC, file_sn, MEMAD, NC, NV, &
              XS, YS, ZS, VERTEX_ARRAY, SILO_DIAGNOSTIC_LIBRARY) BIND(C, NAME="c_silo_init_vof")
           USE, INTRINSIC :: iso_c_binding
           INTEGER(C_INT), VALUE :: NCYC, file_sn, MEMAD, NC, NV
           real(C_DOUBLE), VALUE :: T
           real(C_DOUBLE) :: XS(NV), YS(NV), ZS(NV)
           INTEGER(C_INT) :: VERTEX_ARRAY(NC*3)
           character(kind=C_CHAR), intent(in) :: SILO_DIAGNOSTIC_LIBRARY(*)
         END FUNCTION c_silo_init_vof
        
        INTEGER(C_INT) function c_silo_finalize_vof() &
            BIND(C, NAME="c_silo_finalize_vof")
            USE, INTRINSIC :: iso_c_binding
            implicit none
        END FUNCTION c_silo_finalize_vof

    END INTERFACE

     
END MODULE SILO_FORTRAN_TO_C_INTERFACE
