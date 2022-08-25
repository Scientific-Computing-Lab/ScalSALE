
module silo_diagnostic_module
    use hydro_step_module, only : hydro_step_t
    use time_module, only : time_t
    use silo_fortran_to_c_interface
    use, intrinsic :: iso_c_binding
    implicit none

    type, public :: silo_diagnostic_t
        private
        character (80)               :: file_prefix
        class (hydro_step_t), pointer :: hydro_step
        class (time_t), pointer       :: time
        integer                       :: diagnostics_group_file
        integer                       :: diagnostic_number
        integer                       :: diagnostic_counter

    contains

        procedure, public  :: Init_diagnostic => silo_diagnostic_init

        procedure, public  :: Close_diagnostic => silo_diagnostic_close

        procedure, public :: Apply => silo_diagnostic_apply

    end type silo_diagnostic_t

    interface silo_diagnostic_t
        module procedure Constructor_silo_diagnostic_t
    end interface

contains

      type(silo_diagnostic_t) function Constructor_silo_diagnostic_t()
          implicit none
Constructor_silo_diagnostic_t%diagnostic_number = 0
        end function Constructor_silo_diagnostic_t
  
    subroutine silo_diagnostic_init (this, hydro, time)
        implicit none

        class (silo_diagnostic_t), intent(in out)               :: this
        type(hydro_step_t), pointer, intent(in)                   :: hydro           
        type(time_t), pointer      , intent(in)                   :: time            
   
        this%hydro_step  => hydro
        this%time  => time
        this%diagnostic_counter = 0

    end subroutine silo_diagnostic_init

    subroutine silo_diagnostic_apply (this)
        implicit none
        class (silo_diagnostic_t), intent(in out)               :: this
        character(3) :: num_string_3      
        real(8), dimension(:,:,:), pointer    :: ptr_x, ptr_y, ptr_z 
        real(8), dimension(:,:,:), allocatable   :: tmp_x, tmp_y, tmp_z 
        real(8), dimension(:,:,:), allocatable   :: x_zone, y_zone, z_zone  
        real(8), dimension(:), allocatable    :: tmp_u, tmp_v, tmp_w
        integer :: ierr, counter, i, j, nxp, nyp, nzp,k

        call this%hydro_step%Point_to_mesh_data(ptr_x, ptr_y)
        nxp = size(ptr_x(1:,1,1))-1
        nyp = size(ptr_x(1,1:,1))-1
        nzp = size(ptr_x(1,1,1:))-1
nullify(ptr_x)
        nullify(ptr_y)
        nullify(ptr_z)

        if (nzp < 2) then
            ierr = c_silo_init_2d(this%time%time_passed, &
                this%time%current_cycle, &
                this%diagnostic_counter, &
                2, &
                nxp, &
                nyp, &
                reshape(ptr_x(1:nxp,1:nyp,1),[nxp, nyp]), &
                reshape(ptr_y(1:nxp,1:nyp,1),[nxp, nyp]), &
                trim('pdb') // C_NULL_CHAR)
            allocate(tmp_x(1:nxp-1, 1:nyp-1,1))
            allocate(tmp_y(1:nxp-1, 1:nyp-1,1))

        else
            call this%hydro_step%Point_to_mesh_data(ptr_x, ptr_y, ptr_z)
            allocate(tmp_x(1:nxp, 1:nyp, 1:nzp))
            allocate(tmp_y(1:nxp, 1:nyp, 1:nzp))
            allocate(tmp_z(1:nxp, 1:nyp, 1:nzp))
            do k = 1, nzp
                do j = 1, nyp
                    do i = 1, nxp
                        tmp_x(i,j,k) = ptr_x(i,j,k)
                        tmp_y(i,j,k) = ptr_y(i,j,k)
                        tmp_z(i,j,k) = ptr_z(i,j,k)
                    end do
                end do
            end do
            ierr = c_silo_init_3d(this%time%time_passed, &
                this%time%current_cycle, &
                this%diagnostic_counter, &
                3, &
                nxp, &
                nyp, &
                nzp, &
                tmp_x, &
                tmp_y, &
                tmp_z, &
                trim('pdb') // C_NULL_CHAR)

            allocate(x_zone(1:nxp-1, 1:nyp-1,1:nzp-1))
            allocate(y_zone(1:nxp-1, 1:nyp-1,1:nzp-1))
            allocate(z_zone(1:nxp-1, 1:nyp-1,1:nzp-1))
        end if
nullify(ptr_x)
        nullify(ptr_y)
        nullify(ptr_z)
        call this%hydro_step%Point_to_pressure_data(ptr_x)
        do k = 1, max(1,nzp-1)
            do j = 1, nyp-1
                do i = 1, nxp-1
                    x_zone(i,j,k) = ptr_x(i,j,k)

                end do
            end do
        end do

        ierr = c_silo_write_zone_data_double("P" // C_NULL_CHAR, x_zone)
nullify(ptr_x)
        nullify(ptr_y)
        nullify(ptr_z)
        call this%hydro_step%Point_to_density_data(ptr_x)
        do k = 1, max(1,nzp-1)
            do j = 1, nyp-1
                do i = 1, nxp-1
                    x_zone(i,j,k) = ptr_x(i,j,k)
                end do
            end do
        end do
        ierr = c_silo_write_zone_data_double("RO" // C_NULL_CHAR, x_zone)
nullify(ptr_x)
        nullify(ptr_y)
        nullify(ptr_z)
        call this%hydro_step%Point_to_sie_data(ptr_x)
        do k = 1, max(1,nzp-1)
            do j = 1, nyp-1
                do i = 1, nxp-1
                    x_zone(i,j,k) = ptr_x(i,j,k)
                end do
            end do
        end do
        ierr = c_silo_write_zone_data_double("SIE" // C_NULL_CHAR, x_zone)
nullify(ptr_x)
        nullify(ptr_y)
        nullify(ptr_z)
        call this%hydro_step%Point_to_artificial_viscosity_data(ptr_x)
        do k = 1, max(1,nzp-1)
            do j = 1, nyp-1
                do i = 1, nxp-1
                    x_zone(i,j,k) = ptr_x(i,j,k)
                end do
            end do
        end do
        ierr = c_silo_write_zone_data_double("Q" // C_NULL_CHAR, x_zone)
nullify(ptr_x)
        nullify(ptr_y)
        nullify(ptr_z)
        if (nzp == 1) then
            call this%hydro_step%Point_to_velocity_data(ptr_x, ptr_y)
            do k = 1, max(1,nzp)
                do j = 1, nyp
                    do i = 1, nxp
                        tmp_x(i,j,k) = ptr_x(i,j,k)
                        tmp_y(i,j,k) = ptr_y(i,j,k)
                    end do
                end do
            end do
            ierr = c_silo_write_node_data_double("U" // C_NULL_CHAR, tmp_x)
            ierr = c_silo_write_node_data_double("V" // C_NULL_CHAR, tmp_y)
        end if

        if (nzp > 1) then

        nullify(ptr_x)
        nullify(ptr_y)
        nullify(ptr_z)
            call this%hydro_step%Point_to_velocity_data(ptr_x, ptr_y,ptr_z)
            do k = 1, max(1,nzp)
                do j = 1, nyp
                    do i = 1, nxp
                                            tmp_x(i,j,k) = ptr_x(i,j,k)
                        tmp_y(i,j,k) = ptr_y(i,j,k)
                        tmp_z(i,j,k) = ptr_z(i,j,k)
                    end do
                end do
            end do
                        ierr = c_silo_write_node_data_double("U" // C_NULL_CHAR, tmp_x)
            ierr = c_silo_write_node_data_double("V" // C_NULL_CHAR, tmp_y)
            ierr = c_silo_write_node_data_double("W" // C_NULL_CHAR, tmp_z)
        end if

        deallocate(tmp_x)
        deallocate(tmp_y)
        deallocate(tmp_z)
        deallocate(x_zone)
        deallocate(y_zone)
        deallocate(z_zone)
        ierr = c_silo_finalize()
        this%diagnostic_counter = this%diagnostic_counter + 1
        write (num_string_3, "(I3)") this%diagnostic_counter
1002    format(100000(1PE25.17))



    end subroutine silo_diagnostic_apply

    subroutine silo_diagnostic_close (this)
        implicit none
        class (silo_diagnostic_t), intent(in out)               :: this

    end subroutine silo_diagnostic_close
end module silo_diagnostic_module
