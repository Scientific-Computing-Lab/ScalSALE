
module textual_diagnostic_module
    use hydro_step_module, only : hydro_step_t
    use time_module, only : time_t

    implicit none

    type, public :: textual_diagnostic_t
        private
        character (80)               :: file_prefix
        class (hydro_step_t), pointer :: hydro_step
        class (time_t), pointer       :: time
        integer                       :: diagnostics_group_file
        integer                       :: diagnostic_number
        integer                       :: diagnostic_counter

    contains

        procedure, public  :: Init_diagnostic => textual_diagnostic_init

        procedure, public  :: Close_diagnostic => textual_diagnostic_close

        procedure, public :: Apply => textual_diagnostic_apply

    end type textual_diagnostic_t

    interface textual_diagnostic_t
        module procedure Constructor_textual_diagnostic_t
    end interface


contains

    type(textual_diagnostic_t) function Constructor_textual_diagnostic_t(prefix, diag_group_file, rank)
        implicit none
        character(len=*)      :: prefix
        integer, intent(in)        :: diag_group_file
        integer, intent(in)        :: rank
        integer :: i
        character*80 :: tmp
        character*80 :: tmp2
        Constructor_textual_diagnostic_t%file_prefix = prefix
        Constructor_textual_diagnostic_t%diagnostics_group_file = diag_group_file
        Constructor_textual_diagnostic_t%diagnostic_counter = 0

!        if (rank /= 0) then
!            tmp2 = trim(Constructor_textual_diagnostic_t%file_prefix)
!            open (Constructor_textual_diagnostic_t%diagnostics_group_file, status = 'replace', file = trim(tmp2(:len(trim(tmp2)) -4 ) // "_" &
!                // trim(tmp) // ".txt"))
!        else
!            open (Constructor_textual_diagnostic_t%diagnostics_group_file, status = 'replace', file = trim(Constructor_textual_diagnostic_t%file_prefix))
!        end if
    end function Constructor_textual_diagnostic_t

    subroutine textual_diagnostic_init (this, hydro, time, diag_num, rank)
        implicit none

        class (textual_diagnostic_t), intent(in out)              :: this
        type(hydro_step_t), pointer, intent(in)                   :: hydro
        type(time_t), pointer      , intent(in)                   :: time
        integer                     , intent(in)     :: diag_num
        integer                     , intent(in)     :: rank

        character*80 :: tmp
        character*80 :: tmp2

        if (rank > 9) then
            write(tmp, '(I2)') rank
        else
            write(tmp, '(I1)') rank
        end if

        this%hydro_step  => hydro
        this%time  => time
        this%diagnostic_counter = 0

        if (rank /= 0) then
            tmp2 = trim(this%file_prefix)
            open (this%diagnostics_group_file, status = 'replace', file = trim(tmp2(:len(trim(tmp2)) -4 ) // "_" &
                // trim(tmp) // ".txt"))
        else
            open (this%diagnostics_group_file, status = 'replace', file = trim(this%file_prefix))
        end if


    end subroutine textual_diagnostic_init

    subroutine textual_diagnostic_apply (this)
        implicit none
        class (textual_diagnostic_t), intent(in out)               :: this
        integer :: nx, ny, nz, nzp, nyp, nxp
        real(8), dimension(:,:,:), pointer    :: ptr_x, ptr_y, ptr_z
        integer :: i,j,k
1001    format((I5),2(1PE25.15))
1002    format(100000(1PE25.15))
1003    format(100000(I2))

        nx  = this%hydro_step%nx
        ny  = this%hydro_step%ny
        nz  = this%hydro_step%nz

        nxp = this%hydro_step%nxp
        nyp = this%hydro_step%nyp
        nzp = this%hydro_step%nzp
        select case( trim(this%file_prefix) )

            case ('time.txt')
                write(this%diagnostics_group_file,1001), this%time%current_cycle,this%time%time_passed,this%time%dt

            case ('position_x.txt')
                call this%hydro_step%Point_to_mesh_data(ptr_x, ptr_y)

                write(this%diagnostics_group_file,1002), ptr_x(1:nxp,1:nyp,1:nzp)
            case ('position_y.txt')
                call this%hydro_step%Point_to_mesh_data(ptr_x, ptr_y)
                write(this%diagnostics_group_file,1002), ptr_y(1:nxp,1:nyp,1:nzp)

            case ('position_z.txt')
                if (nzp /= 1) then
                    call this%hydro_step%Point_to_mesh_data(ptr_x, ptr_y, ptr_z)
                    write(this%diagnostics_group_file,1002), ptr_z(1:nxp,1:nyp,1:nzp)
                end if
            case ('velocity_x.txt')
                call this%hydro_step%Point_to_velocity_data(ptr_x, ptr_y)
                write(this%diagnostics_group_file,1002), ptr_x(1:nxp,1:nyp,1:nzp)

            case ('velocity_y.txt')
                call this%hydro_step%Point_to_velocity_data(ptr_x, ptr_y)
                write(this%diagnostics_group_file,1002), ptr_y(1:nxp,1:nyp,1:nzp)

            case ('velocity_z.txt')
                if (nzp /= 1) then
                    call this%hydro_step%Point_to_velocity_data(ptr_x, ptr_y, ptr_z)
                    write(this%diagnostics_group_file,1002), ptr_z(1:nxp,1:nyp,1:nzp)
                end if
            case ('pressure.txt')
                call this%hydro_step%Point_to_pressure_data(ptr_x)
                write(this%diagnostics_group_file,1002), ptr_x(1:nx,1:ny,1:nz)

            case ('density.txt')
                call this%hydro_step%Point_to_density_data(ptr_x)
                write(this%diagnostics_group_file,1002), ptr_x(1:nx,1:ny,1:nz)

            case ('sie.txt')
                call this%hydro_step%Point_to_sie_data(ptr_x)
                write(this%diagnostics_group_file,1002), ptr_x(1:nx,1:ny,1:nz)

            case ('index.txt')
                call this%hydro_step%Point_to_mat_id_data(ptr_x)
                write(this%diagnostics_group_file,1002), ptr_x(1:nx,1:ny,1:nz)

            case default
                write (this%diagnostics_group_file,*) 'The specified file name does ' // &
                    'not correspond to any known variable'
        end select

    end subroutine textual_diagnostic_apply

    subroutine textual_diagnostic_close (this)
        implicit none
        class (textual_diagnostic_t), intent(in out)               :: this

        close(this%diagnostics_group_file)
    end subroutine textual_diagnostic_close

end module textual_diagnostic_module
