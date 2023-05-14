
module diagnostic_module
    use hydro_step_module, only : hydro_step_t
    use time_module, only       : time_t
    !use hdf5
    implicit none
    private
    public :: diagnostic_t, diagnostic_wrapper_t, textual_diagnostic_t, plot_diagnostic_t
        !textual_diagnostic_hdf5_t, Static_hdf5_init, Static_hdf5_close

    character(len=17), parameter :: dir_name = "HDF5_Diagnostics/"
    character(len=9),  parameter :: main_file = "hdf5_file.d"

    type, abstract :: diagnostic_t
        private

        class (hydro_step_t), pointer :: hydro_step             
        class (time_t), pointer       :: time                   
        integer                       :: diagnostics_group_file 
        integer                       :: diagnostic_number      
        integer                       :: diagnostic_counter     
    contains


        procedure, public  :: Init_diagnostic

        procedure, public  :: Close_diagnostic

        procedure, public :: Apply



    end type diagnostic_t

    type diagnostic_wrapper_t
        class(diagnostic_t), pointer :: diag
    end type diagnostic_wrapper_t


    type, extends(diagnostic_t) :: textual_diagnostic_t
        private
        character (80)               :: file_prefix            

    contains


        procedure, public  :: Init_diagnostic => textual_diagnostic_init

        procedure, public  :: Close_diagnostic => textual_diagnostic_close

        procedure, public :: Apply => textual_diagnostic_apply

    end type textual_diagnostic_t


    !    type, extends(diagnostic_t) :: textual_diagnostic_hdf5_t
    !        private
    !        character (80)                 :: file_prefix
    !
    !1        integer(hid_t) :: main_file_id
    !        integer(hid_t) :: gid
    !        integer(hid_t) :: dset_id
    !        integer(hid_t) :: dspace_id

    !    contains


    !        procedure, public  :: Init_diagnostic => textual_diagnostic_hdf5_init
    !
    !        procedure, public  :: Close_diagnostic => textual_diagnostic_hdf5_close
    !
    !        procedure, public :: Apply => textual_diagnostic_hdf5_apply
    !
    !    end type textual_diagnostic_hdf5_t



    !   type, extends(diagnostic_t) :: silo_diagnostic_t
    !        private
    !        character (80)                 :: file_prefix
    !
    !    contains
    !
    !
    !        procedure, public  :: Init_diagnostic => silo_diagnostic_init
    !
    !        procedure, public  :: Close_diagnostic => silo_diagnostic_close
    !
    !        procedure, public :: Apply => silo_diagnostic_apply
    !
    !    end type silo_diagnostic_t

    type, extends (diagnostic_t) :: plot_diagnostic_t

        private
        character (80)                 :: file_prefix            

    contains

        procedure, public  :: Init_diagnostic => plot_diagnostic_init

        procedure, public  :: Close_diagnostic => plot_diagnostic_close

        procedure, public :: Apply => plot_diagnostic_apply

    end type plot_diagnostic_t

    interface

        module subroutine textual_diagnostic_apply (this)
            implicit none
 
            class (textual_diagnostic_t), intent(in out)               :: this

        end subroutine textual_diagnostic_apply

        module subroutine textual_diagnostic_init (this, hydro, time, diag_num, rank)
            implicit none
            class (textual_diagnostic_t), intent(in out)               :: this
            type(hydro_step_t), pointer, intent(in)                   :: hydro           
            type(time_t), pointer, intent(in)                         :: time            
            integer                     , intent(in)     :: diag_num        
            integer                     , intent(in)     :: rank        
        end subroutine textual_diagnostic_init

        module subroutine textual_diagnostic_close (this)
            implicit none
            class (textual_diagnostic_t), intent(in out)               :: this
        end subroutine textual_diagnostic_close


        !        module subroutine textual_diagnostic_hdf5_apply (this)
        !            implicit none
        !
        !            class (textual_diagnostic_hdf5_t), intent(in out)               :: this
        !        end subroutine textual_diagnostic_hdf5_apply
        !
        !        module subroutine textual_diagnostic_hdf5_init (this, hydro, time, diag_num, main_file_id)
        !            implicit none
        !            class (textual_diagnostic_hdf5_t), intent(in out)               :: this
        !            type(hydro_step_t), pointer, intent(in)                   :: hydro
        !            type(time_t), pointer, intent(in)                         :: time
        !            integer                     , intent(in)     :: diag_num
        !            integer(hid_t)              , intent(in)     :: main_file_id
        !        end subroutine textual_diagnostic_hdf5_init

        !        module subroutine textual_diagnostic_hdf5_close (this)
        !            implicit none
        !            class (textual_diagnostic_hdf5_t), intent(in out)               :: this
        !        end subroutine textual_diagnostic_hdf5_close



        !       module subroutine silo_diagnostic_apply (this)
        !           implicit none
        !
        !            class (silo_diagnostic_t), intent(in out)               :: this
        !        end subroutine silo_diagnostic_apply

        !        module subroutine silo_diagnostic_init (this, hydro, time, diag_num)
        !            implicit none
        !            class (silo_diagnostic_t), intent(in out)               :: this
        !            type(hydro_step_t), pointer, intent(in)                   :: hydro
        !            type(time_t), pointer, intent(in)                         :: time
        !            integer                     , intent(in)     :: diag_num
        !        end subroutine silo_diagnostic_init

        !        module subroutine silo_diagnostic_close (this)
        !            implicit none
        !            class (silo_diagnostic_t), intent(in out)               :: this
        !        end subroutine silo_diagnostic_close

        module subroutine plot_diagnostic_apply (this)
            implicit none

            class (plot_diagnostic_t), intent(in out)               :: this
        end subroutine plot_diagnostic_apply

        module subroutine plot_diagnostic_init (this, hydro, time, diag_num)
            implicit none
            class (plot_diagnostic_t), intent(in out)               :: this
            type(hydro_step_t), pointer, intent(in)                   :: hydro           
            type(time_t), pointer, intent(in)                         :: time            
            integer                     , intent(in)     :: diag_num        
        end subroutine plot_diagnostic_init

        module subroutine plot_diagnostic_close (this)
            implicit none
            class (plot_diagnostic_t), intent(in out)               :: this
        end subroutine plot_diagnostic_close

    end interface


    interface textual_diagnostic_t
        module procedure Constructor_textual_diagnostic_t
    end interface

    !    interface textual_diagnostic_hdf5_t
    !        module procedure Constructor_textual_diagnostic_hdf5_t
    !    end interface

    !    interface silo_diagnostic_t
    !        module procedure Constructor_silo_diagnostic_t
    !    end interface

    interface plot_diagnostic_t
        module procedure Constructor_plot_diagnostic_t
    end interface
contains

    subroutine Init_diagnostic (this, hydro, time, diag_num, rank)
        class(diagnostic_t)         , intent(in out) :: this            
        type(hydro_step_t), pointer, intent(in)      :: hydro           
        type(time_t), pointer, intent(in)            :: time            
        integer                     , intent(in)     :: diag_num
        integer                     , intent(in)     :: rank
!
!        character*80 :: tmp
!        character*80 :: tmp2
!
!        if (rank > 9) then
!            write(tmp, '(I2)') rank
!        else
!            write(tmp, '(I1)') rank
!        end if
!
!        this%hydro_step  => hydro
!        this%time  => time
!        this%diagnostic_counter = 0
!        write(*,*)"FILE PREFIX", trim(this%file_prefix)
!
!        if (rank /= 0) then
!            tmp2 = trim(this%file_prefix)
!            open (this%diagnostics_group_file, status = 'replace', file = trim(tmp2(:len(trim(tmp2)) -4 ) // "_" &
!                // trim(tmp) // ".txt"))
!        else
!            open (this%diagnostics_group_file, status = 'replace', file = trim(this%file_prefix))
!        end if

        this%hydro_step  => hydro
        this%time  => time
        this%diagnostic_number = diag_num
    end subroutine


    subroutine Apply (this)
        class (diagnostic_t), intent(in out) :: this

    end subroutine

    subroutine Close_diagnostic (this)
        class (diagnostic_t), intent(in out) :: this         
    end subroutine

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

        if (rank /= 0) then
            tmp2 = trim(Constructor_textual_diagnostic_t%file_prefix)
            open (Constructor_textual_diagnostic_t%diagnostics_group_file, status = 'replace', file = trim(tmp2(:len(trim(tmp2)) -4 ) // "_" &
                // trim(tmp) // ".txt"))
        else
            open (Constructor_textual_diagnostic_t%diagnostics_group_file, status = 'replace', file = trim(Constructor_textual_diagnostic_t%file_prefix))
        end if
    end function Constructor_textual_diagnostic_t

    !    type(textual_diagnostic_hdf5_t) function Constructor_textual_diagnostic_hdf5_t(prefix, diag_group_file)
    !        implicit none
    !        character(len=*)      :: prefix
    !        integer, intent(in)        :: diag_group_file
    !        integer :: i

    !       Constructor_textual_diagnostic_hdf5_t%file_prefix = prefix
    !       Constructor_textual_diagnostic_hdf5_t%diagnostics_group_file = diag_group_file
    !   end function Constructor_textual_diagnostic_hdf5_t

    !  type(silo_diagnostic_t) function Constructor_silo_diagnostic_t()
    !      implicit none
    !
    !    end function Constructor_silo_diagnostic_t

    type(plot_diagnostic_t) function Constructor_plot_diagnostic_t()
        implicit none

    end function Constructor_plot_diagnostic_t

 !   function Static_hdf5_init() result (file_id)
 !       implicit none

!        integer(hid_t) :: file_id       
  !      integer :: error

        !call h5open_f(error)

        !call h5fcreate_f(dir_name // main_file, H5F_ACC_TRUNC_F, file_id, error)

   ! end function Static_hdf5_init

   ! function Static_hdf5_close (file_id) result (error)
   !     implicit none
!
!        integer(hid_t) :: file_id       
!        integer     ::   error 
!
!        call h5fclose_f(file_id, error)
!
!        call h5close_f(error)
!
!    end function Static_hdf5_close

end module diagnostic_module
