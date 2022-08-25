
module pressure_module
    use cell_quantity_module, only : cell_quantity_t
    use cell_boundary_condition_module, only : cell_boundary_condition_t, cell_bc_wrapper_t
    use boundary_parameters_module      , only : boundary_parameters_t

    implicit none
    private
    public :: pressure_t

    type, extends(cell_quantity_t) :: pressure_t
        private


    contains

        procedure, public :: Clean_pressure

        procedure, public :: Write_pressure
        procedure, public :: Write_quantity_abstract => Write_pressure

        procedure, public :: Read_pressure
        procedure, public :: Read_quantity_abstract => Read_pressure

    end type


    interface pressure_t

        module procedure Constructor
    end interface pressure_t

contains

    type(pressure_t) function Constructor(initial_val, d1, d2, d3, bc, bc_params)
        real(8)                                           , intent(in) :: initial_val 
        integer                                           , intent(in) :: d1           
        integer                                           , intent(in) :: d2           
        integer                                           , intent(in) :: d3           
        type(cell_bc_wrapper_t), dimension(:), pointer    , intent(in) :: bc           
        type(boundary_parameters_t), pointer, intent(in) :: bc_params

        call Constructor%Init_cell_quantity_init_val (initial_val, d1, d2, d3, bc, bc_params)
    end function


    subroutine Clean_pressure (this)
        class (pressure_t), intent(in out) :: this   

        call this%Clean_cell_quantity ()
    end subroutine Clean_pressure

    subroutine Write_pressure(this, unit, iostat, iomsg)
        class (pressure_t), intent(in) :: this  
        integer,      intent(in)     :: unit
        integer,      intent(out)    :: iostat
        character(*), intent(in out)  :: iomsg

#ifdef DEBUG
        write(*,*) '@@@ in Write_pressure @@@'
#endif

        call this%Write_cell_quantity(unit, iostat=iostat, iomsg=iomsg)

#ifdef DEBUG
        write(*,*) '@@@ end Write_pressure @@@'
#endif

    end subroutine Write_pressure

    subroutine Read_pressure(this, unit, iostat, iomsg)
        class (pressure_t), intent(in out) :: this  
        integer,      intent(in)     :: unit
        integer,      intent(out)    :: iostat
        character(*), intent(in out)  :: iomsg

#ifdef DEBUG
        write(*,*) '@@@ in Read_pressure @@@'
#endif

        call this%Read_cell_quantity(unit, iostat=iostat, iomsg=iomsg)

#ifdef DEBUG
        write(*,*) '@@@ end Read_pressure @@@'
#endif

    end subroutine Read_pressure



end module pressure_module
