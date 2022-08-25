
module temperature_module
   use cell_quantity_module, only : cell_quantity_t
   use cell_boundary_condition_module, only : cell_boundary_condition_t, cell_bc_wrapper_t
   use boundary_parameters_module      , only : boundary_parameters_t

   implicit none
   private
   public :: temperature_t

   type, extends(cell_quantity_t) :: temperature_t
      private


   contains


      procedure, public :: Clean_temperature

      procedure, public :: Write_temperature
      procedure, public :: Write_quantity_abstract => Write_temperature

      procedure, public :: Read_temperature
      procedure, public :: Read_quantity_abstract => Read_temperature

   end type


   interface temperature_t

      module procedure Constructor
   end interface temperature_t

contains

   type(temperature_t) function Constructor(temperature_init, d1, d2, d3, bc, bc_params)

      real(8)                                           , intent(in) :: temperature_init

      integer                                           , intent(in) :: d1 
      integer                                           , intent(in) :: d2 
      integer                                           , intent(in) :: d3 
      type(cell_bc_wrapper_t), dimension(:), pointer    , intent(in) :: bc           
      type(boundary_parameters_t), pointer, intent(in) :: bc_params

      call Constructor%Init_cell_quantity_init_val (temperature_init, d1, d2, d3, bc, bc_params)
   end function


   subroutine Clean_temperature (this)
      class (temperature_t), intent(in out) :: this   

      call this%Clean_cell_quantity ()
   end subroutine Clean_temperature

   subroutine Write_temperature(this, unit, iostat, iomsg)
      class (temperature_t), intent(in) :: this  
      integer,      intent(in)     :: unit
      integer,      intent(out)    :: iostat
      character(*), intent(in out)  :: iomsg

#ifdef DEBUG
      write(*,*) '@@@ in Write_temperature @@@'
#endif

      call this%Write_cell_quantity(unit, iostat=iostat, iomsg=iomsg)

#ifdef DEBUG
      write(*,*) '@@@ end Write_temperature @@@'
#endif

   end subroutine Write_temperature

   subroutine Read_temperature(this, unit, iostat, iomsg)
      class (temperature_t), intent(in out) :: this  
      integer,      intent(in)     :: unit
      integer,      intent(out)    :: iostat
      character(*), intent(in out)  :: iomsg

#ifdef DEBUG
      write(*,*) '@@@ in Read_temperature @@@'
#endif

      call this%Read_cell_quantity(unit, iostat=iostat, iomsg=iomsg)

#ifdef DEBUG
      write(*,*) '@@@ end Read_temperature @@@'
#endif

   end subroutine Read_temperature


end module temperature_module
