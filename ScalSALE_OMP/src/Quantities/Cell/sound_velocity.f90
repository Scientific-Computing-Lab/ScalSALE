
module sound_velocity_module
   use cell_quantity_module, only : cell_quantity_t
   use cell_boundary_condition_module, only : cell_boundary_condition_t, cell_bc_wrapper_t
   use boundary_parameters_module      , only : boundary_parameters_t

   implicit none
   private
   public :: sound_velocity_t

   type, extends(cell_quantity_t) :: sound_velocity_t
      private


   contains


      procedure, public :: Clean_sound_velocity

      procedure, public :: Write_sound_velocity
      procedure, public :: Write_quantity_abstract => Write_sound_velocity

      procedure, public :: Read_sound_velocity
      procedure, public :: Read_quantity_abstract => Read_sound_velocity

   end type


   interface sound_velocity_t

      module procedure Constructor_init_val
   end interface sound_velocity_t

contains

   type(sound_velocity_t) function Constructor_init_val(initial_val, d1, d2, d3, bc, bc_params)
      real(8)                                           , intent(in) :: initial_val
      integer                                           , intent(in) :: d1           
      integer                                           , intent(in) :: d2           
      integer                                           , intent(in) :: d3           
      type(cell_bc_wrapper_t), dimension(:), pointer    , intent(in) :: bc           
      type(boundary_parameters_t), pointer, intent(in) :: bc_params

      call Constructor_init_val%Init_cell_quantity_init_val (initial_val, d1, d2, d3, bc, bc_params)
   end function Constructor_init_val


   subroutine Clean_sound_velocity (this)
      class (sound_velocity_t), intent(in out) :: this   

      call this%Clean_cell_quantity ()
   end subroutine Clean_sound_velocity

   subroutine Write_sound_velocity(this, unit, iostat, iomsg)
      class (sound_velocity_t), intent(in) :: this  
      integer,      intent(in)     :: unit
      integer,      intent(out)    :: iostat
      character(*), intent(in out)  :: iomsg

#ifdef DEBUG
      write(*,*) '@@@ in Write_sound_velocity @@@'
#endif

      call this%Write_cell_quantity(unit, iostat=iostat, iomsg=iomsg)

#ifdef DEBUG
      write(*,*) '@@@ end Write_sound_velocity @@@'
#endif

   end subroutine Write_sound_velocity

   subroutine Read_sound_velocity(this, unit, iostat, iomsg)
      class (sound_velocity_t), intent(in out) :: this  
      integer,      intent(in)     :: unit
      integer,      intent(out)    :: iostat
      character(*), intent(in out)  :: iomsg

#ifdef DEBUG
      write(*,*) '@@@ in Read_sound_velocity @@@'
#endif

      call this%Read_cell_quantity(unit, iostat=iostat, iomsg=iomsg)

#ifdef DEBUG
      write(*,*) '@@@ end Read_sound_velocity @@@'
#endif

   end subroutine Read_sound_velocity


end module sound_velocity_module
