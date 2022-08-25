
module acceleration_module
   use vertex_quantity_module, only : vertex_quantity_t
   use vertex_boundary_condition_module, only : vertex_boundary_condition_t, vertex_bc_wrapper_t
   use general_utils_module, only : Get_axises_num
   use boundary_parameters_module      , only : boundary_parameters_t

   implicit none
   private
   public :: acceleration_t

   type, extends(vertex_quantity_t) :: acceleration_t
      private


   contains


      procedure, public :: Clean_acceleration

      procedure, public :: Write_acceleration
      procedure, public :: Write_quantity_abstract => Write_acceleration

      procedure, public :: Read_acceleration
      procedure, public :: Read_quantity_abstract => Read_acceleration

   end type


   interface acceleration_t

      module procedure Constructor
   end interface acceleration_t

contains

   type(acceleration_t) function Constructor(initial_data, d1, d2, d3, bc, bc_params)
      real(8)                                           , intent(in) :: initial_data 
      integer                                           , intent(in) :: d1           
      integer                                           , intent(in) :: d2           
      integer                                           , intent(in) :: d3           
      type (vertex_bc_wrapper_t) , dimension(:), pointer, intent(in out) :: bc           
      type(boundary_parameters_t), pointer, intent(in) :: bc_params

      integer                                                        :: axises_num

      axises_num = Get_axises_num (d1, d2, d3)
      call Constructor%Init_vertex_quantity_init_val (initial_data, d1, d2, d3, axises_num, bc, bc_params)
   end function


   subroutine Clean_acceleration (this)
      class (acceleration_t), intent(in out) :: this 

      call this%Clean_vertex_quantity ()
   end subroutine Clean_acceleration

   subroutine Write_acceleration(this, unit, iostat, iomsg)
      class (acceleration_t), intent(in) :: this  
      integer,      intent(in)     :: unit
      integer,      intent(out)    :: iostat
      character(*), intent(in out)  :: iomsg

#ifdef DEBUG
      write(*,*) '@@@ in Write_acceleration @@@'
#endif

      call this%Write_vertex_quantity(unit, iostat=iostat, iomsg=iomsg)

#ifdef DEBUG
      write(*,*) '@@@ end Write_acceleration @@@'
#endif

   end subroutine Write_acceleration

   subroutine Read_acceleration(this, unit, iostat, iomsg)
      class (acceleration_t), intent(in out) :: this  
      integer,      intent(in)     :: unit
      integer,      intent(out)    :: iostat
      character(*), intent(in out)  :: iomsg

#ifdef DEBUG
      write(*,*) '@@@ in Read_acceleration @@@'
#endif

      call this%Read_vertex_quantity(unit, iostat=iostat, iomsg=iomsg)

#ifdef DEBUG
      write(*,*) '@@@ end Read_acceleration @@@'
#endif

   end subroutine Read_acceleration

end module acceleration_module
