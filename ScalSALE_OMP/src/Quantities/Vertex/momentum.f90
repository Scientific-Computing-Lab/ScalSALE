
module momentum_module
   use vertex_quantity_module, only : vertex_quantity_t
   use vertex_boundary_condition_module, only : vertex_boundary_condition_t, vertex_bc_wrapper_t
   use general_utils_module, only : Get_axises_num
   use boundary_parameters_module      , only : boundary_parameters_t

   implicit none
   private
   public :: momentum_t

   type, extends(vertex_quantity_t) :: momentum_t
      private


   contains

      procedure, public :: Clean_momentum
   end type


   interface momentum_t

      module procedure Constructor
   end interface momentum_t

contains

   type(momentum_t) function Constructor(initial_data, nx, ny, nz, bc, bc_params)
      real(8)                                           , intent(in) :: initial_data 
      integer                                           , intent(in) :: nx           
      integer                                           , intent(in) :: ny           
      integer                                           , intent(in) :: nz           
      type (vertex_bc_wrapper_t) , dimension(:), pointer, intent(in) :: bc           
      type(boundary_parameters_t), pointer, intent(in) :: bc_params

      integer                                                        :: axises_num

      axises_num = Get_axises_num (nx, ny, nz)
      call Constructor%Init_vertex_quantity_init_val (initial_data, nx, ny, nz, axises_num, bc, bc_params)
   end function


   subroutine Clean_momentum (this)
      class (momentum_t), intent(in out) :: this 

      call this%Clean_vertex_quantity ()
   end subroutine Clean_momentum


end module momentum_module
