module material_quantity_module
   use data_4d_module, only : data_4d_t
   use cell_boundary_condition_module, only : cell_boundary_condition_t, cell_bc_wrapper_t
   use quantity_module
   use boundary_parameters_module      , only : boundary_parameters_t

   implicit none
   private
   public :: material_quantity_t

   type, extends(quantity_t) :: material_quantity_t
      private

      integer, public :: nmats
      type(cell_bc_wrapper_t), dimension(:), pointer, public :: boundary_conditions

   contains

!      procedure, public :: Init_material_quantity_init_val

!      procedure, public :: Init_material_quantity_copy

      procedure, public :: Clean_material_quantity

      procedure, public :: Apply_boundary

   end type

    interface material_quantity_t
      module procedure Constructor
      module procedure Constructor1

   end interface material_quantity_t

contains

   type(material_quantity_t) function Constructor(initial_val, d1, d2, d3,d4, bc, bc_params)
      implicit none
      real(8)                                , intent(in)     :: initial_val
      integer                                , intent(in)     :: d1
      integer                                , intent(in)     :: d2

      integer                                , intent(in)     :: d3
      integer                                , intent(in)     :: d4
      type(cell_bc_wrapper_t), dimension(:), pointer, intent(in) :: bc
      type(boundary_parameters_t), pointer, intent(in) :: bc_params
      integer                                                 :: i

      call Constructor%Init_quantity_init_val_4d (initial_val, d1, d2, d3,d4, 1, bc_params)
      Constructor%boundary_conditions => bc
      Constructor%nmats = d4
   end function


   type(material_quantity_t) function Constructor1( initial_val, d1, d2, d3,d4)
      implicit none
      real(8)                                , intent(in)     :: initial_val
      integer                                , intent(in)     :: d1
      integer                                , intent(in)     :: d2

      integer                                , intent(in)     :: d3
      integer                                , intent(in)     :: d4
!      type(cell_bc_wrapper_t), dimension(:), pointer, intent(in) :: bc
!      type(boundary_parameters_t), pointer, intent(in) :: bc_params
      integer                                                 :: i

      call Constructor1%Init_quantity_no_bc (initial_val, d1, d2, d3,d4, 1)
!      Constructor1%boundary_conditions => bc
      Constructor1%nmats = d4
   end function


   subroutine Init_material_quantity_copy(this, copy)
      implicit none
      class (material_quantity_t), intent(in out) :: this
      class (material_quantity_t), intent(in)     :: copy

   end subroutine

   subroutine Clean_material_quantity (this)
      class (material_quantity_t), intent(in out) :: this

      deallocate (this%boundary_conditions)
      call this%Clean_quantity ()
   end subroutine Clean_material_quantity


   subroutine Apply_boundary(this, is_blocking)
      implicit none
      class (material_quantity_t)   , intent(in out) :: this
      logical, optional :: is_blocking

      integer                                    :: i, edge
      integer, dimension(:)     , pointer        :: ed_num
      logical :: is_blocking_local

      if (.not. present(is_blocking)) then
         is_blocking_local = .True.
      else
         is_blocking_local = is_blocking
      end if


      call this%boundary_params%Point_to_edges (ed_num)
      do i = 1, size(ed_num)
         edge = ed_num(i)
         call this%boundary_conditions(edge)%bc%Calculate (this, edge)
      end do

      if (is_blocking_local) then
         call this%Exchange_virtual_space_blocking()
      else
         call this%Exchange_virtual_space_nonblocking()
      end if

   end subroutine Apply_boundary


end module material_quantity_module
