
module cell_quantity_module
   use data_module, only : data_t
   use cell_boundary_condition_module, only : cell_boundary_condition_t, cell_bc_wrapper_t
   use quantity_module
   use boundary_parameters_module      , only : boundary_parameters_t

   implicit none
   private
   public :: cell_quantity_t

   type, abstract, extends(quantity_t) :: cell_quantity_t
      private


      type(cell_bc_wrapper_t), dimension(:), pointer, public :: boundary_conditions

   contains


      procedure, public :: Init_cell_quantity_init_arr

      procedure, public :: Init_cell_quantity_init_val

      procedure, public :: Init_cell_quantity_no_init

      procedure, public :: Init_cell_quantity_copy

      procedure, public :: Clean_cell_quantity

      procedure, public :: Apply_boundary

      procedure, public :: Write_cell_quantity

      procedure, public :: Read_cell_quantity

   end type

contains

   subroutine Init_cell_quantity_init_arr(this, initial_data, d1, d2, d3, bc, bc_params)
      implicit none
      class (cell_quantity_t)                , intent(in out) :: this          
      real(8), dimension(:,:,:)              , intent(in)     :: initial_data  
      integer                                , intent(in)     :: d1            
      integer                                , intent(in)     :: d2            
      integer                                , intent(in)     :: d3            
      type(cell_bc_wrapper_t), dimension(:), pointer, intent(in)     :: bc           
      type(boundary_parameters_t), pointer, intent(in) :: bc_params
      integer                                                 :: i

      call this%Init_quantity_init_arr (initial_data, d1, d2, d3, 1, bc_params)
      this%boundary_conditions => bc

   end subroutine

   subroutine Init_cell_quantity_init_val(this, initial_val, d1, d2, d3, bc, bc_params)
      implicit none
      class (cell_quantity_t)                , intent(in out) :: this          
      real(8)                                , intent(in)     :: initial_val   
      integer                                , intent(in)     :: d1            
      integer                                , intent(in)     :: d2            
      integer                                , intent(in)     :: d3            
      type(cell_bc_wrapper_t), dimension(:), pointer, intent(in) :: bc           
      type(boundary_parameters_t), pointer, intent(in) :: bc_params
      integer                                                 :: i

      call this%Init_quantity_init_val (initial_val, d1, d2, d3, 1, bc_params)
      this%boundary_conditions => bc
   end subroutine

   subroutine Init_cell_quantity_no_init(this, d1, d2, d3, bc, bc_params)
      implicit none
      class (cell_quantity_t)                , intent(in out) :: this          
      integer                                , intent(in)     :: d1            
      integer                                , intent(in)     :: d2            
      integer                                , intent(in)     :: d3            
      type(cell_bc_wrapper_t), dimension(:), pointer, intent(in) :: bc           
      type(boundary_parameters_t), pointer, intent(in) :: bc_params
      integer                                                 :: i

      call this%Init_quantity_no_init (d1, d2, d3, 1, bc_params)
      this%boundary_conditions => bc
   end subroutine

   subroutine Init_cell_quantity_copy(this, copy)
      implicit none
      class (cell_quantity_t), intent(in out) :: this  
      class (cell_quantity_t), intent(in)     :: copy  

   end subroutine

   subroutine Clean_cell_quantity (this)
      class (cell_quantity_t), intent(in out) :: this  

      deallocate (this%boundary_conditions)
      call this%Clean_quantity ()
   end subroutine Clean_cell_quantity


   subroutine Apply_boundary(this, is_blocking)
      implicit none
      class (cell_quantity_t)   , intent(in out) :: this
      logical, optional :: is_blocking 

      real(8), dimension (:,:,:), pointer        :: values
      integer                                    :: i, edge
      integer, dimension(:)     , pointer        :: ed_num
      logical :: is_blocking_local

      if (.not. present(is_blocking)) then
         is_blocking_local = .True.  
      else
         is_blocking_local = is_blocking   
      end if

      call this%Point_to_data (1, values)
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

   subroutine Write_cell_quantity(this, unit, iostat, iomsg)
      class (cell_quantity_t), intent(in) :: this  
      integer,      intent(in)     :: unit
      integer,      intent(out)    :: iostat
      character(*), intent(in out)  :: iomsg

#ifdef DEBUG
      write(*,*) '@@@ in Write_cell_quantity @@@'
#endif

      call this%Write_quantity(unit, iostat=iostat, iomsg=iomsg)


#ifdef DEBUG

      write(*,*) '@@@ end Write_cell_quantity @@@'
#endif

   end subroutine Write_cell_quantity

   subroutine Read_cell_quantity(this, unit, iostat, iomsg)
      class (cell_quantity_t), intent(in out) :: this  
      integer,      intent(in)     :: unit
      integer,      intent(out)    :: iostat
      character(*), intent(in out)  :: iomsg

#ifdef DEBUG
      write(*,*) '@@@ in Read_cell_quantity @@@'
#endif

      call this%Read_quantity(unit, iostat=iostat, iomsg=iomsg)


#ifdef DEBUG

      write(*,*) '@@@ end Read_cell_quantity @@@'
#endif

   end subroutine Read_cell_quantity

end module cell_quantity_module
