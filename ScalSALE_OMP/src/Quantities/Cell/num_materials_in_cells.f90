
module num_materials_in_cells_module
   use cell_quantity_module, only : cell_quantity_t
   use cell_boundary_condition_module, only : cell_boundary_condition_t, cell_bc_wrapper_t
   use materials_in_cells_module       , only : materials_in_cells_t
   use boundary_parameters_module      , only : boundary_parameters_t

   implicit none
   private
   public :: num_materials_in_cells_t

   type, extends(cell_quantity_t) :: num_materials_in_cells_t
      private


   contains

      procedure, public :: Clean_num_materials_in_cells

      procedure, public :: Write_num_materials_in_cells
      procedure, public :: Write_quantity_abstract => Write_num_materials_in_cells

      procedure, public :: Read_num_materials_in_cells
      procedure, public :: Read_quantity_abstract => Read_num_materials_in_cells

   end type


   interface num_materials_in_cells_t

      module procedure Constructor
      module procedure Constructor_mat_cells

   end interface num_materials_in_cells_t

contains

   type(num_materials_in_cells_t) function Constructor(initial_val, d1, d2, d3, bc, bc_params)
      real(8)                                           , intent(in) :: initial_val  
      integer                                           , intent(in) :: d1           
      integer                                           , intent(in) :: d2           
      integer                                           , intent(in) :: d3           
      type(cell_bc_wrapper_t), dimension(:), pointer    , intent(in) :: bc           
      type(boundary_parameters_t), pointer, intent(in) :: bc_params


      call Constructor%Init_cell_quantity_init_val (initial_val, d1, d2, d3, bc, bc_params)
   end function

   type(num_materials_in_cells_t) function Constructor_mat_cells(mat_cells, d1, d2, d3, bc, bc_params)
      type(materials_in_cells_t), pointer               , intent(in) :: mat_cells  
      integer                                           , intent(in) :: d1           
      integer                                           , intent(in) :: d2           
      integer                                           , intent(in) :: d3           
      type(cell_bc_wrapper_t), dimension(:), pointer    , intent(in) :: bc           
      type(boundary_parameters_t), pointer, intent(in) :: bc_params

      real(8), dimension (:, :, :), pointer                             ::  mat_cell  
      real(8), dimension (d1, d2, d3)                                   ::  init_values  
      integer :: i, j, k

      call mat_cells%Point_to_data (mat_cell)
      init_values = 0d0
      do k = 1, d3
         do j = 1, d2
            do i = 1, d1
               if (mat_cell(i, j, k) /= 0) then
                  init_values(i, j, k) = init_values(i, j, k) + 1
               end if
            end do
         end do
      end do

      call Constructor_mat_cells%Init_cell_quantity_init_arr (init_values, d1, d2, d3, bc, bc_params)
   end function


   subroutine Clean_num_materials_in_cells (this)
      class (num_materials_in_cells_t), intent(in out) :: this 

      call this%Clean_cell_quantity ()
   end subroutine Clean_num_materials_in_cells

   subroutine Write_num_materials_in_cells(this, unit, iostat, iomsg)
      class (num_materials_in_cells_t), intent(in) :: this  
      integer,      intent(in)     :: unit
      integer,      intent(out)    :: iostat
      character(*), intent(in out)  :: iomsg

#ifdef DEBUG
      write(*,*) '@@@ in Write_num_materials_in_cells @@@'
#endif

      call this%Write_cell_quantity(unit, iostat=iostat, iomsg=iomsg)

#ifdef DEBUG
      write(*,*) '@@@ end Write_num_materials_in_cells @@@'
#endif

   end subroutine Write_num_materials_in_cells

   subroutine Read_num_materials_in_cells(this, unit, iostat, iomsg)
      class (num_materials_in_cells_t), intent(in out) :: this  
      integer,      intent(in)     :: unit
      integer,      intent(out)    :: iostat
      character(*), intent(in out)  :: iomsg

#ifdef DEBUG
      write(*,*) '@@@ in Read_num_materials_in_cells @@@'
#endif

      call this%Read_cell_quantity(unit, iostat=iostat, iomsg=iomsg)

#ifdef DEBUG
      write(*,*) '@@@ end Read_num_materials_in_cells @@@'
#endif

   end subroutine Read_num_materials_in_cells

end module num_materials_in_cells_module

