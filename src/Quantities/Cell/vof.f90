
module vof_module
   use cell_quantity_module, only : cell_quantity_t
   use cell_boundary_condition_module, only : cell_boundary_condition_t, cell_bc_wrapper_t
   use materials_in_cells_module       , only : materials_in_cells_t
   use boundary_parameters_module      , only : boundary_parameters_t

   implicit none
   private
   public :: vof_t

   type, extends(cell_quantity_t) :: vof_t
      private


   contains

      procedure, public :: Clean_vof

      procedure, public :: Write_vof
      procedure, public :: Write_quantity_abstract => Write_vof

      procedure, public :: Read_vof
      procedure, public :: Read_quantity_abstract => Read_vof

   end type


   interface vof_t

      module procedure Constructor_mat

      module procedure Constructor_total

      module procedure Constructor_reg

   end interface vof_t

contains

   type(vof_t) function Constructor_mat(d1, d2, d3, bc, bc_params, mat_cells, mat_index)
      implicit none
      integer                                           , intent(in)    :: d1           
      integer                                           , intent(in)    :: d2           
      integer                                           , intent(in)    :: d3           
      type(cell_bc_wrapper_t), dimension(:), pointer    , intent(in):: bc           
      type(materials_in_cells_t)                        , intent(inout) :: mat_cells
      integer                                           , intent(in)    :: mat_index
      type(boundary_parameters_t), pointer, intent(in) :: bc_params

      real(8), dimension (:, :, :), allocatable                         :: init_values  
      real(8), dimension (:, :, :), pointer                             ::  mat_cell  
      integer                                                           :: i, j, k
      allocate(init_values(d1,d2,d3))
      init_values = 0d0
      write(*,*) "done init vof"
      call mat_cells%Point_to_data(mat_cell)
      do k = 1, d3
         do j = 1, d2
            do i = 1, d1
               if (mat_cell(i, j, k) == mat_index .and. mat_index /= 0 ) then
                  init_values(i, j, k) = 1d0
               end if
            end do
         end do
      end do

      call Constructor_mat%Init_cell_quantity_init_arr (init_values, d1, d2, d3, bc, bc_params)
      deallocate(init_values)
   end function

   type(vof_t) function Constructor_reg(val, d1, d2, d3, bc, bc_params)
      implicit none
      integer                                           , intent(in)    :: d1           
      integer                                           , intent(in)    :: d2           
      integer                                           , intent(in)    :: d3           
      type(cell_bc_wrapper_t), dimension(:), pointer    , intent(in):: bc           
      real(8)                                           , intent(in)    :: val
      type(boundary_parameters_t), pointer, intent(in) :: bc_params


      call Constructor_reg%Init_cell_quantity_init_val (val, d1, d2, d3, bc, bc_params)

   end function


   type(vof_t) function Constructor_total(d1, d2, d3, bc, bc_params, mat_cells)
      integer                                           , intent(in)    :: d1           
      integer                                           , intent(in)    :: d2           
      integer                                           , intent(in)    :: d3           
      type(cell_bc_wrapper_t), dimension(:), pointer    , intent(in)    :: bc           
      type(materials_in_cells_t)                        , intent(inout) :: mat_cells
      type(boundary_parameters_t), pointer, intent(in) :: bc_params

      real(8), dimension (:, :, :), allocatable                         :: init_values  
      real(8), dimension (:, :, :), pointer                             ::  mat_cell  

      integer                                              :: i, j, k

      allocate(init_values(d1,d2,d3))

      init_values = 0d0

      call mat_cells%Point_to_data(mat_cell)

      do k = 1, d3
         do j = 1, d2
            do i = 1, d1
               if (mat_cell(i, j, k) /= 0) then
                  init_values(i, j, k) = 1d0
               end if
            end do
         end do
      end do

      call Constructor_total%Init_cell_quantity_init_arr (init_values, d1, d2, d3, bc, bc_params)
   end function


   subroutine Clean_vof (this)
      class (vof_t), intent(in out) :: this 

      call this%Clean_cell_quantity ()
   end subroutine Clean_vof

   subroutine Write_vof(this, unit, iostat, iomsg)
      class (vof_t), intent(in) :: this  
      integer,      intent(in)     :: unit
      integer,      intent(out)    :: iostat
      character(*), intent(in out)  :: iomsg

#ifdef DEBUG
      write(*,*) '@@@ in Write_vof @@@'
#endif

      call this%Write_cell_quantity(unit, iostat=iostat, iomsg=iomsg)

#ifdef DEBUG
      write(*,*) '@@@ end Write_vof @@@'
#endif

   end subroutine Write_vof

   subroutine Read_vof(this, unit, iostat, iomsg)
      class (vof_t), intent(in out) :: this  
      integer,      intent(in)     :: unit
      integer,      intent(out)    :: iostat
      character(*), intent(in out)  :: iomsg

#ifdef DEBUG
      write(*,*) '@@@ in Read_vof @@@'
#endif

      call this%Read_cell_quantity(unit, iostat=iostat, iomsg=iomsg)

#ifdef DEBUG
      write(*,*) '@@@ end Read_vof @@@'
#endif

   end subroutine Read_vof

end module vof_module

