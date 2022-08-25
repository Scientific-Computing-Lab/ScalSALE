
module density_module
   use cell_quantity_module, only : cell_quantity_t
   use cell_boundary_condition_module, only : cell_boundary_condition_t, cell_bc_wrapper_t
   use materials_in_cells_module       , only : materials_in_cells_t
   use boundary_parameters_module      , only : boundary_parameters_t

   implicit none
   private
   public :: density_t

   type, extends(cell_quantity_t) :: density_t
      private


   contains


      procedure, public :: Clean_density

      procedure, public :: Write_density

      procedure, public :: Write_quantity_abstract => Write_density

      procedure, public :: Read_density
      procedure, public :: Read_quantity_abstract => Read_density

   end type


   interface density_t

      module procedure Constructor

      module procedure Constructor_materials_arr

      module procedure Copy_constructor
   end interface density_t

contains

   type(density_t) function Constructor(initial_val, d1, d2, d3, bc, bc_params)
      implicit none
      real(8)                                           , intent(in) :: initial_val 
      integer                                           , intent(in) :: d1          
      integer                                           , intent(in) :: d2          
      integer                                           , intent(in) :: d3          
      type(cell_bc_wrapper_t), dimension(:), pointer    , intent(in) :: bc           
      type(boundary_parameters_t), pointer, intent(in) :: bc_params


      call Constructor%Init_cell_quantity_init_val (initial_val, d1, d2, d3, bc, bc_params)
   end function



   type(density_t) function Constructor_materials_arr(init_density, mat_cells, d1, d2, d3, bc, bc_params)
      implicit none
      real(8), dimension (:), allocatable               , intent(in)    :: init_density 
      type(materials_in_cells_t)                        , intent(inout) :: mat_cells    
      integer                                           , intent(in)    :: d1           
      integer                                           , intent(in)    :: d2           
      integer                                           , intent(in)    :: d3           
      type(cell_bc_wrapper_t), dimension(:), pointer    , intent(in)    :: bc           
      real(8), dimension(:, :, :), allocatable                          :: init_values  
      integer                                                           :: i, j, k
      real(8), dimension (:, :, :), pointer                             :: mat_cell  
      type(boundary_parameters_t), pointer, intent(in) :: bc_params

      allocate(init_values(d1,d2,d3))
      call mat_cells%Point_to_data(mat_cell)

      init_values = 0d0
      if (d3 == 1) then
            do j = 1, d2 - 1
               do i = 1, d1 - 1
                  if (mat_cell(i, j, 1) /= 0) then
                     init_values(i, j ,1) = init_density(int(mat_cell(i, j, 1)))
                  end if
               end do
            end do
      else
         do k = 1, d3 - 1
            do j = 1, d2 - 1
               do i = 1, d1 - 1
                  if (mat_cell(i, j, k) /= 0) then
                     init_values(i, j ,k) = init_density(int(mat_cell(i, j, k)))
                  end if
               end do
            end do
         end do
      end if
      call Constructor_materials_arr%Init_cell_quantity_init_arr (init_values, d1, d2, d3, bc, bc_params)
      deallocate(init_values)
   end function

   type(density_t) function Copy_constructor(copy)
      type (density_t), intent(in) :: copy  


      call Copy_constructor%Init_cell_quantity_copy (copy)
   end function


   subroutine Clean_density (this)
      class (density_t), intent(in out) :: this 

      call this%Clean_cell_quantity ()
   end subroutine Clean_density

   subroutine Write_density(this, unit, iostat, iomsg)
      class (density_t), intent(in) :: this  
      integer,      intent(in)     :: unit
      integer,      intent(out)    :: iostat
      character(*), intent(in out)  :: iomsg

#ifdef DEBUG
      write(*,*) '@@@ in Write_density @@@'
#endif

      call this%Write_cell_quantity(unit, iostat=iostat, iomsg=iomsg)

#ifdef DEBUG
      write(*,*) '@@@ end Write_density @@@'
#endif

   end subroutine Write_density

   subroutine Read_density(this, unit, iostat, iomsg)
      class (density_t), intent(in out) :: this  
      integer,      intent(in)     :: unit
      integer,      intent(out)    :: iostat
      character(*), intent(in out)  :: iomsg

#ifdef DEBUG
      write(*,*) '@@@ in Read_density @@@'
#endif

      call this%Read_cell_quantity(unit, iostat=iostat, iomsg=iomsg)

#ifdef DEBUG
      write(*,*) '@@@ end Read_density @@@'
#endif

   end subroutine Read_density

end module density_module
