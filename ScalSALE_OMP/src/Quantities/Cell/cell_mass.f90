
module cell_mass_module
   use cell_quantity_module, only : cell_quantity_t
   use cell_boundary_condition_module, only : cell_boundary_condition_t, cell_bc_wrapper_t
   use volume_module, only : volume_t
   use density_module, only : density_t
   use boundary_parameters_module      , only : boundary_parameters_t

   implicit none
   private
   public :: cell_mass_t

   type, extends(cell_quantity_t) :: cell_mass_t
      private


   contains


      procedure, public :: Clean_cell_mass
      procedure, public :: Calculate_cell_mass

      procedure, public :: Write_cell_mass
      procedure, public :: Write_quantity_abstract => Write_cell_mass

      procedure, public :: Read_cell_mass
      procedure, public :: Read_quantity_abstract => Read_cell_mass

   end type


   interface cell_mass_t

      module procedure Constructor_init_val

      module procedure Constructor_vol_density

      module procedure Copy_constructor
   end interface cell_mass_t

contains


   type(cell_mass_t) function Constructor_init_val(initial_val, d1, d2, d3, bc, bc_params)
      implicit none
      real(8)                                           , intent(in) :: initial_val 
      integer                                           , intent(in) :: d1            
      integer                                           , intent(in) :: d2            
      integer                                           , intent(in) :: d3            
      type(cell_bc_wrapper_t), dimension(:), pointer    , intent(in) :: bc           
      type(boundary_parameters_t), pointer, intent(in) :: bc_params


      call Constructor_init_val%Init_cell_quantity_init_val (initial_val, d1, d2, d3, bc, bc_params)
   end function


   type(cell_mass_t) function Constructor_vol_density(total_volume, total_density, d1, d2, d3, bc, bc_params)
      implicit none
      type (volume_t)                                   , intent(in out) :: total_volume  
      type (density_t)                                  , intent(in out) :: total_density 
      integer                                           , intent(in)     :: d1            
      integer                                           , intent(in)     :: d2            
      integer                                           , intent(in)     :: d3            
      type(cell_bc_wrapper_t), dimension(:), pointer    , intent(in)     :: bc           
      type(boundary_parameters_t), pointer, intent(in) :: bc_params

      real(8), dimension (:, :, :), allocatable                          :: init_values   
      real(8), dimension (:,:,:), pointer                                :: density       
      real(8), dimension (:,:,:), pointer                                :: vol           
      integer                                                            :: i, j , k

      allocate(init_values(d1,d2,d3))
      call total_volume%Point_to_data(1, vol)
      call total_density%Point_to_data(1, density)
      do i = 1, d1
         do j = 1, d2
            do k = 1, d3
               init_values(i, j ,k) = vol(i, j ,k) * density(i, j, k)
            end do
         end do
      end do

      call Constructor_vol_density%Init_cell_quantity_init_arr (init_values, d1, d2, d3, bc, bc_params)
      deallocate(init_values)
   end function

   type(cell_mass_t) function Copy_constructor(copy)
      type (cell_mass_t), intent(in) :: copy 

      call Copy_constructor%Init_cell_quantity_copy (copy)
   end function


   subroutine Clean_cell_mass (this)
      class (cell_mass_t), intent(in out) :: this 

      call this%Clean_cell_quantity ()
   end subroutine Clean_cell_mass

   subroutine Calculate_cell_mass(this, vol, density)
      implicit none
      class (cell_mass_t)                               , intent(in out)  :: this
      type (volume_t)                                   , intent(in out)  :: vol  
      type (density_t)                                  , intent(in out)  :: density 

      integer                                                           :: i, j , k
      real(8), dimension(:, :, :), pointer                              :: rho, volume, cm
      integer :: nxp, nzp, nyp
      nxp = this%d1 + 1
      nyp = this%d2 + 1
      if (this%d3 == 1) then
         nzp = this%d3
      else
         nzp = this%d3 + 1
      end if

      call density%Point_to_data(rho)
      call vol%Point_to_data(volume)
call this%Point_to_data(cm)

      do k = 1, nzp
         do j = 1, nyp
            do i = 1, nxp
               cm(i, j ,k) = volume(i, j ,k) * rho(i, j, k)
            end do
         end do
      end do

   end subroutine Calculate_cell_mass

   subroutine Write_cell_mass(this, unit, iostat, iomsg)
      class (cell_mass_t), intent(in) :: this  
      integer,      intent(in)     :: unit
      integer,      intent(out)    :: iostat
      character(*), intent(in out)  :: iomsg

#ifdef DEBUG
      write(*,*) '@@@ in Write_cell_mass @@@'
#endif

      call this%Write_cell_quantity(unit, iostat=iostat, iomsg=iomsg)

#ifdef DEBUG
      write(*,*) '@@@ end Write_cell_mass @@@'
#endif

   end subroutine Write_cell_mass

   subroutine Read_cell_mass(this, unit, iostat, iomsg)
      class (cell_mass_t), intent(in out) :: this  
      integer,      intent(in)     :: unit
      integer,      intent(out)    :: iostat
      character(*), intent(in out)  :: iomsg

#ifdef DEBUG
      write(*,*) '@@@ in Read_cell_mass @@@'
#endif
      call this%Read_cell_quantity(unit, iostat=iostat, iomsg=iomsg)

#ifdef DEBUG
      write(*,*) '@@@ end Read_cell_mass @@@'
#endif

   end subroutine Read_cell_mass

end module cell_mass_module
