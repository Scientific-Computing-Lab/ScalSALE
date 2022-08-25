
module volume_module
   use cell_quantity_module, only : cell_quantity_t
   use cell_boundary_condition_module, only : cell_boundary_condition_t, cell_bc_wrapper_t
   use boundary_parameters_module      , only : boundary_parameters_t

   implicit none
   private
   public :: volume_t

   type, extends(cell_quantity_t) :: volume_t
      private


   contains


      procedure, public :: Clean_volume

      procedure, public :: Calculate

      procedure, public :: Write_volume
      procedure, public :: Write_quantity_abstract => Write_volume

      procedure, public :: Read_volume
      procedure, public :: Read_quantity_abstract => Read_volume

   end type


   interface volume_t

      module procedure Constructor
   end interface volume_t

contains

   type(volume_t) function Constructor(initial_val, d1, d2, d3, bc, bc_params)
      implicit none
      real(8)                                           , intent(in) :: initial_val 
      integer                                           , intent(in) :: d1          
      integer                                           , intent(in) :: d2          
      integer                                           , intent(in) :: d3          
      type(cell_bc_wrapper_t), dimension(:), pointer    , intent(in) :: bc           
      type(boundary_parameters_t), pointer, intent(in) :: bc_params

      call Constructor%Init_cell_quantity_init_val (initial_val, d1, d2, d3, bc, bc_params)
   end function

   subroutine Calculate(this, coordinates, sw_vertex_mass_optional, cyl_optional)
      use quantity_module, only: quantity_t
      use geometry_module, only: x_proj, Polygon_volume
      use general_utils_module, only : Get_axises_num
      implicit none
      class (volume_t)         , intent(in out) :: this   
      class (quantity_t)       , intent(in out) :: coordinates   
      real(8), optional        , intent(in)     :: cyl_optional
      integer, optional        , intent(in)     :: sw_vertex_mass_optional

      real(8), dimension(:,:,:), pointer  :: vol    
      real(8), dimension(:,:,:), pointer  :: x      
      real(8), dimension(:,:,:), pointer  :: y      
      real(8), dimension(:,:,:), pointer  :: z      
      real(8), dimension(4)                           :: x_vert 
      real(8), dimension(4)                           :: y_vert

      integer                             :: nx, ny, nz  

      integer :: i, j, k   
      integer :: is_neg    
      integer :: ip,jp,kp
      real(8) :: x1, x2, x3, x4, x5, x6, x7, x8  
      real(8) :: y1, y2, y3, y4, y5, y6, y7, y8  
      real(8) :: z1, z2, z3, z4, z5, z6, z7, z8  
      real(8) :: cyl
      integer :: sw_vertex_mass
      integer                                         :: dimension



      if (.not. present(cyl_optional)) then
         cyl = 0d0
      else
         cyl = cyl_optional
      end if

      if (.not. present(sw_vertex_mass_optional)) then
         sw_vertex_mass = 1
      else
         sw_vertex_mass = sw_vertex_mass_optional
      end if

      dimension = Get_axises_num (coordinates%d1, coordinates%d2, coordinates%d3)
      call this%Point_to_data(vol)
      nx = this%d1
      ny = this%d2

      if (dimension == 2) then
         call coordinates%Point_to_data(x, y)

         do j = 1, ny
                  do i = 1, nx
                     x_vert(1) = x(i    , j    , 1)
                     x_vert(2) = x(i + 1, j    , 1)
                     x_vert(3) = x(i + 1, j + 1, 1)
                     x_vert(4) = x(i    , j + 1, 1)
                     y_vert(1) = y(i    , j    , 1)
                     y_vert(2) = y(i + 1, j    , 1)
                     y_vert(3) = y(i + 1, j + 1, 1)
                     y_vert(4) = y(i    , j + 1, 1)
                     call Polygon_volume(4, x_vert, y_vert, vol(i, j, 1), cyl)
                  end do
               end do
   else
      call coordinates%Point_to_data(x, y, z)

     is_neg = 0
     nz = this%d3
         do k = 1, nz
            do j = 1, ny
               do i = 1, nx
                  ip = i + 1
                  jp = j + 1
                  kp = k + 1
                  x1 = x(i , j , k  )
                  x2 = x(ip, j , k  )
                  x3 = x(ip, jp, k  )
                  x4 = x(i , jp, k  )
                  x5 = x(i , j , kp )
                  x6 = x(ip, j , kp )
                  x7 = x(ip, jp, kp )
                  x8 = x(i , jp, kp )

                  y1 = y(i , j , k )
                  y2 = y(ip, j , k )
                  y3 = y(ip, jp, k )
                  y4 = y(i , jp, k )
                  y5 = y(i , j , kp)
                  y6 = y(ip, j , kp)
                  y7 = y(ip, jp, kp)
                  y8 = y(i , jp, kp)

                  z1 = z(i , j , k  )
                  z2 = z(ip, j , k  )
                  z3 = z(ip, jp, k  )
                  z4 = z(i , jp, k  )
                  z5 = z(i , j , kp)
                  z6 = z(ip, j , kp)
                  z7 = z(ip, jp, kp)
                  z8 = z(i , jp, kp)

                  vol(i, j, k) = (-x_proj(y2, z2, y3, z3, y4, z4, y5, z5, y6, z6, y8, z8) * x1 &
                                  -x_proj(y3, z3, y4, z4, y1, z1, y6, z6, y7, z7, y5, z5) * x2 &
                                  -x_proj(y4, z4, y1, z1, y2, z2, y7, z7, y8, z8, y6, z6) * x3 &
                                  -x_proj(y1, z1, y2, z2, y3, z3, y8, z8, y5, z5, y7, z7) * x4 &
                                  +x_proj(y6, z6, y7, z7, y8, z8, y1, z1, y2, z2, y4, z4) * x5 &
                                  +x_proj(y7, z7, y8, z8, y5, z5, y2, z2, y3, z3, y1, z1) * x6 &
                                  +x_proj(y8, z8, y5, z5, y6, z6, y3, z3, y4, z4, y2, z2) * x7 &
                                  +x_proj(y5, z5, y6, z6, y7, z7, y4, z4, y1, z1, y3, z3) * x8 ) / 12d0



               end do
            end do
         end do






      end if
   end subroutine Calculate

   subroutine Clean_volume (this)
      class (volume_t), intent(in out) :: this   

      call this%Clean_cell_quantity ()
   end subroutine Clean_volume

   subroutine Write_volume(this, unit, iostat, iomsg)
      class (volume_t), intent(in) :: this  
      integer,      intent(in)     :: unit
      integer,      intent(out)    :: iostat
      character(*), intent(in out)  :: iomsg

#ifdef DEBUG
      write(*,*) '@@@ in Write_volume @@@'
#endif

      call this%Write_cell_quantity(unit, iostat=iostat, iomsg=iomsg)

#ifdef DEBUG
      write(*,*) '@@@ end Write_volume @@@'
#endif

   end subroutine Write_volume

   subroutine Read_volume(this, unit, iostat, iomsg)
      class (volume_t), intent(in out) :: this  
      integer,      intent(in)     :: unit
      integer,      intent(out)    :: iostat
      character(*), intent(in out)  :: iomsg

#ifdef DEBUG
      write(*,*) '@@@ in Read_volume @@@'
#endif

      call this%Read_cell_quantity(unit, iostat=iostat, iomsg=iomsg)

#ifdef DEBUG
      write(*,*) '@@@ end Read_volume @@@'
#endif

   end subroutine Read_volume

end module volume_module
