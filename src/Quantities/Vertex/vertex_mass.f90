
module vertex_mass_module
   use vertex_quantity_module          , only : vertex_quantity_t
   use vertex_boundary_condition_module, only : vertex_boundary_condition_t, vertex_bc_wrapper_t
   use geometry_module                 , only : Quad_area, Tetrahederon_volume
   use cell_mass_module                , only : cell_mass_t
   use density_module                  , only : density_t
   use coordinates_module                  , only : coordinates_t
   use parallel_parameters_module, only : parallel_parameters_t
   use boundary_parameters_module      , only : boundary_parameters_t

   implicit none
   private
   public :: vertex_mass_t

   type, extends(vertex_quantity_t) :: vertex_mass_t
      private


   contains


      procedure, public :: Clean_vertex_mass
      procedure, public :: Calculate_vertex_mass_2d
      procedure, public :: Calculate_vertex_mass_3d
   end type


   interface vertex_mass_t

      module procedure Constructor

   end interface vertex_mass_t

contains



   type(vertex_mass_t) function Constructor(init_data, nxp, nyp, nzp, bc, bc_params)
      real(8)                                           , intent(in) :: init_data
      integer                                           , intent(in) :: nxp           
      integer                                           , intent(in) :: nyp           
      integer                                           , intent(in) :: nzp           
      type (vertex_bc_wrapper_t) , dimension(:), pointer, intent(in) :: bc           
      type(boundary_parameters_t), pointer, intent(in) :: bc_params

      call Constructor%Init_vertex_quantity_init_val (init_data, nxp, nyp, nzp, 1, bc, bc_params)

   end function

   subroutine Calculate_vertex_mass_2d(this, coordinates, tot_density, c_mass, wilkins_scheme, cyl)
      use constants_module, only : THIRD
      implicit none
      class (vertex_mass_t), intent(in out) :: this
      type(cell_mass_t)    , intent(in)     :: c_mass
      type(density_t)      , intent(in out) :: tot_density
      type(coordinates_t)      , intent(in out) :: coordinates

      integer, intent(in) :: wilkins_scheme 
      real(8), intent(in) :: cyl

      real(8), dimension(:, :, :), pointer :: cell_mass  
      real(8), dimension(:, :, :), pointer :: x          
      real(8), dimension(:, :, :), pointer :: y          
      real(8), dimension(:, :, :), pointer :: density    

      real(8), dimension(:, :, :), pointer :: vertex_mass

      integer, dimension(4) :: i_2_add = (/1,  0, -1,  0/)  
      integer, dimension(4) :: j_2_add = (/0,  1,  0, -1/)  
      integer, dimension(4) :: i_3_add = (/0, -1,  0,  1/)  
      integer, dimension(4) :: j_3_add = (/1,  0, -1,  0/)  
      integer, dimension(4) :: i_add   = (/0, -1, -1,  0/)  
      integer, dimension(4) :: j_add   = (/0,  0, -1, -1/)  

      real(8) :: x1 
      real(8) :: y1 
      real(8) :: x2 
      real(8) :: y2 
      real(8) :: x3 
      real(8) :: y3 
      real(8) :: r_factor_1 
      real(8) :: r_factor_2 
      real(8) :: r_factor_3 
      real(8) :: triag_volume 
      real(8) :: omcyl

      integer :: i, j
      integer :: k          
      integer :: nxp,nyp

      nxp = this%d1
      nyp = this%d2
      omcyl = 1d0 - cyl

      call coordinates%Point_to_data (x, y)
      call this%Point_to_data (1, vertex_mass)
      call tot_density%Point_to_data(density)
      call c_mass%Point_to_data (cell_mass)

         vertex_mass(1:nxp, 1:nyp, 1) = 0d0


      do j = 1, nyp
         do i = 1, nxp
            x1 = x(i, j, 1)
            y1 = y(i, j, 1)
            r_factor_1 = x1 * cyl + omcyl
            do k = 1, 4
               if (cell_mass(i + i_add(k), j + j_add(k), 1) > 0d0) then
                  x2 = x(i + i_2_add(k), j + j_2_add(k), 1)
                  y2 = y(i + i_2_add(k), j + j_2_add(k), 1)
                  r_factor_2 = x2 * cyl + omcyl
                  x3 = x(i + i_3_add(k), j + j_3_add(k), 1)
                  y3 = y(i + i_3_add(k), j + j_3_add(k), 1)
                  r_factor_3 = x3 * cyl + omcyl
                  triag_volume = (r_factor_1 + r_factor_2 + r_factor_3) * ((x2 - x1) * (y3 - y1) - (y2 - y1) * (x3 - x1))* &
                         THIRD * 0.25d0
                  vertex_mass(i, j, 1) = vertex_mass(i, j, 1) + triag_volume * density(i + i_add(k), j + j_add(k), 1)
               end if
            end do
         end do
      end do

   end subroutine Calculate_vertex_mass_2d



   subroutine Calculate_vertex_mass_3d(this, coordinates, tot_density, c_mass)
      use geometry_module, only : Tetrahederon_volume
      implicit none
      class(vertex_mass_t), intent(in out) :: this
      type(cell_mass_t)   , intent(in)     :: c_mass
      type(density_t)     , intent(in)     :: tot_density
      class(coordinates_t)  , intent(in out) :: coordinates

      real(8), dimension(:, :, :), pointer :: vertex_mass  
      real(8), dimension(:, :, :), pointer :: cell_mass    
      real(8), dimension(:, :, :), pointer :: x            
      real(8), dimension(:, :, :), pointer :: y            
      real(8), dimension(:, :, :), pointer :: z            
      real(8), dimension(:, :, :), pointer :: density      
      integer, save :: cyc = 0
      integer, dimension(:), pointer :: i_virt, j_virt, k_virt
      integer, dimension(8, 3) :: i_tetr = reshape([-1, 1, 0, 1, 0, 0, 0, 1 , &
                                                     0, 0,-1, 0, 0, 0, 0, 0 , &
                                                     0, 0, 0, 0,-1, 1,-1, 0 ], [8, 3])  
      integer, dimension(8, 3) :: j_tetr = reshape([ 0, 0, 1, 0, 0,-1, 1, 0, &
                                                    -1, 0, 0, 1,-1, 0, 0, 0, &
                                                     0,-1, 0, 0, 0, 0, 0, 1 ], [8, 3]) 
      integer, dimension(8, 3) :: k_tetr = reshape([ 0, 0, 0, 0, 1, 0, 0, 0, &
                                                     0,-1, 0, 0, 0, 1, 1, 1, &
                                                    -1, 0,-1,-1, 0, 0, 0, 0 ], [8, 3])  

      integer :: i, j, k        
      integer :: t              
      integer :: ii, jj, kk     
      integer :: i1, i2, i3  
      integer :: j1, j2, j3  
      integer :: virt_nxp, virt_nyp, virt_nzp
      integer :: k1, k2, k3  
      logical :: wall_x_top, wall_x_bot, wall_y_top, wall_y_bot, wall_z_top, wall_z_bot
      integer :: nxp,nyp,nzp

      nxp = this%d1
      nyp = this%d2
      nzp = this%d3

      wall_x_top = this%parallel_params%is_wall_x_top
      wall_x_bot = this%parallel_params%is_wall_x_bot
      wall_y_top = this%parallel_params%is_wall_y_top
      wall_y_bot = this%parallel_params%is_wall_y_bot
      wall_z_top = this%parallel_params%is_wall_z_top
      wall_z_bot = this%parallel_params%is_wall_z_bot
      virt_nxp = this%parallel_params%virt_nxp
      virt_nyp = this%parallel_params%virt_nyp
      virt_nzp = this%parallel_params%virt_nzp
      call this%parallel_params%Point_to_virtual_array(i_virt, j_virt, k_virt)


      call this       %Point_to_data(1, vertex_mass)
      call c_mass     %Point_to_data(cell_mass)
      call tot_density%Point_to_data(density)
      call coordinates%Point_to_data(x, y, z)

         vertex_mass(1:nxp, 1:nyp , 1:nzp) = 0d0
         do k = 1, nzp
            do j = 1, nyp
               do i = 1, nxp
                  do t = 1, 8
                     kk = k - (8 - t) / 4
                     jj = j - mod(t + 1, 4) / 2
                     ii = i - mod(t, 2)
                     if (cell_mass(ii, jj, kk) /= 0d0) then
                        i1 = i + i_tetr(t, 1)
                        i2 = i + i_tetr(t, 2)
                        i3 = i + i_tetr(t, 3)
                        j1 = j + j_tetr(t, 1)
                        j2 = j + j_tetr(t, 2)
                        j3 = j + j_tetr(t, 3)
                        k1 = k + k_tetr(t, 1)
                        k2 = k + k_tetr(t, 2)
                        k3 = k + k_tetr(t, 3)

                        if ( .false. .eqv. ( (i_virt(ii) == 0        .and. wall_x_bot .eqv. .true.) .or. &
                             (i_virt(ii) == virt_nxp .and. wall_x_top .eqv. .true.) .or. &
                             (j_virt(jj) == 0        .and. wall_y_bot .eqv. .true.) .or. &
                             (j_virt(jj) == virt_nyp .and. wall_y_top .eqv. .true.) .or. &
                             (k_virt(kk) == 0        .and. wall_z_bot .eqv. .true.) .or. &
                             (k_virt(kk) == virt_nzp .and. wall_z_top .eqv. .true.) )       ) then
                        	vertex_mass(i, j, k) = vertex_mass(i, j, k) + density(ii, jj, kk) * &
                            	                   Tetrahederon_volume(x(i  , j  , k  ), y(i  , j  , k  ), z(i  , j  , k  ), &
                            	                                       x(i1, j1, k1), y(i1, j1, k1), z(i1, j1, k1), &
                            	                                       x(i2, j2, k2), y(i2, j2, k2), z(i2, j2, k2), &
                            	                                       x(i3, j3, k3), y(i3, j3, k3), z(i3, j3, k3))
						     else 
								cycle
                             end if
                     end if
                  end do
                  vertex_mass(i, j, k) = vertex_mass(i, j, k) / 8d0
               end do
            end do
         end do






cyc = cyc + 1
   end subroutine Calculate_vertex_mass_3d


   subroutine Clean_vertex_mass (this)
      class (vertex_mass_t), intent(in out) :: this 

      call this%Clean_vertex_quantity ()
   end subroutine Clean_vertex_mass


end module vertex_mass_module
