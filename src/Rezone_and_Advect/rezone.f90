
module rezone_module
   use data_module     , only : data_t
   use coordinates_module     , only : coordinates_t
   use velocity_module , only : velocity_t
   use mesh_base_module, only : mesh_base_t!, mesh_t, mesh_3d_t
use mesh_module, only : mesh_t
use mesh_3d_module, only : mesh_3d_t
   use volume_module   , only : volume_t
   use vertex_mass_module              , only : vertex_mass_t
   use vertex_boundary_condition_module, only : vertex_boundary_condition_t, vertex_bc_wrapper_t
   use cell_boundary_condition_module  , only : cell_boundary_condition_t  , cell_bc_wrapper_t
   use boundary_parameters_module, only : boundary_parameters_t

   implicit none
   private
   public :: rezone_t

   type :: rezone_t
      private


      class (mesh_base_t) , pointer, public :: mesh        
      type (velocity_t)   , pointer, public :: velocity    
      type (vertex_mass_t), pointer, public :: vertex_mass 

      type (coordinates_t), pointer, public :: material_coordinates

      type(velocity_t)          , pointer, public :: mesh_velocity    
      type (volume_t)            , pointer, public :: material_volume  

      real(8), public :: velocity_limit 
      real(8), public :: mass_threshold 

      integer, public :: rezone_type    
      integer :: nyp, nxp, nzp  
      integer :: dimension      
      integer :: cyc_delete      

   contains


      procedure, public :: Clean_rezone

      procedure, public :: Calculate_rezone_2d

      procedure, public :: Calculate_rezone_3d

      procedure, public :: Set_communication

      procedure, public :: Point_to_coordinates_2d

      procedure, public :: Point_to_coordinates_3d

      generic,   public    :: Point_to_data  =>     &
                              Point_to_coordinates_2d,   &
                              Point_to_coordinates_3d


      procedure, public :: Point_to_velocities

      procedure, public :: Point_to_volume

   end type


   interface rezone_t
      module procedure Constructor
   end interface rezone_t

contains


   type(rezone_t) function Constructor(rezone_type, nxp, nyp, nzp, bc_c, bc_v, mesh_p, velocity_p, vertex_mass_p, bc_params)
      integer                                           , intent(in) :: rezone_type
      integer                                           , intent(in) :: nxp       
      integer                                           , intent(in) :: nyp       
      integer                                           , intent(in) :: nzp       
      type (vertex_bc_wrapper_t) , dimension(:), pointer, intent(inout) :: bc_v    
      type (cell_bc_wrapper_t) , dimension(:), pointer, intent(inout) :: bc_c    
      class (mesh_base_t)       , pointer, intent(inout) :: mesh_p        
      type (velocity_t)   , pointer, intent(in) :: velocity_p    
      type (vertex_mass_t), pointer, intent(in) :: vertex_mass_p 
      type(boundary_parameters_t), pointer, intent(inout) :: bc_params

      type(data_t)        :: xm, ym, zm

      Constructor%rezone_type = rezone_type

      Constructor%dimension = 3
      Constructor%nxp       = nxp
      Constructor%nyp       = nyp
      Constructor%nzp       = nzp

      Constructor%velocity_limit = 4d20
      Constructor%cyc_delete = 1
      Constructor%mass_threshold = 1d20

      allocate (Constructor%material_coordinates)
      allocate(Constructor%mesh_velocity)
      if (nzp == 1) then
         Constructor%mesh_velocity = velocity_t(0d0, nxp+1, nyp+1, 1, 2, bc_v, bc_params)
         Constructor%material_coordinates = coordinates_t(0d0, nxp + 1, nyp + 1, 1, mesh_p%coordinates%boundary_conditions&
                                                          , bc_params)
      else
         Constructor%mesh_velocity = velocity_t(0d0, nxp+1, nyp+1, nzp+1, 3, bc_v, bc_params)
         Constructor%material_coordinates = coordinates_t(0d0, nxp + 1, nyp + 1, nzp + 1, mesh_p%coordinates%boundary_conditions&
                                                          , bc_params)
      end if
      allocate(Constructor%material_volume)
      Constructor%material_volume = volume_t(0d0, nxp, nyp, nzp, bc_c, bc_params)

      Constructor%mesh => mesh_p

      Constructor%velocity => velocity_p
      Constructor%vertex_mass => vertex_mass_p

   end function

   subroutine Calculate_rezone_2d(this, dt)
      use geometry_module, only : Mirror_image
      implicit none
      class (rezone_t)   , intent(in out) :: this

      real(8), intent(in) :: dt

      real(8), dimension(:, :, :), pointer :: x               
      real(8), dimension(:, :, :), pointer :: y               
      real(8), dimension(:, :, :), pointer :: material_x      
      real(8), dimension(:, :, :), pointer :: material_y      
      real(8), dimension(:, :, :), pointer :: velocity_x      
      real(8), dimension(:, :, :), pointer :: velocity_y      
      real(8), dimension(:, :, :), pointer :: mesh_velocity_x 
      real(8), dimension(:, :, :), pointer :: mesh_velocity_y 
      real(8), dimension(:, :, :), pointer :: vertex_mass     
      real(8), dimension(:, :, :), pointer :: mat_vol     


      real(8) :: velocity_sq 
      real(8) :: factor      

      integer :: i, i1, i2, j, j1, j2




      call this%mesh%Point_to_data(x, y)
      call this%velocity%Point_to_data(velocity_x, velocity_y)
      call this%material_coordinates%Point_to_data (material_x, material_y)
      call this%mesh_velocity%Point_to_data(mesh_velocity_x, mesh_velocity_y)
      call this%vertex_mass%Point_to_data(vertex_mass)
      call this%material_volume%Point_to_data(mat_vol)




      material_x = x + velocity_x * dt
      material_y = y + velocity_y * dt

      call this%material_coordinates%Calculate(dt, this%velocity, coords=this%mesh%coordinates)
      call this%material_coordinates%Apply_boundary(this%material_coordinates%data)

      mesh_velocity_x = velocity_x
      mesh_velocity_y = velocity_y











      do j = 1, this%nyp
         do i = 1, this%nxp
            velocity_sq = velocity_y(i, j, 1) * velocity_y(i, j, 1) + velocity_x(i, j, 1) * velocity_x(i, j, 1)
            if ((vertex_mass(i, j, 1) < this%mass_threshold) .and. (velocity_sq > this%velocity_limit)) then
               factor = sqrt(this%velocity_limit / velocity_sq)
               velocity_x(i, j, 1) = velocity_x(i, j, 1) * factor
               velocity_y(i, j, 1) = velocity_y(i, j, 1) * factor
            end if
         end do
      end do

      select case (this%rezone_type)
         case (0)


            mesh_velocity_x = velocity_x
            mesh_velocity_y = velocity_y

         case (1)
            mesh_velocity_x = 0d0
            mesh_velocity_y = 0d0



         case default

      end select

      call this%mesh_velocity%Apply_boundary(this%mesh%coordinates%data)

      call this%material_volume%Calculate(this%material_coordinates, cyl_optional=this%mesh%cyl)

      return
   end subroutine Calculate_rezone_2d


   subroutine Calculate_rezone_3d(this, dt)
      use omp_lib
      implicit none
      class (rezone_t)   , intent(in out) :: this

      real(8), intent(in) :: dt

      real(8), dimension(:, :, :), pointer :: velocity_x      
      real(8), dimension(:, :, :), pointer :: velocity_y      
      real(8), dimension(:, :, :), pointer :: velocity_z      
      real(8), dimension(:, :, :), pointer :: x      
      real(8), dimension(:, :, :), pointer :: y      
      real(8), dimension(:, :, :), pointer :: z      
      real(8), dimension(:, :, :), pointer :: mesh_velocity_x 
      real(8), dimension(:, :, :), pointer :: mesh_velocity_y 
      real(8), dimension(:, :, :), pointer :: mesh_velocity_z 
      real(8), dimension(:, :, :), pointer :: material_x 
      real(8), dimension(:, :, :), pointer :: material_y 
      real(8), dimension(:, :, :), pointer :: material_z 
      real(8), dimension(:, :, :), pointer :: vertex_mass     
      real(8), dimension(:, :, :), pointer :: mat_vol          

      real(8) :: velocity_sq 
      real(8) :: factor      
      integer :: i, j, k

      call this%velocity%Point_to_data(velocity_x, velocity_y, velocity_z)
      call this%mesh%Point_to_data(x, y, z)
      call this%material_coordinates%Point_to_data (material_x, material_y, material_z)
      call this%material_volume%Point_to_data(mat_vol)
      call this%mesh_velocity%Point_to_data(mesh_velocity_x, mesh_velocity_y, mesh_velocity_z)
      call this%vertex_mass%Point_to_data(vertex_mass)








      do k = 0, this%nzp + 1
         do j = 0, this%nyp + 1
            do i = 0, this%nxp + 1
               material_x(i, j, k) = x(i, j, k) + velocity_x(i, j, k) * dt
               material_y(i, j, k) = y(i, j, k) + velocity_y(i, j, k) * dt
               material_z(i, j, k) = z(i, j, k) + velocity_z(i, j, k) * dt
            end do
         end do
      end do

      call this%material_volume%Calculate(this%material_coordinates)
      call this%material_volume%Exchange_virtual_space_nonblocking()






      do k = 1, this%nzp
         do j = 1, this%nyp
            do i = 1, this%nxp
               velocity_sq = velocity_z(i, j, k) * velocity_z(i, j, k) &
                  + velocity_y(i, j, k) * velocity_y(i, j, k) &
                  + velocity_x(i, j, k) * velocity_x(i, j, k)

               if ((vertex_mass(i, j, k) < this%mass_threshold) .and. (velocity_sq > this%velocity_limit)) then
                  factor = sqrt(this%velocity_limit / velocity_sq)
                  velocity_x(i, j, k) = velocity_x(i, j, k) * factor
                  velocity_y(i, j, k) = velocity_y(i, j, k) * factor
                  velocity_z(i, j, k) = velocity_z(i, j, k) * factor
               end if
            end do
         end do
      end do

      select case (this%rezone_type)
         case(0)


            do k = 1, this%nzp
               do j = 1, this%nyp
                  do i = 1, this%nxp
                     mesh_velocity_x(i, j, k) = velocity_x(i, j, k)
                     mesh_velocity_y(i, j, k) = velocity_y(i, j, k)
                     mesh_velocity_z(i, j, k) = velocity_z(i, j, k)
                  end do
               end do
            end do

         case(1)
            do k = 1, this%nzp
               do j = 1, this%nyp
                  do i = 1, this%nxp
                     mesh_velocity_x(i, j, k) = 0
                     mesh_velocity_y(i, j, k) = 0
                     mesh_velocity_z(i, j, k) = 0
                  end do
               end do
            end do

         case default

      end select

      call this%mesh_velocity%Apply_boundary(this%mesh%coordinates%data)
      call this%material_volume%Exchange_end()

      this%cyc_delete = this%cyc_delete + 1

      return
   end subroutine Calculate_rezone_3d

   subroutine Point_to_coordinates_2d (this, ptr_x, ptr_y)
      implicit none
      class (rezone_t)                         , intent(in out) :: this    
      real(8)       , dimension(:,:,:), pointer, intent(out)    :: ptr_x   
      real(8)       , dimension(:,:,:), pointer, intent(out)    :: ptr_y   

      call this%material_coordinates%data(1)%Point_to_data (ptr_x)
      call this%material_coordinates%data(2)%Point_to_data (ptr_y)

   end subroutine Point_to_coordinates_2d


   subroutine Point_to_coordinates_3d (this, ptr_x, ptr_y, ptr_z)
      implicit none
      class (rezone_t)                         , intent(in out) :: this    
      real(8)       , dimension(:,:,:), pointer, intent(out)    :: ptr_x   
      real(8)       , dimension(:,:,:), pointer, intent(out)    :: ptr_y   
      real(8)       , dimension(:,:,:), pointer, intent(out)    :: ptr_z   

      call this%material_coordinates%data(1)%Point_to_data (ptr_x)
      call this%material_coordinates%data(2)%Point_to_data (ptr_y)
      call this%material_coordinates%data(3)%Point_to_data (ptr_z)

   end subroutine Point_to_coordinates_3d


   subroutine Point_to_velocities(this, ptr_1, ptr_2)
      implicit none
      class (rezone_t)                         , intent(in out) :: this    
      real(8)       , dimension(:,:,:), pointer, intent(out)    :: ptr_1   
      real(8)       , dimension(:,:,:), pointer, intent(out)    :: ptr_2   

      call this%mesh_velocity%Point_to_data (ptr_1, ptr_2)

   end subroutine Point_to_velocities


   subroutine Point_to_volume(this, ptr)
      implicit none
      class (rezone_t)                         , intent(in out) :: this    
      real(8)       , dimension(:,:,:), pointer, intent(out)    :: ptr   

      call this%material_volume%Point_to_data (ptr)

   end subroutine Point_to_volume

   subroutine Clean_rezone (this)
      implicit none
      class (rezone_t), intent(in out) :: this 
      integer                        :: i

      deallocate (this%material_coordinates)

      call this%mesh_velocity%Clean_velocity()

   end subroutine Clean_rezone

   subroutine Set_communication(this, comm, comm_params_cell, comm_params_vertex)
      use communication_module, only : communication_t
      use communication_parameters_module, only : communication_parameters_t
      implicit none
      class (rezone_t)            :: this 
      type(communication_t), pointer            :: comm
      type(communication_parameters_t), pointer :: comm_params_cell, comm_params_vertex

      call this%material_coordinates%Set_communication(comm, comm_params_vertex)
      call this%mesh_velocity%Set_communication(comm, comm_params_vertex)
      call this%material_volume%Set_communication(comm, comm_params_cell)
   end subroutine Set_communication

end module rezone_module
