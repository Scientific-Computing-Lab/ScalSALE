
module time_module
   use mesh_base_module          , only : mesh_base_t
   use velocity_module           , only : velocity_t
   use vof_module                , only : vof_t
   use vertex_mass_module        , only : vertex_mass_t
   use datafile_module           , only : datafile_t
   use geometry_module           , only : Distance_2
   use communication_module      , only : communication_t
   use parallel_parameters_module, only: parallel_parameters_t
   implicit none
   private

   type, public         :: time_t
      private
      class(mesh_base_t), pointer :: mesh     

      real(8), public        :: dt            
      real(8), public        :: dt_old        
      integer, public        :: current_cycle 
      real(8), public        :: time_passed   
      real(8), public        :: finish_time   
      real(8)                :: t0
      real(8), public        :: dt_rez        
      real(8)                :: dt_grad       
      real(8)                :: dt_grow       
      real(8)                :: dt_min        
      real(8)                :: dt_max        
      real(8), public        :: dt_cour       
      real(8), public        :: dt_mid        
      integer, public        :: i_cour        
      integer, public        :: j_cour        
      integer, public        :: k_cour        
      integer, public        :: i_rez         
      integer, public        :: j_rez         
      integer, public        :: k_rez         
      integer, public        :: i_grad        
      integer, public        :: j_grad        
      integer, public        :: k_grad        
      real(8), public        :: dt_cour_fac   
      real(8), public        :: emf           
      real(8), public        :: dt_factor     

      type(communication_t), pointer, public :: mpi_comm
      type(parallel_parameters_t), pointer, public :: parallel_params
      character(len=:), allocatable, public :: dt_reason     


   contains

      procedure, public  :: Calculate_dt

      procedure, public  :: Calculate_vel_grad_dt_2d

      procedure, public  :: Calculate_vel_grad_dt_3d

      procedure, public  :: Calculate_rezone_constraint_3d

      procedure, public  :: Calculate_rezone_constraint_2d

      procedure, public  :: Update_dt

      procedure, public  :: Update_time

      procedure, public  :: Should_continue

      procedure, public :: Write_time
      generic :: write(unformatted) => Write_time

      procedure, public :: Read_time
      generic :: read(unformatted) => Read_time


   end type time_t

   interface time_t
      module procedure Constructor
   end interface time_t

contains

   type(time_t) function Constructor(df, mpi_comm, parallel_params)
      type(datafile_t), intent (in) :: df
      type(communication_t),pointer, optional,intent (inout) :: mpi_comm
      type(parallel_parameters_t),pointer,optional, intent (inout) :: parallel_params

      Constructor%dt_max = df%dt_max
      Constructor%finish_time = df%finish_time
      Constructor%time_passed = df%t0
      Constructor%current_cycle = 0
      Constructor%dt = df%dt0
      Constructor%emf = df%emf
      Constructor%dt_grad = 1d20
      Constructor%dt_rez  = 1d20
      Constructor%dt_cour = 1d20
      Constructor%dt_cour_fac = df%dt_cour_fac
      Constructor%dt_factor = df%dt_factor
      if (present(mpi_comm)) then
         Constructor%mpi_comm => mpi_comm
      end if
      if (present(parallel_params)) then
         Constructor%parallel_params => parallel_params
      end if

   end function


   subroutine Update_dt(this, new_dt)
      class (time_t) :: this     
      real(8)        :: new_dt   

      this%dt = min(this%dt, new_dt)
   end subroutine Update_dt

   subroutine Update_time(this)
      class (time_t), intent(in out)         :: this 

      this%time_passed = this%time_passed + this%dt
      this%current_cycle = this%current_cycle + 1

   end subroutine Update_time

   logical function Should_continue(this)
      class (time_t), intent(in out)         :: this 

      Should_continue = .true.
      if (this%time_passed > this%finish_time) Should_continue = .false.
   end function Should_continue


   subroutine Calculate_dt(this,mesh, velocity, mesh_velocity,vertex_mass, vof, emfm)
      implicit none
      class (time_t)                      , intent(in out):: this     
      class(mesh_base_t)  , pointer       , intent(in)    :: mesh     
      type (velocity_t)   , pointer       , intent(in)    :: velocity 
      type (vof_t)        , pointer       , intent(in)    :: vof 
      type (velocity_t)   , pointer       , intent(in)    :: mesh_velocity     
      type (vertex_mass_t), pointer       , intent(inout) :: vertex_mass           
      real(8)                             , intent(in)    :: emfm

      real(8), dimension(:),allocatable :: coords_reason
      integer :: ierr


      integer :: i_dt_reason   
      integer :: j_dt_reason   
      integer :: k_dt_reason   







      if (mesh%dimension == 2) then
         call this%Calculate_vel_grad_dt_2d      (mesh, velocity, vof)
         call this%Calculate_rezone_constraint_2d(mesh, vertex_mass, velocity, mesh_velocity, emfm)
      else if(mesh%dimension == 3) then
         call this%Calculate_vel_grad_dt_3d      (mesh, velocity)
         call this%Calculate_rezone_constraint_3d(mesh, vertex_mass, velocity, mesh_velocity, emfm)
      end if

      if (this%current_cycle == 0) then
         this%dt_grow = this%dt 
      else
         this%dt_grow = 1.1d0 * this%dt
      end if

      this%dt_old = this%dt


      this%dt = min(this%dt_grow, this%dt_max, this%dt_grad, this%dt_cour, this%dt_rez) 



      i_dt_reason = 0
      j_dt_reason = 0
      k_dt_reason = 0

      if (this%dt == this%dt_max)  this%dt_reason = "Max"
      if (this%dt == this%dt_grow) this%dt_reason = "Grow"
      if (this%dt == this%dt_grad) then
         this%dt_reason = "Vel Grad"
         i_dt_reason    = this%i_grad
         j_dt_reason    = this%j_grad
         k_dt_reason    = this%k_grad
      end if
      if (this%dt == this%dt_rez ) then
         this%dt_reason = "Rezone"
         i_dt_reason    = this%i_rez
         j_dt_reason    = this%j_rez
         k_dt_reason    = this%k_rez
      end if
      if (this%dt == this%dt_cour) then
         this%dt_reason   = "Courant"
         i_dt_reason = this%i_cour
         j_dt_reason = this%j_cour
         k_dt_reason = this%k_cour
      end if





      allocate(coords_reason(3))
      coords_reason(1) = i_dt_reason
      coords_reason(2) = j_dt_reason
      coords_reason(3) = k_dt_reason
      call this%mpi_comm%Send_get_array_by_min_val(this%dt, coords_reason)
      deallocate(coords_reason)

      this%dt_mid = 0.5d0 * (this%dt + this%dt_old)

      if (this%current_cycle == 1) then      
         this%dt_min = this%dt * 1d-4
      endif

      if (this%dt < this%dt_min) then 




      end if


   end subroutine Calculate_dt

   subroutine Calculate_vel_grad_dt_2d(this, mesh, velocity, vof)
      class (time_t)               , intent(in out):: this     
      class(mesh_base_t), pointer  , intent(in)    :: mesh     
      type (velocity_t), pointer   , intent(in)    :: velocity 
      type (vof_t)     , pointer   , intent(in)    :: vof      

      integer :: nx, ny  

      real(8), dimension(:,:,:), pointer :: x     
      real(8), dimension(:,:,:), pointer :: y     
      real(8), dimension(:,:,:), pointer :: velocity_x       
      real(8), dimension(:,:,:), pointer :: total_vof        
      real(8), dimension(:,:,:), pointer :: velocity_y       

      integer :: i, j          
      real(8) :: dist_sq_i     
      real(8) :: dist_sq_j     
      real(8) :: vel_diff_sq_i 
      real(8) :: vel_diff_sq_j 
      real(8) :: vel_grad      
      real(8) :: vel_grad_max  

      call mesh%Point_to_data(x, y)
      call velocity%Point_to_data(velocity_x, velocity_y)
      call vof%Point_to_data(total_vof)

      vel_grad_max = 0d0

      this%dt_grad = 1d20
      nx = velocity%d1 - 1
      ny = velocity%d2 - 1

      do j = 1, ny
         do i = 1, nx

            if (total_vof(i,j,1) < this%emf) cycle

            dist_sq_i = 1d0 / Distance_2(x(i, j, 1), y(i, j, 1), x(i+1, j, 1), y(i+1, j, 1))
            dist_sq_j = 1d0 / Distance_2(x(i, j, 1), y(i, j, 1), x(i, j+1, 1), y(i, j+1, 1))

            vel_diff_sq_i = (velocity_x(i+1, j, 1) - velocity_x(i, j, 1)) ** 2 + (velocity_y(i+1, j, 1) - velocity_y(i, j, 1)) ** 2  
            vel_diff_sq_j = (velocity_x(i, j+1, 1) - velocity_x(i, j, 1)) ** 2 + (velocity_y(i, j+1, 1) - velocity_y(i, j, 1)) ** 2  

            vel_grad = sqrt(max(vel_diff_sq_i * dist_sq_i, vel_diff_sq_j * dist_sq_j)) 

            if (vel_grad >= vel_grad_max) then 
               vel_grad_max  = vel_grad
               this%i_grad = i
               this%j_grad = j
               this%k_grad = 1
            end if

         end do
      end do

      if (vel_grad_max /= 0) this%dt_grad = this%dt_factor / vel_grad_max 

   end subroutine Calculate_vel_grad_dt_2d

   subroutine Calculate_vel_grad_dt_3d(this, mesh, velocity)
      class (time_t)               , intent(in out):: this     
      class(mesh_base_t), pointer  , intent(in)    :: mesh     
      type (velocity_t), pointer   , intent(in)    :: velocity 

      integer :: nx, ny, nz  

      real(8), dimension(:,:,:), pointer :: x     
      real(8), dimension(:,:,:), pointer :: y     
      real(8), dimension(:,:,:), pointer :: z     
      real(8), dimension(:,:,:), pointer :: velocity_x       
      real(8), dimension(:,:,:), pointer :: velocity_y       
      real(8), dimension(:,:,:), pointer :: velocity_z       


      integer :: i, j, k       
      real(8) :: vel_grad      
      real(8) :: vel_grad_max  
      real(8) :: vel_diff_coor_diff_i  
      real(8) :: vel_diff_coor_diff_j  
      real(8) :: vel_diff_coor_diff_k  

      call mesh%Point_to_data(x, y, z)
      call velocity%Point_to_data(velocity_x, velocity_y, velocity_z)

      vel_grad_max = 0d0

      this%dt_grad = 1d20

      nx = velocity%d1 - 1
      ny = velocity%d2 - 1
      nz = velocity%d3 - 1

      do k = 1, nz
         do j = 1, ny
            do i = 1, nx

               vel_diff_coor_diff_i = ((velocity_x(i+1, j, k) - velocity_x(i, j, k)) * (x(i+1, j, k) - x(i, j, k)) &
                                      +(velocity_y(i+1, j, k) - velocity_y(i, j, k)) * (y(i+1, j, k) - y(i, j, k)) &
                                      +(velocity_z(i+1, j, k) - velocity_z(i, j, k)) * (z(i+1, j, k) - z(i, j, k))) / &
                                      ((x(i+1, j, k) - x(i, j, k)) ** 2 &
                                      +(y(i+1, j, k) - y(i, j, k)) ** 2 &
                                      +(z(i+1, j, k) - z(i, j, k)) ** 2)

               vel_diff_coor_diff_j = ((velocity_x(i, j+1, k) - velocity_x(i, j, k)) * (x(i, j+1, k) - x(i, j, k)) &
                                      +(velocity_y(i, j+1, k) - velocity_y(i, j, k)) * (y(i, j+1, k) - y(i, j, k)) &
                                      +(velocity_z(i, j+1, k) - velocity_z(i, j, k)) * (z(i, j+1, k) - z(i, j, k))) / &
                                      ((x(i, j+1, k) - x(i, j, k)) ** 2 &
                                      +(y(i, j+1, k) - y(i, j, k)) ** 2 &
                                      +(z(i, j+1, k) - z(i, j, k)) ** 2)

               vel_diff_coor_diff_k = ((velocity_x(i, j, k+1) - velocity_x(i, j, k)) * (x(i, j, k+1) - x(i, j, k)) &
                                      +(velocity_y(i, j, k+1) - velocity_y(i, j, k)) * (y(i, j, k+1) - y(i, j, k)) &
                                      +(velocity_z(i, j, k+1) - velocity_z(i, j, k)) * (z(i, j, k+1) - z(i, j, k))) / &
                                      ((x(i, j, k+1) - x(i, j, k)) ** 2 &
                                      +(y(i, j, k+1) - y(i, j, k)) ** 2 &
                                      +(z(i, j, k+1) - z(i, j, k)) ** 2)

               vel_grad = max(abs(vel_diff_coor_diff_i), abs(vel_diff_coor_diff_j), abs(vel_diff_coor_diff_k))

               if (vel_grad >= vel_grad_max) then 

                  vel_grad_max  = vel_grad
                  this%i_grad = i
                  this%j_grad = j
                  this%k_grad = k
               end if

            end do
         end do
      end do

      if (vel_grad_max /= 0) this%dt_grad = this%dt_factor / vel_grad_max 
   end subroutine Calculate_vel_grad_dt_3d

   subroutine Calculate_rezone_constraint_2d(this, mesh, vert_mass, velocity, mesh_velocity, emfm)
      implicit none
      class(time_t)                    , intent(inout)    :: this     
      class(mesh_base_t)      , pointer, intent(in)       :: mesh     
      type (velocity_t)       , pointer, intent(in)       :: velocity     
      type (velocity_t)       , pointer, intent(in)       :: mesh_velocity     
      type (vertex_mass_t)    , pointer, intent(inout)    :: vert_mass           
      real(8)                          , intent(in)       :: emfm

      integer :: i, j, ii 
      real(8) :: vel_diff_x   
      real(8) :: vel_diff_y   
      real(8) :: vel_diff     
      real(8) :: vel_grad     
      real(8) :: vel_rez_max  
      real(8), dimension(4) :: dx  
      real(8), dimension(4) :: dy  
      real(8), dimension(4) :: inv_dist_sq  
      real(8), dimension(:, :, :), pointer :: x               
      real(8), dimension(:, :, :), pointer :: y               
      real(8), dimension(:, :, :), pointer :: velocity_x      
      real(8), dimension(:, :, :), pointer :: velocity_y      
      real(8), dimension(:, :, :), pointer :: mesh_velocity_x 
      real(8), dimension(:, :, :), pointer :: mesh_velocity_y 
      real(8), dimension(:, :, :), pointer :: vertex_mass           
      integer :: nxp, nyp

      call mesh         %Point_to_data(x, y)
      call velocity     %Point_to_data(velocity_x, velocity_y)
      call mesh_velocity%Point_to_data(mesh_velocity_x, mesh_velocity_y)
      call vert_mass    %Point_to_data(vertex_mass)

      nxp = vert_mass%d1
      nyp = vert_mass%d2

      vel_rez_max = -1d20

      do j = 1, nyp
         do i = 1, nxp
            if(j /= 1 .and. vertex_mass(i, j - 1, 1) >= emfm) then

               if (i == nxp) then
                  dx(1) = 0d0
                  dy(1) = 0d0
                  inv_dist_sq(1) = 0d0
               else
                  dx(1) = x(i + 1, j - 1, 1) - x(i, j - 1, 1)
                  dy(1) = y(i + 1, j - 1, 1) - y(i, j - 1, 1)
                  inv_dist_sq(1) = 1d0 / (dx(1) * dx(1) + dy(1) * dy(1))
               end if

               if (i == 1) then
                  dx(2) = 0d0
                  dy(2) = 0d0
                  inv_dist_sq(2) = 0d0
               else
                  dx(2) = x(i - 1, j - 1, 1) - x(i, j - 1, 1)
                  dy(2) = y(i - 1, j - 1, 1) - y(i, j - 1, 1)
                  inv_dist_sq(2) = 1d0 / (dx(2) * dx(2) + dy(2) * dy(2))
               end if

               dx(3) = x(i, j, 1) - x(i, j - 1, 1)
               dy(3) = y(i, j, 1) - y(i, j - 1, 1)
               inv_dist_sq(3) = 1d0 / (dx(3) * dx(3) + dy(3) * dy(3))

               dx(4) = 0d0
               dy(4) = 0d0
               inv_dist_sq(4) = 0d0


                     if (j /= 2) then

                        dx(4) = x(i, j - 2, 1) - x(i, j - 1, 1)
                        dy(4) = y(i, j - 2, 1) - y(i, j - 1, 1)
                        inv_dist_sq(4) = 1d0 / (dx(4) * dx(4) + dy(4) * dy(4))
                     end if


                 vel_diff_x = mesh_velocity_x(i, j - 1, 1) - velocity_x(i, j - 1, 1)
                 vel_diff_y = mesh_velocity_y(i, j - 1, 1) - velocity_y(i, j - 1, 1)


               do ii = 1, 4
                  vel_diff = (vel_diff_x * dx(ii) + vel_diff_y * dy(ii)) * sqrt(inv_dist_sq(ii))
                  vel_grad = sqrt(vel_diff * vel_diff * inv_dist_sq(ii))
                  if (vel_rez_max <= vel_grad) then
                     vel_rez_max = vel_grad
                     this%i_rez  = i
                     this%j_rez  = j - 1
                  end if
               end do
            end if
         end do
      end do

      this%dt_rez = this%dt_factor / (vel_rez_max + 1d-30)
      this%dt_rez = max(this%dt_rez, 1d-14)

   end subroutine Calculate_rezone_constraint_2d

   subroutine Calculate_rezone_constraint_3d(this, mesh, vert_mass, velocity, mesh_velocity, emfm)
      implicit none
      class(time_t)                    , intent(inout)    :: this     
      class(mesh_base_t)      , pointer, intent(in)       :: mesh     
      type (velocity_t)       , pointer, intent(in)       :: velocity     
      type (velocity_t)       , pointer, intent(in)       :: mesh_velocity     
      type (vertex_mass_t)    , pointer, intent(inout)    :: vert_mass           
      real(8)                          , intent(in)       :: emfm

      real(8), dimension(:, :, :), pointer :: x               
      real(8), dimension(:, :, :), pointer :: y               
      real(8), dimension(:, :, :), pointer :: z               
      real(8), dimension(:, :, :), pointer :: velocity_x      
      real(8), dimension(:, :, :), pointer :: velocity_y      
      real(8), dimension(:, :, :), pointer :: velocity_z      
      real(8), dimension(:, :, :), pointer :: mesh_velocity_x 
      real(8), dimension(:, :, :), pointer :: mesh_velocity_y 
      real(8), dimension(:, :, :), pointer :: mesh_velocity_z 
      real(8), dimension(:, :, :), pointer :: vertex_mass     

      real(8) :: vel_diff_x   
      real(8) :: vel_diff_y   
      real(8) :: vel_diff_z   
      real(8) :: vel_grad     
      real(8) :: vel_rez_max  

      real(8), parameter :: eps = 1d-20  

      real(8) :: vel_grad_i_p  
      real(8) :: vel_grad_i_m  
      real(8) :: vel_grad_j_p  
      real(8) :: vel_grad_j_m  
      real(8) :: vel_grad_k_p  
      real(8) :: vel_grad_k_m  

      real(8) :: len_sq        
      integer :: nxp, nyp, nzp

      integer :: i, j, k 
      logical :: is_wall_x_top, is_wall_x_bot,is_wall_y_top, is_wall_y_bot,is_wall_z_top,is_wall_z_bot

      is_wall_x_top = this%parallel_params%is_wall_x_top
      is_wall_x_bot = this%parallel_params%is_wall_x_bot

      is_wall_y_top = this%parallel_params%is_wall_y_top
      is_wall_y_bot = this%parallel_params%is_wall_y_bot

      is_wall_z_top = this%parallel_params%is_wall_z_top
      is_wall_z_bot = this%parallel_params%is_wall_z_bot


      call velocity         %Point_to_data(velocity_x, velocity_y, velocity_z)
      call mesh             %Point_to_data(x, y, z)
      call vert_mass%Point_to_data(vertex_mass)
      call mesh_velocity  %Point_to_data(mesh_velocity_x, mesh_velocity_y, mesh_velocity_z)

      nxp = vert_mass%d1
      nyp = vert_mass%d2
      nzp = vert_mass%d3





      vel_rez_max = -1d20
      do k = 1, nzp
         do j = 1, nyp
            do i = 1, nxp
               if(vertex_mass(i, j, k) < emfm) cycle
               vel_diff_x = velocity_x(i, j, k) - mesh_velocity_x(i, j, k)
               vel_diff_y = velocity_y(i, j, k) - mesh_velocity_y(i, j, k)
               vel_diff_z = velocity_z(i, j, k) - mesh_velocity_z(i, j, k)

               if (i < nxp .or. is_wall_x_top .eqv. .false.) then
                  len_sq = (x(i+1, j, k) - x(i, j, k)) ** 2 &
                         + (y(i+1, j, k) - y(i, j, k)) ** 2 &
                         + (z(i+1, j, k) - z(i, j, k)) ** 2 + 1d-8

                  vel_grad_i_p = abs(vel_diff_x * (x(i+1, j, k) - x(i, j, k)) &
                                    +vel_diff_y * (y(i+1, j, k) - y(i, j, k)) &
                                    +vel_diff_z * (z(i+1, j, k) - z(i, j, k)) ) / len_sq
               else
                  vel_grad_i_p = eps
               end if

               if (j < nyp .or. is_wall_y_top .eqv. .false.) then
                  len_sq = (x(i, j+1, k) - x(i, j, k)) ** 2 &
                         + (y(i, j+1, k) - y(i, j, k)) ** 2 &
                         + (z(i, j+1, k) - z(i, j, k)) ** 2 + 1d-8

                  vel_grad_j_p = abs(vel_diff_x * (x(i, j+1, k) - x(i, j, k)) &
                                    +vel_diff_y * (y(i, j+1, k) - y(i, j, k)) &
                                    +vel_diff_z * (z(i, j+1, k) - z(i, j, k)) ) / len_sq
               else
                  vel_grad_j_p = eps
               end if

               if (k < nzp .or. is_wall_z_top .eqv. .false.) then
                  len_sq = (x(i, j, k+1) - x(i, j, k)) ** 2 &
                         + (y(i, j, k+1) - y(i, j, k)) ** 2 &
                         + (z(i, j, k+1) - z(i, j, k)) ** 2 + 1d-8

                  vel_grad_k_p = abs(vel_diff_x * (x(i, j, k+1) - x(i, j, k)) &
                                    +vel_diff_y * (y(i, j, k+1) - y(i, j, k)) &
                                    +vel_diff_z * (z(i, j, k+1) - z(i, j, k)) ) / len_sq
               else
                  vel_grad_k_p = eps
               end if

               if (i > 1 .or. is_wall_x_bot .eqv. .false.) then
                  len_sq = (x(i-1, j, k) - x(i, j, k)) ** 2 &
                         + (y(i-1, j, k) - y(i, j, k)) ** 2 &
                         + (z(i-1, j, k) - z(i, j, k)) ** 2 + 1d-8

                  vel_grad_i_m = abs(vel_diff_x * (x(i-1, j, k) - x(i, j, k)) &
                                    +vel_diff_y * (y(i-1, j, k) - y(i, j, k)) &
                                    +vel_diff_z * (z(i-1, j, k) - z(i, j, k)) ) / len_sq
               else
                  vel_grad_i_m = eps
               end if

               if (j > 1 .or. is_wall_y_bot .eqv. .false.) then
                  len_sq = (x(i, j-1, k) - x(i, j, k)) ** 2 &
                         + (y(i, j-1, k) - y(i, j, k)) ** 2 &
                         + (z(i, j-1, k) - z(i, j, k)) ** 2 + 1d-8

                  vel_grad_j_m = abs(vel_diff_x * (x(i, j-1, k) - x(i, j, k)) &
                                    +vel_diff_y * (y(i, j-1, k) - y(i, j, k)) &
                                    +vel_diff_z * (z(i, j-1, k) - z(i, j, k)) ) / len_sq
               else
                  vel_grad_j_m = eps
               end if

               if (k > 1 .or. is_wall_z_bot .eqv. .false.) then
                  len_sq = (x(i, j, k-1) - x(i, j, k)) ** 2 &
                         + (y(i, j, k-1) - y(i, j, k)) ** 2 &
                         + (z(i, j, k-1) - z(i, j, k)) ** 2 + 1d-8

                  vel_grad_k_m = abs(vel_diff_x * (x(i, j, k-1) - x(i, j, k)) &
                                    +vel_diff_y * (y(i, j, k-1) - y(i, j, k)) &
                                    +vel_diff_z * (z(i, j, k-1) - z(i, j, k)) ) / len_sq
               else
                  vel_grad_k_m = eps
               end if

               vel_grad = max(vel_grad_i_m, vel_grad_j_m, vel_grad_k_m, vel_grad_i_p, vel_grad_j_p, vel_grad_k_p)
               if (vel_rez_max < vel_grad) then
                  vel_rez_max = vel_grad
                  this%i_rez  = i
                  this%j_rez  = j
                  this%k_rez  = k
               end if
            end do
         end do
      end do




         this%dt_rez = this%dt_factor / (vel_rez_max + 1d-30)
         this%dt_rez = max(this%dt_rez, 1d-14)



   end subroutine Calculate_rezone_constraint_3d


   subroutine Write_time(this, unit, iostat, iomsg)
      class (time_t), intent(in) :: this  
      integer,      intent(in)     :: unit
      integer,      intent(out)    :: iostat
      character(*), intent(in out)  :: iomsg

#ifdef DEBUG
      write(*,*) '@@@ in Write_time @@@'
#endif

      write(unit, iostat=iostat, iomsg=iomsg) &
         len(this%dt_reason)

      write(unit, iostat=iostat, iomsg=iomsg) &
         this%dt, &
         this%dt_old, &
         this%current_cycle, &
         this%time_passed, &
         this%dt_rez, &
         this%dt_grad, &
         this%dt_grow, &
         this%dt_min, &
         this%dt_max, &
         this%dt_cour, &
         this%dt_mid, &
         this%i_cour, &
         this%j_cour, &
         this%k_cour, &
         this%i_rez, &
         this%j_rez, &
         this%k_rez, &
         this%i_grad, &
         this%j_grad, &
         this%k_grad, &
         this%dt_reason

#ifdef DEBUG
      write(*,*) '@@@ end Write_time @@@'

     write(*,*) &
         'len_dt_reason', &
         len(this%dt_reason)

      write(*,*) &
         'dt', &
         this%dt, &
         'dt_old', &
         this%dt_old, &
         'current_cycle', &
         this%current_cycle, &
         'time_passed', &
         this%time_passed, &
         'dt_rez', &
         this%dt_rez, &
         'dt_grad', &
         this%dt_grad, &
         'dt_grow', &
         this%dt_grow, &
         'dt_min', &
         this%dt_min, &
         'dt_max', &
         this%dt_max, &
         'dt_cour', &
         this%dt_cour, &
         'dt_mid', &
         this%dt_mid, &
         'i_cour', &
         this%i_cour, &
         'j_cour', &
         this%j_cour, &
         'k_cour', &
         this%k_cour, &
         'i_rez', &
         this%i_rez, &
         'j_rez', &
         this%j_rez, &
         'k_rez', &
         this%k_rez, &
         'i_grad', &
         this%i_grad, &
         'j_grad', &
         this%j_grad, &
         'k_grad', &
         this%k_grad, &
         'dt_reason', &
         this%dt_reason, &
         '###'
#endif

   end subroutine Write_time

   subroutine Read_time(this, unit, iostat, iomsg)
      class (time_t), intent(in out) :: this  
      integer,      intent(in)     :: unit
      integer,      intent(out)    :: iostat
      character(*), intent(in out)  :: iomsg
      integer :: len_dt_reason

#ifdef DEBUG
      write(*,*) '@@@ in Read_time @@@'
#endif

      read(unit, iostat=iostat, iomsg=iomsg) &
         len_dt_reason

      allocate(character(len=len_dt_reason) :: this%dt_reason)

      read(unit, iostat=iostat, iomsg=iomsg) &
         this%dt, &
         this%dt_old, &
         this%current_cycle, &
         this%time_passed, &
         this%dt_rez, &
         this%dt_grad, &
         this%dt_grow, &
         this%dt_min, &
         this%dt_max, &
         this%dt_cour, &
         this%dt_mid, &
         this%i_cour, &
         this%j_cour, &
         this%k_cour, &
         this%i_rez, &
         this%j_rez, &
         this%k_rez, &
         this%i_grad, &
         this%j_grad, &
         this%k_grad, &
         this%dt_reason

#ifdef DEBUG
      write(*,*) &
         'len_dt_reason', &
         len_dt_reason

      write(*,*) &
         'dt', &
         this%dt, &
         'dt_old', &
         this%dt_old, &
         'current_cycle', &
         this%current_cycle, &
         'time_passed', &
         this%time_passed, &
         'dt_rez', &
         this%dt_rez, &
         'dt_grad', &
         this%dt_grad, &
         'dt_grow', &
         this%dt_grow, &
         'dt_min', &
         this%dt_min, &
         'dt_max', &
         this%dt_max, &
         'dt_cour', &
         this%dt_cour, &
         'dt_mid', &
         this%dt_mid, &
         'i_cour', &
         this%i_cour, &
         'j_cour', &
         this%j_cour, &
         'k_cour', &
         this%k_cour, &
         'i_rez', &
         this%i_rez, &
         'j_rez', &
         this%j_rez, &
         'k_rez', &
         this%k_rez, &
         'i_grad', &
         this%i_grad, &
         'j_grad', &
         this%j_grad, &
         'k_grad', &
         this%k_grad, &
         'dt_reason', &
         this%dt_reason, &
         '###'

      write(*,*) '@@@ end Read_time @@@'
#endif

   end subroutine Read_time

end module time_module
