
module coordinates_module
   use vertex_quantity_module, only : vertex_quantity_t
   use vertex_boundary_condition_module, only : vertex_boundary_condition_t, vertex_bc_wrapper_t
   use general_utils_module, only : Get_axises_num
   use boundary_parameters_module      , only : boundary_parameters_t

   implicit none
   private
   public :: coordinates_t

   type, extends(vertex_quantity_t) :: coordinates_t
      private

      integer, private :: dimension
   contains


      procedure, public :: Calculate
      procedure, public :: Clean_coordinates
      procedure, public :: Apply_Boundary => Apply_coordinates_boundary

      procedure, public :: Write_coordinates
      procedure, public :: Write_quantity_abstract => Write_coordinates

      procedure, public :: Read_coordinates
      procedure, public :: Read_quantity_abstract => Read_coordinates

   end type


   interface coordinates_t

      module procedure Constructor
   end interface coordinates_t

contains

   type(coordinates_t) function Constructor(initial_data, d1, d2, d3, bc, bc_params)
      real(8)                                           , intent(in) :: initial_data 
      integer                                           , intent(in) :: d1           
      integer                                           , intent(in) :: d2           
      integer                                           , intent(in) :: d3           
      type (vertex_bc_wrapper_t) , dimension(:), pointer, intent(in out) :: bc           
      type(boundary_parameters_t), pointer, intent(in) :: bc_params

      integer                                                        :: axises_num

      axises_num = Get_axises_num (d1, d2, d3)
      Constructor%dimension = axises_num
      call Constructor%Init_vertex_quantity_init_val (initial_data, d1, d2, d3, axises_num, bc, bc_params)
   end function


   subroutine Clean_coordinates (this)
      class (coordinates_t), intent(in out) :: this 

      call this%Clean_vertex_quantity ()
   end subroutine Clean_coordinates

   subroutine Apply_coordinates_boundary(this, coordinates_values, is_blocking)
      use data_module, only : data_t
      implicit none
      class (coordinates_t)        , intent(in out) :: this
      type(data_t), dimension(:), pointer,intent(inout)        :: coordinates_values
      logical, optional :: is_blocking 
      integer, dimension(:)            , pointer        :: ed_num
      integer, dimension(:,:)          , pointer        :: b_map
      integer                                           :: i, b_size, edge
      integer                                           :: ind_i, ind_j, ind_k
      logical :: is_blocking_local

      if (.not. present(is_blocking)) then
         is_blocking_local = .True.  
      else
         is_blocking_local = is_blocking   
      end if

       call this%boundary_conditions(1)%bc% Calculate (this &
                                    ,   this%data, 0d0, 1)


        if (is_blocking_local) then
!            call this%Exchange_virtual_space_blocking()
        else
!            call this%Exchange_virtual_space_nonblocking()
        end if
   end subroutine Apply_coordinates_boundary

   subroutine Calculate(this, dt, velocity, coords)
      use vertex_quantity_module, only : vertex_quantity_t
      implicit none
      class (coordinates_t)                       , intent(inout)  :: this
      class (vertex_quantity_t)          , intent(inout)  :: velocity
      type (coordinates_t)   , pointer , optional, intent(inout)  :: coords
      real(8)                                     , intent(in)     :: dt

      integer :: i,j,k
      real(8), dimension(:, :, :), pointer :: velocity_x      
      real(8), dimension(:, :, :), pointer :: velocity_y      
      real(8), dimension(:, :, :), pointer :: velocity_z      
      real(8), dimension(:, :, :), pointer :: x      
      real(8), dimension(:, :, :), pointer :: y      
      real(8), dimension(:, :, :), pointer :: z      
      real(8), dimension(:, :, :), pointer :: x_tag      
      real(8), dimension(:, :, :), pointer :: y_tag      
      real(8), dimension(:, :, :), pointer :: z_tag      

      integer                                                        :: dimension


      if (this%dimension == 2) then
         call velocity%Point_to_data(velocity_x, velocity_y)
         call this    %Point_to_data(x         , y)

         if (.not. present(coords)) then
            call this  %Point_to_data(x_tag         , y_tag)
         else
            call coords%Point_to_data(x_tag         , y_tag)
         end if

         do j = 1, this%d2
            do i = 1, this%d1
               x(i, j, 1) = x_tag(i, j, 1) + dt * velocity_x(i, j, 1)
               y(i, j, 1) = y_tag(i, j, 1) + dt * velocity_y(i, j, 1)
            end do
         end do
      else
         call velocity%Point_to_data(velocity_x, velocity_y, velocity_z)
         call this    %Point_to_data(x         , y         , z)

         if (.not. present(coords)) then
            call this  %Point_to_data(x_tag, y_tag, z_tag)
         else
            call coords%Point_to_data(x_tag, y_tag, z_tag)
         end if
         do k = 1,this%d3
            do j = 1, this%d2
               do i = 1, this%d1
                  x(i, j, k) = x_tag(i, j, k) + dt * velocity_x(i, j, k)
                  y(i, j, k) = y_tag(i, j, k) + dt * velocity_y(i, j, k)
                  z(i, j, k) = z_tag(i, j, k) + dt * velocity_z(i, j, k)
               end do
            end do
         end do
      end if
   end subroutine Calculate



   subroutine Write_coordinates(this, unit, iostat, iomsg)
      class (coordinates_t), intent(in) :: this  
      integer,      intent(in)     :: unit
      integer,      intent(out)    :: iostat
      character(*), intent(in out)  :: iomsg

#ifdef DEBUG
      write(*,*) '@@@ in Write_coordinates @@@'
#endif

      call this%Write_vertex_quantity(unit, iostat=iostat, iomsg=iomsg)

#ifdef DEBUG
      write(*,*) '@@@ end Write_coordinates @@@'
#endif

   end subroutine Write_coordinates

   subroutine Read_coordinates(this, unit, iostat, iomsg)
      class (coordinates_t), intent(in out) :: this  
      integer,      intent(in)     :: unit
      integer,      intent(out)    :: iostat
      character(*), intent(in out)  :: iomsg

#ifdef DEBUG
      write(*,*) '@@@ in Read_coordinates @@@'
#endif

      call this%Read_vertex_quantity(unit, iostat=iostat, iomsg=iomsg)

#ifdef DEBUG
      write(*,*) '@@@ end Read_coordinates @@@'
#endif

   end subroutine Read_coordinates

end module coordinates_module
