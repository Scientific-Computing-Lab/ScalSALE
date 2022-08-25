
module data_module
   use communication_module, only : communication_t
   use communication_parameters_module, only : communication_parameters_t
   use parallel_parameters_module, only: parallel_parameters_t
   implicit none
   private
   public :: data_t


   type :: data_t
      private


      real(8), dimension(:,:,:), pointer, public            :: values
      integer, public                                       :: nx   
      integer, public                                       :: ny   
      integer, public                                       :: nz   
      type(communication_t), pointer                             :: communication
      type(communication_parameters_t), pointer :: communication_parameters
      type (parallel_parameters_t)   , public, pointer :: parallel_params      


      real(8), dimension(:), allocatable     :: send_buf
      real(8), dimension(:), allocatable     :: recv_buf
      integer :: request
   contains


      procedure, public :: Point_to_data => Ptr_data

      procedure, public :: Clean_data

      procedure, public :: Set_communication

      procedure, public :: Debug_check_nan

      procedure, public :: Exchange_virtual_space_blocking


      procedure, private ::Set_send_buf

      procedure, private ::Get_recv_buf


      procedure, public :: Exchange_virtual_space_nonblocking

      procedure, public :: Exchange_end

      procedure, public :: debug_print

      procedure, public :: Write_data
      generic :: write(unformatted) => Write_data


      procedure, public :: Read_data
      generic :: read(unformatted) => Read_data

   end type data_t

   public :: Get_copy

   interface data_t
      module procedure Constructor_init_arr
      module procedure Constructor_no_init
      module procedure Constructor_init_val

      module procedure Constructor_is_parallel
      module procedure Constructor_init_arr_is_parallel
      module procedure Constructor_no_init_is_parallel
      module procedure Constructor_init_val_is_parallel

   end interface data_t

contains

   type(data_t) function Constructor_init_arr(initial_data, d1, d2, d3)
      implicit none
      real(8), dimension(:,:,:), intent(in) :: initial_data 
      integer                  , intent(in) :: d1           
      integer                  , intent(in) :: d2           
      integer                  , intent(in) :: d3           
      integer, dimension(3) :: vals_shape

         allocate (Constructor_init_arr%values (0:d1, 0:d2, 0:d3))
      Constructor_init_arr%values(1:, 1:, 1:) =  initial_data
      Constructor_init_arr%nx = d1
      Constructor_init_arr%ny = d2
      Constructor_init_arr%nz = d3
      if (d1 == -1 .or. d2 == -1 .or. d3 == -1) write(*,*)"data is -1"
   end function

   type(data_t) function Constructor_init_val(initial_val, d1, d2, d3)
      implicit none
      real(8)           , intent(in) :: initial_val  
      integer           , intent(in) :: d1           
      integer           , intent(in) :: d2           
      integer           , intent(in) :: d3           

      allocate (Constructor_init_val%values (0:d1, 0:d2, 0:d3))
      Constructor_init_val%values = initial_val
      Constructor_init_val%nx = d1
      Constructor_init_val%ny = d2
      Constructor_init_val%nz = d3
      if (d1 == -1 .or. d2 == -1 .or. d3 == -1) write(*,*)"data is -1"


   end function

   type(data_t) function Constructor_no_init(d1, d2, d3)
      implicit none
      integer           , intent(in) :: d1           
      integer           , intent(in) :: d2           
      integer           , intent(in) :: d3           

      allocate (Constructor_no_init%values (0:d1, 0:d2, 0:d3))
      Constructor_no_init%nx = d1
      Constructor_no_init%ny = d2
      Constructor_no_init%nz = d3
            if (d1 == -1 .or. d2 == -1 .or. d3 == -1) write(*,*)"data is -1"

   end function

   type(data_t) function Constructor_init_arr_is_parallel(initial_data, d1, d2, d3, is_parallel)
      implicit none
      real(8), dimension(:,:,:), intent(in) :: initial_data 
      integer                  , intent(in) :: d1           
      integer                  , intent(in) :: d2           
      integer                  , intent(in) :: d3           
      logical                               :: is_parallel

   end function

   type(data_t) function Constructor_no_init_is_parallel(d1, d2, d3, is_parallel)
      implicit none
      integer           , intent(in) :: d1           
      integer           , intent(in) :: d2           
      integer           , intent(in) :: d3           
      logical                        :: is_parallel

   end function

   type(data_t) function Constructor_init_val_is_parallel(initial_val, d1, d2, d3, is_parallel)
      implicit none
      real(8)           , intent(in) :: initial_val  
      integer           , intent(in) :: d1           
      integer           , intent(in) :: d2           
      integer           , intent(in) :: d3           
      logical                        :: is_parallel

   end function

   type(data_t) function Constructor_is_parallel(is_parallel)
      implicit none
      logical, optional :: is_parallel
      if (.not. present(is_parallel) .or. is_parallel .eqv. .true.) then
      else
      end if

   end function



   subroutine Ptr_data (this, ptr)
      class (data_t)                    , intent(in)  :: this 
      real(8), dimension(:,:,:), pointer, intent(out) :: ptr  

      ptr => this%values
   end subroutine Ptr_data

   function Get_copy (this)
      class (data_t)       , intent(in)  :: this 
      real(8), dimension(:,:,:), pointer :: Get_copy  

      Get_copy = this%values
   end function Get_copy

   subroutine Clean_data (this)
      class (data_t), intent(in out) :: this  

      deallocate (this%values)
   end subroutine Clean_data


   subroutine Set_communication (this, comm, comm_params)
      class (data_t), intent(in out) :: this  
      type(communication_t), pointer            :: comm
      type(communication_parameters_t), pointer :: comm_params
      integer, dimension(3) :: vals_shape
      integer :: d1,d2,d3

      this%communication => comm
      this%communication_parameters => comm_params

      this%parallel_params => this%communication%parallel_params
      if (this%communication%is_parallel .eqv. .true.) then
         vals_shape = shape(this%values)
         d1 = vals_shape(1) - 2
         d2 = vals_shape(2) - 2
         d3 = vals_shape(3) - 2
         allocate(this%send_buf(0:2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+4*(d2+2)+4*(d1+2)+7))
         allocate(this%recv_buf(0:2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+4*(d2+2)+4*(d1+2)+7))
      end if
   end subroutine Set_communication


   subroutine Exchange_virtual_space_blocking (this, ghost_width)
      class (data_t), intent(in out) :: this  
      real(8), dimension(:,:,:), allocatable :: send_buf, recv_buf
      integer, dimension(3) :: vals_shape
      integer, optional :: ghost_width
      integer :: ghost_width_local
      integer :: ghost_width_local_x
      integer :: i,j,k
      if (.not. present(ghost_width)) then
         ghost_width_local = 1   
      else
         ghost_width_local = 1   
      end if
      if (this%communication%is_parallel .eqv. .true.) then
         call this%Set_send_buf()
         call this%communication%Send_recv_neighbors_diag (this%communication_parameters, this%send_buf, this%recv_buf)
         call this%Get_recv_buf()
      end if


   end subroutine Exchange_virtual_space_blocking




   subroutine Exchange_virtual_space_nonblocking (this, ghost_width)
      class (data_t), intent(in out) :: this  
      integer, optional :: ghost_width
      integer :: ghost_width_local

      if (.not. present(ghost_width)) then
         ghost_width_local = 1   
      else
         ghost_width_local = 1   
      end if
      if (this%communication%is_parallel .eqv. .true.) then
         call this%Set_send_buf()
         call this%communication%Send_neighbors_diag (this%communication_parameters,&
                                                      this%send_buf, this%recv_buf, this%request)

      end if

   end subroutine Exchange_virtual_space_nonblocking

   subroutine Exchange_end (this)
      class (data_t), intent(in out) :: this  

      if (this%communication%is_parallel .eqv. .true.) then
         call this%communication%Wait_recv_neighbors_diag (this%communication_parameters,&
                                                           this%send_buf, this%recv_buf, this%request)
         call this%Get_recv_buf()
      end if

   end subroutine Exchange_end


   subroutine Get_recv_buf(this)
      class (data_t), intent(in out) :: this  
      integer, dimension(3) :: vals_shape
      integer :: d1, d2, d3, x, y, z

      vals_shape = shape(this%values)
      d1 = vals_shape(1) - 2
      d2 = vals_shape(2) - 2
      d3 = vals_shape(3) - 2

      x = this%parallel_params%my_coords(1)
      y = this%parallel_params%my_coords(2)
      z = this%parallel_params%my_coords(3)

      if (x+1 /= this%parallel_params%npx+1) then
         this%values(d1+1, 0:d2+1, 0:d3+1) = reshape(this%recv_buf(0 : (d2+2)*(d3+2)-1), (/d2+2, d3+2/))
      end if

      if (x-1 /= 0) then
         this%values(0, 0:d2+1, 0:d3+1) = reshape(this%recv_buf((d2+2)*(d3+2) : 2*(d2+2)*(d3+2)-1), (/d2+2, d3+2/))
      end if

      if (y+1 /= this%parallel_params%npy+1) then
         this%values(0:d1+1, d2+1, 0:d3+1) = reshape(this%recv_buf(2*(d2+2)*(d3+2) : 2*(d2+2)*(d3+2)+(d1+2)*(d3+2)-1), (/d1+2, d3+2/))
      end if
      if (y-1 /= 0) then
         this%values(0:d1+1, 0, 0:d3+1) = reshape(this%recv_buf(2*(d2+2)*(d3+2)+(d1+2)*(d3+2) : 2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)-1), (/d1+2, d3+2/))
      end if
      if (z+1 /= this%parallel_params%npz+1) then
         this%values(0:d1+1, 0:d2+1, d3+1) = reshape(this%recv_buf(2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2) : 2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+(d1+2)*(d2+2)-1), (/d1+2, d2+2/))

      end if
      if (z-1 /= 0) then
         this%values(0:d1+1, 0:d2+1, 0) = reshape(this%recv_buf(2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+(d1+2)*(d2+2) : 2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)-1), (/d1+2, d2+2/))
      end if

      if (x+1 /= this%parallel_params%npx+1 .and. y+1 /= this%parallel_params%npy+1) then
         this%values(d1+1, d2+1, 0:d3+1) = reshape(this%recv_buf(2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2) : 2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+(d3+2)-1), (/d3+2/))
      end if
      if (x+1 /= this%parallel_params%npx+1 .and. y-1 /= 0) then
         this%values(d1+1, 0, 0:d3+1) = reshape(this%recv_buf(2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+(d3+2) : 2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+2*(d3+2)-1), (/d3+2/))
      end if
      if (x-1 /= 0 .and. y+1 /= this%parallel_params%npy+1) then
         this%values(0, d2+1, 0:d3+1) = reshape(this%recv_buf(2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+2*(d3+2) : 2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+3*(d3+2)-1), (/d3+2/))
      end if
      if (x-1 /= 0 .and. y-1 /= 0) then
         this%values(0, 0, 0:d3+1) = reshape(this%recv_buf(2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+3*(d3+2) : 2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)-1), (/d3+2/))
      end if
      if (x+1 /= this%parallel_params%npx+1 .and. z+1 /= this%parallel_params%npz+1) then
         this%values(d1+1, 0:d2+1, d3+1) = reshape(this%recv_buf(2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2) : 2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+(d2+2)-1), (/d2+2/))
      end if
      if (x+1 /= this%parallel_params%npx+1 .and. z-1 /= 0) then
         this%values(d1+1, 0:d2+1, 0) = reshape(this%recv_buf(2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+(d2+2) : 2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+2*(d2+2)-1), (/d2+2/))
      end if
      if (x-1 /= 0 .and. z+1 /= this%parallel_params%npz+1) then
         this%values(0, 0:d2+1, d3+1) = reshape(this%recv_buf(2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+2*(d2+2) : 2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+3*(d2+2)-1), (/d2+2/))
      end if
      if (x-1 /= 0 .and. z-1 /= 0) then
         this%values(0, 0:d2+1, 0) = reshape(this%recv_buf(2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+3*(d2+2) : 2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+4*(d2+2)-1), (/d2+2/))
      end if
      if (y+1 /= this%parallel_params%npy+1 .and. z+1 /= this%parallel_params%npz+1) then
         this%values(0:d1+1, d2+1, d3+1) = reshape(this%recv_buf(2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+4*(d2+2) : 2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+4*(d2+2)+(d1+2)-1), (/d1+2/))
      end if
      if (y+1 /= this%parallel_params%npy+1 .and. z-1 /= 0) then
         this%values(0:d1+1, d2+1, 0) = reshape(this%recv_buf(2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+4*(d2+2)+(d1+2) : 2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+4*(d2+2)+2*(d1+2)-1), (/d1+2/))
      end if
      if (y-1 /= 0 .and. z+1 /= this%parallel_params%npz+1) then
         this%values(0:d1+1, 0, d3+1) = reshape(this%recv_buf(2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+4*(d2+2)+2*(d1+2) : 2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+4*(d2+2)+3*(d1+2)-1), (/d1+2/))
      end if
      if (y-1 /= 0 .and. z-1 /= 0) then
         this%values(0:d1+1, 0, 0) = reshape(this%recv_buf(2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+4*(d2+2)+3*(d1+2) : 2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+4*(d2+2)+4*(d1+2)-1), (/d1+2/))
      end if



      if (x+1 /= this%parallel_params%npx+1 .and. y+1 /= this%parallel_params%npy+1 .and. z+1 /= this%parallel_params%npz+1) then
         this%values(d1+1, d2+1, d3+1) = this%recv_buf(2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+4*(d2+2)+4*(d1+2))
      end if
      if (x+1 /= this%parallel_params%npx+1 .and. y+1 /= this%parallel_params%npy+1 .and. z-1 /= 0) then
         this%values(d1+1, d2+1, 0) = this%recv_buf(2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+4*(d2+2)+4*(d1+2)+1)
      end if
      if (x+1 /= this%parallel_params%npx+1 .and. y-1 /= 0 .and. z+1 /= this%parallel_params%npz+1) then
         this%values(d1+1, 0, d3+1) = this%recv_buf(2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+4*(d2+2)+4*(d1+2)+2)
      end if
      if (x+1 /= this%parallel_params%npx+1 .and. y-1 /= 0 .and. z-1 /= 0) then
         this%values(d1+1, 0, 0) = this%recv_buf(2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+4*(d2+2)+4*(d1+2)+3)
      end if
      if (x-1 /= 0 .and. y+1 /= this%parallel_params%npy+1 .and. z+1 /= this%parallel_params%npz+1) then
         this%values(0, d2+1, d3+1) = this%recv_buf(2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+4*(d2+2)+4*(d1+2)+4)
      end if
      if (x-1 /= 0 .and. y+1 /= this%parallel_params%npy+1 .and. z-1 /= 0) then
         this%values(0, d2+1, 0) = this%recv_buf(2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+4*(d2+2)+4*(d1+2)+5)
      end if
      if (x-1 /= 0 .and. y-1 /= 0 .and. z+1 /= this%parallel_params%npz+1) then
         this%values(0, 0, d3+1) = this%recv_buf(2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+4*(d2+2)+4*(d1+2)+6)
      end if
      if (x-1 /= 0 .and. y-1 /= 0 .and. z-1 /= 0) then
         this%values(0, 0, 0) = this%recv_buf(2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+4*(d2+2)+4*(d1+2)+7)
      end if
!
!      if (this%parallel_params%my_rank == 0) then
!      write(*,*)this%parallel_params%my_rank,this%values
!end if
   end subroutine Get_recv_buf


   subroutine Set_send_buf(this)
      class (data_t), intent(in out) :: this  
      integer, dimension(3) :: vals_shape
      integer :: d1, d2, d3, x, y, z, offset

      vals_shape = shape(this%values)
      d1 = vals_shape(1) - 2
      d2 = vals_shape(2) - 2
      d3 = vals_shape(3) - 2

      x = this%communication%parallel_params%my_coords(1)
      y = this%communication%parallel_params%my_coords(2)
      z = this%communication%parallel_params%my_coords(3)
      offset = this%communication_parameters%dim_offset

      if (x+1 /= this%parallel_params%npx+1) then
         this%send_buf(0 : (d2+2)*(d3+2)-1) = reshape(this%values(d1 - offset, 0:d2+1, 0:d3+1), (/d2+2 * d3+2/))
      end if

      if (x-1 /= 0) then
         this%send_buf((d2+2)*(d3+2) : 2*(d2+2)*(d3+2)-1) = reshape(this%values(1 + offset, 0:d2+1, 0:d3+1), (/d2+2 * d3+2/))
      end if

      if (y+1 /= this%parallel_params%npy+1) then
         this%send_buf(2*(d2+2)*(d3+2) : 2*(d2+2)*(d3+2)+(d1+2)*(d3+2)-1) = reshape(this%values(0:d1+1, d2 - offset, 0:d3+1), (/d1+2 * d3+2/))
      end if
      if (y-1 /= 0) then
         this%send_buf(2*(d2+2)*(d3+2)+(d1+2)*(d3+2) : 2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)-1) = reshape(this%values(0:d1+1, 1 + offset, 0:d3+1), (/d1+2 * d3+2/))
      end if
      if (z+1 /= this%parallel_params%npz+1) then
         this%send_buf(2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2) : 2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+(d1+2)*(d2+2)-1) = reshape(this%values(0:d1+1, 0:d2+1, d3 - offset), (/d1+2 * d2+2/))
      end if
      if (z-1 /= 0) then
         this%send_buf(2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+(d1+2)*(d2+2) : 2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)-1) = reshape(this%values(0:d1+1, 0:d2+1, 1 + offset), (/d1+2 * d2+2/))

      end if

      if (x+1 /= this%parallel_params%npx+1 .and. y+1 /= this%parallel_params%npy+1) then
         this%send_buf(2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2) : 2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+(d3+2)-1) = this%values(d1 - offset, d2 - offset, 0:d3+1)
      end if
      if (x+1 /= this%parallel_params%npx+1 .and. y-1 /= 0) then
         this%send_buf(2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+(d3+2) : 2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+2*(d3+2)-1) = this%values(d1 - offset, 1 + offset, 0:d3+1)
      end if
      if (x-1 /= 0 .and. y+1 /= this%parallel_params%npy+1) then
         this%send_buf(2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+2*(d3+2) : 2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+3*(d3+2)-1) = this%values(1 + offset, d2 - offset, 0:d3+1)
      end if
      if (x-1 /= 0 .and. y-1 /= 0) then
         this%send_buf(2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+3*(d3+2) : 2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)-1) = this%values(1+offset, 1+offset, 0:d3+1)
      end if
      if (x+1 /= this%parallel_params%npx+1 .and. z+1 /= this%parallel_params%npz+1) then
         this%send_buf(2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2) : 2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+(d2+2)-1) = this%values(d1 - offset, 0:d2+1, d3 - offset)
      end if
      if (x+1 /= this%parallel_params%npx+1 .and. z-1 /= 0) then
         this%send_buf(2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+(d2+2) : 2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+2*(d2+2)-1) = this%values(d1 - offset, 0:d2+1, 1 + offset)
      end if
      if (x-1 /= 0 .and. z+1 /= this%parallel_params%npz+1) then
         this%send_buf(2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+2*(d2+2) : 2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+3*(d2+2)-1) = this%values(1 + offset, 0:d2+1, d3 - offset)
      end if
      if (x-1 /= 0 .and. z-1 /= 0) then
         this%send_buf(2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+3*(d2+2) : 2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+4*(d2+2)-1) = this%values(1 + offset, 0:d2+1, 1 + offset)
      end if
      if (y+1 /= this%parallel_params%npy+1 .and. z+1 /= this%parallel_params%npz+1) then
         this%send_buf(2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+4*(d2+2) : 2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+4*(d2+2)+(d1+2)-1) = this%values(0:d1+1, d2 - offset, d3 - offset)
      end if
      if (y+1 /= this%parallel_params%npy+1 .and. z-1 /= 0) then
         this%send_buf(2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+4*(d2+2)+(d1+2) : 2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+4*(d2+2)+2*(d1+2)-1) = this%values(0:d1+1, d2 - offset, 1 + offset)
      end if
      if (y-1 /= 0 .and. z+1 /= this%parallel_params%npz+1) then
         this%send_buf(2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+4*(d2+2)+2*(d1+2) : 2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+4*(d2+2)+3*(d1+2)-1) = this%values(0:d1+1, 1 + offset, d3 - offset)
      end if
      if (y-1 /= 0 .and. z-1 /= 0) then
         this%send_buf(2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+4*(d2+2)+3*(d1+2) : 2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+4*(d2+2)+4*(d1+2)-1) = this%values(0:d1+1, 1 + offset, 1+offset)
      end if



      if (x+1 /= this%parallel_params%npx+1 .and. y+1 /= this%parallel_params%npy+1 .and. z+1 /= this%parallel_params%npz+1) then
         this%send_buf(2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+4*(d2+2)+4*(d1+2)) = this%values(d1- offset, d2-offset, d3-offset)
      end if
      if (x+1 /= this%parallel_params%npx+1 .and. y+1 /= this%parallel_params%npy+1 .and. z-1 /= 0) then
         this%send_buf(2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+4*(d2+2)+4*(d1+2)+1) = this%values(d1-offset, d2-offset, 1 + offset)
      end if
      if (x+1 /= this%parallel_params%npx+1 .and. y-1 /= 0 .and. z+1 /= this%parallel_params%npz+1) then
         this%send_buf(2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+4*(d2+2)+4*(d1+2)+2) = this%values(d1-offset, 1+offset, d3-offset)
      end if
      if (x+1 /= this%parallel_params%npx+1 .and. y-1 /= 0 .and. z-1 /= 0) then
         this%send_buf(2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+4*(d2+2)+4*(d1+2)+3) = this%values(d1-offset, 1+offset, 1+offset)
      end if
      if (x-1 /= 0 .and. y+1 /= this%parallel_params%npy+1 .and. z+1 /= this%parallel_params%npz+1) then
         this%send_buf(2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+4*(d2+2)+4*(d1+2)+4) = this%values(1+offset, d2-offset, d3-offset)
      end if
      if (x-1 /= 0 .and. y+1 /= this%parallel_params%npy+1 .and. z-1 /= 0) then
         this%send_buf(2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+4*(d2+2)+4*(d1+2)+5) = this%values(1+offset, d2-offset, 1+offset)
      end if
      if (x-1 /= 0 .and. y-1 /= 0 .and. z+1 /= this%parallel_params%npz+1) then
         this%send_buf(2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+4*(d2+2)+4*(d1+2)+6) = this%values(1+offset, 1+offset, d3-offset)
      end if
      if (x-1 /= 0 .and. y-1 /= 0 .and. z-1 /= 0) then
         this%send_buf(2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+4*(d2+2)+4*(d1+2)+7) = this%values(1+offset, 1+offset, 1+offset)
      end if
   end subroutine Set_send_buf


   subroutine Debug_check_nan(this, caller)
      implicit none
      class (data_t), intent(in out) :: this  
      CHARACTER(*) caller

      integer :: i,j,k
      do i = 0, this%nx
         do j = 0, this%ny
            do k = 0, this%nz
               if (this%values(i,j,k) /= this%values(i,j,k)) then
                  write(*,*) "NaN in ", caller, ": ",i,j,k
                  return
               end if
            end do
         end do
      end do
   end subroutine Debug_check_nan

   subroutine debug_print(this, caller, flag)
      implicit none
      class (data_t), intent(in out) :: this  
      integer, optional :: flag
      CHARACTER(*) caller
      integer :: width
      integer :: i,j,k
      if (.not. present(flag)) then
         width = 0   
      else
         width = flag   
      end if

      write(69,*) "---- data_t:", caller, "----"
      do i = width, this%nx - width
         do j = width, this%ny - width
            do k = width, this%nz - width
                  write(69,*) i,j,k,this%values(i,j,k)
            end do
         end do
      end do
   end subroutine debug_print

   subroutine Write_data(this, unit, iostat, iomsg)
      class (data_t), intent(in) :: this  
      integer,      intent(in)     :: unit
      integer,      intent(out)    :: iostat
      character(*), intent(in out)  :: iomsg

#ifdef DEBUG
      write(*,*) '@@@ in Write_data @@@'
#endif

      write(unit, iostat=iostat, iomsg=iomsg) &
         shape(this%values)

      write(unit, iostat=iostat, iomsg=iomsg) &
         this%values, &
         this%nx, &
         this%ny, &
         this%nz

#ifdef DEBUG
      write(*,*) &
         'shape_values', &
         shape(this%values)
      write(*,*) &
         'values', &
         this%values, &
         'nx', &
         this%nx, &
         'ny', &
         this%ny, &
         'nz', &
         this%nz, &
         '###'

      write(*,*) '@@@ end Write_data @@@'
#endif

   end subroutine Write_data

   subroutine Read_data(this, unit, iostat, iomsg)
      class (data_t), intent(in out) :: this  
      integer,      intent(in)     :: unit
      integer,      intent(out)    :: iostat
      character(*), intent(in out)  :: iomsg
      integer, dimension(3) :: shape_values

#ifdef DEBUG
      write(*,*) "@@@ in Read_data @@@"
#endif

      read(unit, iostat=iostat, iomsg=iomsg) &
         shape_values

      deallocate(this%values)
      allocate(this%values(0:shape_values(1) - 1, 0:shape_values(2)-1, 0:shape_values(3) - 1))

      read(unit, iostat=iostat, iomsg=iomsg) &
         this%values, &
         this%nx, &
         this%ny, &
         this%nz

#ifdef DEBUG
      write(*,*) &
         'shape_values', &
         shape(this%values)
      write(*,*) &
         'values', &
         this%values, &
         'nx', &
         this%nx, &
         'ny', &
         this%ny, &
         'nz', &
         this%nz, &
         '###'

      write(*,*) "@@@ end Read_data @@@"
#endif

   end subroutine Read_data
end module data_module
