
module communication_module
   use communication_parameters_module, only : communication_parameters_t
   use mpi
   use communication_utils_module, only : MPI_User_fn_array_by_min_val, MPI_User_fn_array_by_max_val
   use parallel_parameters_module      , only : parallel_parameters_t

   implicit none
   private
   public :: communication_t


   type :: communication_t
      private
      logical, public :: is_parallel
      type (parallel_parameters_t)   , public,pointer :: parallel_params      

   contains


      procedure, public :: Send_blocking

      procedure, public :: Send_nonblocking

      procedure, public :: Send_recv_neighbors

      procedure, public :: Send_recv_neighbors_diag

      procedure, public :: Send_neighbors_diag

      procedure, public :: Wait_recv_neighbors_diag

      procedure, public :: Send_get_max

      procedure, public :: Send_get_min

      procedure, public :: Send_get_array_by_min_val

      procedure, public :: Send_get_array_by_max_val

      procedure, public :: Send_get_sum

      procedure, public :: Destroy

   end type

   interface communication_t
      module procedure Constructor
      module procedure Constructor_test

   end interface communication_t

   interface Get_rank
      module procedure rank
   end interface Get_rank

contains

   type(communication_t) function Constructor_test(use_mpi)
      implicit none
      logical :: use_mpi

      Constructor_test%is_parallel = use_mpi


   end function

   type(communication_t) function Constructor(parallel_p)
      implicit none
      type (parallel_parameters_t) , pointer, intent(in) :: parallel_p     

      Constructor%is_parallel = parallel_p%is_parallel
      Constructor%parallel_params => parallel_p


   end function

   subroutine Send_blocking (this, communication_parameters)
      class (communication_t)             , intent(in)      :: this 
      type(communication_parameters_t)                      :: communication_parameters

   end subroutine Send_blocking

   subroutine Send_nonblocking (this, communication_parameters)
      class (communication_t)             , intent(in)   :: this 
      type(communication_parameters_t)    , pointer      :: communication_parameters

      if (this%is_parallel .eqv. .true. .and. associated(communication_parameters)) then
      end if
   end subroutine Send_nonblocking

   subroutine Send_recv_neighbors (this, send_buf, recv_buf, comm_params)
      class (communication_t)             , intent(in)    :: this 
      real(8), dimension(:,:,:), allocatable, intent(in)     :: send_buf
      real(8), dimension(:,:,:), allocatable, intent(in out) :: recv_buf
      type(communication_parameters_t) , pointer, intent(in out) :: comm_params
      integer :: ierr

      if (this%is_parallel .eqv. .true.) then

         call MPI_neighbor_alltoallw(send_buf, &
                                     comm_params%sizes_send, &
                                     comm_params%start_send, &
                                     comm_params%types_send, &
                                     recv_buf, &
                                     comm_params%sizes_recv, &
                                     comm_params%start_recv, &
                                     comm_params%types_recv, &
                                     comm_params%comm, &
                                     ierr &
                                    )

      end if
   end subroutine Send_recv_neighbors


   subroutine Send_recv_neighbors_diag (this, comm_params, send_buf, recv_buf)
      class (communication_t)             , intent(in)       :: this 
      real(8), dimension(:), allocatable, intent(in out) :: send_buf
      real(8), dimension(:), allocatable, intent(in out) :: recv_buf
      integer, dimension(3) :: vals_shape
      integer  :: d1, d2, d3
      integer  :: x, y, z, i,j,k,c
      type(communication_parameters_t) , pointer, intent(in out) :: comm_params
      integer :: ierr
      integer :: my_rank
      integer :: result_len
      character (mpi_max_error_string) :: err_string



      if (this%is_parallel .eqv. .true.) then


         call MPI_neighbor_alltoallw(send_buf, &
                                     comm_params%sizes_recv, &
                                     comm_params%start_recv, &
                                     comm_params%types_recv, &
                                     recv_buf, &
                                     comm_params%sizes_recv, &
                                     comm_params%start_recv, &
                                     comm_params%types_recv, &
                                     comm_params%comm, &
                                     ierr &
                                    )




      end if

   end subroutine Send_recv_neighbors_diag


   subroutine Send_neighbors_diag (this, comm_params, send_buf, recv_buf, request)
      class (communication_t)             , intent(in)       :: this 
      real(8), dimension(:), allocatable, intent(in out) :: send_buf
      real(8), dimension(:), allocatable, intent(in out)     :: recv_buf
      integer, intent(in out)                                :: request

      type(communication_parameters_t) , pointer, intent(in out) :: comm_params
      integer :: ierr
      integer :: my_rank
      integer :: result_len
      character (mpi_max_error_string) :: err_string


      if (this%is_parallel .eqv. .true.) then

         call MPI_Ineighbor_alltoallw(send_buf, &
                                     comm_params%sizes_recv, &
                                     comm_params%start_recv, &
                                     comm_params%types_recv, &
                                     recv_buf, &
                                     comm_params%sizes_recv, &
                                     comm_params%start_recv, &
                                     comm_params%types_recv, &
                                     comm_params%comm, &
                                     request, &
                                     ierr &
                                    )
      end if

   end subroutine Send_neighbors_diag


   subroutine Wait_recv_neighbors_diag (this, comm_params, send_buf, recv_buf, request)
      class (communication_t)             , intent(in)    :: this 
      real(8), dimension(:), allocatable, intent(in out) :: send_buf
      real(8), dimension(:), allocatable, intent(in out)     :: recv_buf
      integer, intent(in out)                                :: request

      type(communication_parameters_t) , pointer, intent(in out) :: comm_params
      integer :: ierr
      integer :: my_rank
      integer :: result_len
      character (mpi_max_error_string) :: err_string


      if (this%is_parallel .eqv. .true.) then

         call MPI_wait(request, MPI_status_ignore, ierr)


      end if
   end subroutine Wait_recv_neighbors_diag


   subroutine Send_get_max (this, send_val, recv_val)
      class (communication_t)             , intent(in)    :: this 
      real(8), intent(in)  :: send_val
      real(8), intent(in out) :: recv_val
      integer :: ierr
      integer :: rank, i, j  


      if (this%is_parallel .eqv. .true.) then

         call MPI_allreduce(send_val, recv_val, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
      end if
   end subroutine Send_get_max


   subroutine Send_get_min (this, send_val, recv_val)
      class (communication_t)             , intent(in)    :: this 
      real(8), intent(in)  :: send_val
      real(8), intent(in out) :: recv_val
      integer :: ierr
      integer :: rank, i, j  


      if (this%is_parallel .eqv. .true.) then

         call MPI_allreduce(send_val, recv_val, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
      end if
   end subroutine Send_get_min


   subroutine Send_get_array_by_min_val (this, val, array)
      class (communication_t)             , intent(in)    :: this 
      real(8), intent(in out)  :: val
      real(8), dimension(:), allocatable, intent(in out)  :: array
      real(8), dimension(:), allocatable  :: send_array_with_val, recv_array_with_val
      integer :: op_min_val, ierr

      integer :: rank, i, j  

      if (this%is_parallel .eqv. .true.) then

         call MPI_op_create(MPI_User_fn_array_by_min_val, .true., op_min_val, ierr) 

         allocate(send_array_with_val(size(array)+1))
         allocate(recv_array_with_val(size(array)+1))
         send_array_with_val(1) = val
         send_array_with_val(2:) = array
         recv_array_with_val(1) = 1d30
         recv_array_with_val(2:) = 1d30

         call MPI_allreduce(send_array_with_val, recv_array_with_val, size(array)+1,&
          MPI_DOUBLE_PRECISION, op_min_val, MPI_COMM_WORLD, ierr)


         val = recv_array_with_val(1)
         array = recv_array_with_val(2:)

         deallocate(send_array_with_val)
         deallocate(recv_array_with_val)
      end if
   end subroutine Send_get_array_by_min_val

   subroutine Send_get_array_by_max_val (this, val, array)
      class (communication_t)             , intent(in)    :: this 
      real(8), intent(in out)  :: val
      real(8), dimension(:), allocatable, intent(inout)  :: array
      real(8), dimension(:), allocatable  :: send_array_with_val, recv_array_with_val
      integer :: op_max_val, ierr
      real(8) :: temp

      integer :: rank, i, j  

      if (this%is_parallel .eqv. .true.) then

         rank = Get_rank(MPI_COMM_WORLD)

         call MPI_op_create(MPI_User_fn_array_by_max_val, .true., op_max_val, ierr) 

         allocate(send_array_with_val(size(array)+1))
         allocate(recv_array_with_val(size(array)+1))
         send_array_with_val(1) = val
         send_array_with_val(2:) = array
         recv_array_with_val(1) = -1d30
         recv_array_with_val(2:) = -1d30


         call MPI_allreduce(send_array_with_val, recv_array_with_val, size(array)+1&
         , MPI_DOUBLE_PRECISION, op_max_val, MPI_COMM_WORLD, ierr)

         val = recv_array_with_val(1)
         array = recv_array_with_val(2:)

         deallocate(send_array_with_val)
         deallocate(recv_array_with_val)

      end if
   end subroutine Send_get_array_by_max_val


   subroutine Send_get_sum (this, send_val, recv_val)
      class (communication_t)             , intent(in)    :: this 
      real(8), intent(in)  :: send_val
      real(8), intent(in out) :: recv_val
      integer :: ierr
      integer :: rank, i, j  

      if (this%is_parallel .eqv. .true.) then

         call MPI_allreduce(send_val, recv_val, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      end if
   end subroutine Send_get_sum

   subroutine Destroy (this)
      class (communication_t)       , intent(in)  :: this 

   end subroutine Destroy


   type(integer) function rank(comm)
      implicit none
      integer :: comm, ierr
      call MPI_comm_rank(comm, rank, ierr)
   end function



end module communication_module
