module communication_parameters_module
  use mpi
  implicit none
  private
  public :: communication_parameters_t
  type :: communication_parameters_t
     integer, pointer                          :: mpi_params 
     integer                                   :: communicate_to
     integer                                   :: tag
     integer                                   :: comm
     integer                                   :: d1, d2, d3
     integer, public                           :: dim_offset
     integer, dimension (:), allocatable                       :: types_send, types_recv
     integer, dimension (:), allocatable                      :: sizes_send, sizes_recv
     integer (kind=MPI_address_kind), dimension(:), allocatable :: start_send, start_recv

  end type
  interface communication_parameters_t
     module procedure Constructor
     module procedure Constructor_with_offset_and_diag
     module procedure Constructor_4d_with_offset_and_diag
  end interface communication_parameters_t
contains


  type(communication_parameters_t) function Constructor_with_offset_and_diag(d1, d2, d3, nb_num, comm, npx, npy, npz, coords, dim_offset)
     implicit none
     integer, intent(in)                          :: d1, d2, d3
     integer, intent(in)                          :: nb_num
     integer, intent(in)                          :: comm
     integer, intent(in)                          :: npx, npy, npz
     integer, dimension(3), intent(in)            :: coords
     integer, intent(in)                          :: dim_offset
     integer, dimension (3) :: global_sizes, recv_global_sizes
     integer :: ierror, ierror2, ierror3
     integer :: rear_front_t, left_right_t, up_down_t, recv_rear_front_t, recv_left_right_t, recv_up_down_t
     integer :: x_diag_t, y_diag_t, z_diag_t, recv_x_diag_t, recv_y_diag_t, recv_z_diag_t
     integer :: i, j ,k, recv_sum, x,y,z, nb
     integer (kind=MPI_address_kind) :: lb, sizeofreal8
     integer :: d_ndim
     integer, dimension(2)                       :: send_type_left_right, send_type_rear_front, send_type_down_up
     integer, dimension(2)                       :: recv_type_left_right, recv_type_rear_front, recv_type_down_up
     integer, dimension(2)                     :: send_size_left_right, send_size_rear_front, send_size_down_up
     integer, dimension(2)                     :: recv_size_left_right, recv_size_rear_front, recv_size_down_up
     integer (kind=MPI_address_kind), dimension(2) :: send_start_left_right, send_start_rear_front, send_start_down_up
     real(8) :: r8
     integer :: rank 
     call MPI_comm_rank(comm, rank, ierror) 
      r8 = 0d0

     Constructor_with_offset_and_diag%dim_offset = dim_offset
     Constructor_with_offset_and_diag%comm = comm
     Constructor_with_offset_and_diag%d1 = d1
     Constructor_with_offset_and_diag%d2 = d2
     Constructor_with_offset_and_diag%d3 = d3

     global_sizes = (/d1+2, d2+2, d3+2/)
     recv_global_sizes = (/d1+2, d2+2, d3+2/)

     call MPI_type_create_subarray(3, global_sizes, (/1, d2+2, d3+2/), (/0, 0, 0/), &
        MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, left_right_t, ierror)
     call MPI_type_create_subarray(3, global_sizes, (/d1+2, 1, d3+2/), (/0, 0, 0/), &
        MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, rear_front_t, ierror)
     call MPI_type_create_subarray(3, global_sizes, (/d1+2, d2+2, 1/), (/0, 0, 0/), &
        MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, up_down_t, ierror)

     call MPI_type_create_subarray(3, global_sizes, (/1, 1, d3+2/), (/0, 0, 0/), &
        MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, z_diag_t, ierror)
     call MPI_type_create_subarray(3, global_sizes, (/1, d2+2, 1/), (/0, 0, 0/), &
        MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, y_diag_t, ierror)
     call MPI_type_create_subarray(3, global_sizes, (/d1+2, 1, 1/), (/0, 0, 0/), &
        MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, x_diag_t, ierror)


     call MPI_type_commit(left_right_t, ierror)
     call MPI_type_commit(rear_front_t, ierror)
     call MPI_type_commit(up_down_t, ierror)

     call MPI_type_commit(x_diag_t, ierror)
     call MPI_type_commit(y_diag_t, ierror)
     call MPI_type_commit(z_diag_t, ierror)

     allocate(Constructor_with_offset_and_diag%types_send(nb_num))
     allocate(Constructor_with_offset_and_diag%types_recv(nb_num))
     allocate(Constructor_with_offset_and_diag%sizes_send(nb_num))
     allocate(Constructor_with_offset_and_diag%sizes_recv(nb_num))
     allocate(Constructor_with_offset_and_diag%start_send(nb_num))
     allocate(Constructor_with_offset_and_diag%start_recv(nb_num))

     x = coords(1)
     y = coords(2)
     z = coords(3)
     nb = 0

     if (x+1 /= npx+1) then
        nb = nb + 1
        Constructor_with_offset_and_diag%types_send(nb) = left_right_t       
        Constructor_with_offset_and_diag%types_recv(nb) = MPI_DOUBLE_PRECISION
        Constructor_with_offset_and_diag%sizes_send(nb) = 1
        Constructor_with_offset_and_diag%sizes_recv(nb) = (d2+2)*(d3+2)
        Constructor_with_offset_and_diag%start_send(nb) = (d1 - dim_offset) * sizeof(r8)
        Constructor_with_offset_and_diag%start_recv(nb) = (0) * sizeof(r8)
     end if
     if (x-1 /= 0) then
        nb = nb + 1
        Constructor_with_offset_and_diag%types_send(nb) = left_right_t       
        Constructor_with_offset_and_diag%types_recv(nb) = MPI_DOUBLE_PRECISION
        Constructor_with_offset_and_diag%sizes_send(nb) = 1
        Constructor_with_offset_and_diag%sizes_recv(nb) = (d2+2)*(d3+2)
        Constructor_with_offset_and_diag%start_send(nb) = (1 + dim_offset) * sizeof(r8)
        Constructor_with_offset_and_diag%start_recv(nb) = ((d2+2)*(d3+2)) * sizeof(r8)
    end if
      if (y+1 /= npy+1) then
         nb = nb + 1
         Constructor_with_offset_and_diag%types_send(nb) = rear_front_t        
         Constructor_with_offset_and_diag%types_recv(nb) = MPI_DOUBLE_PRECISION
         Constructor_with_offset_and_diag%sizes_send(nb) = 1
         Constructor_with_offset_and_diag%sizes_recv(nb) = (d1+2)*(d3+2)
         Constructor_with_offset_and_diag%start_send(nb) = ((d2 - dim_offset) * (d1+2))  * sizeof(r8)
         Constructor_with_offset_and_diag%start_recv(nb) = (2*(d2+2)*(d3+2)) * sizeof(r8)
      end if
      if (y-1 /= 0) then
         nb = nb + 1
         Constructor_with_offset_and_diag%types_send(nb) = rear_front_t        
         Constructor_with_offset_and_diag%types_recv(nb) = MPI_DOUBLE_PRECISION
         Constructor_with_offset_and_diag%sizes_send(nb) = 1
         Constructor_with_offset_and_diag%sizes_recv(nb) = (d1+2)*(d3+2)
         Constructor_with_offset_and_diag%start_send(nb) = ((1 + dim_offset) * (d1+2)) * sizeof(r8)
         Constructor_with_offset_and_diag%start_recv(nb) = (2*(d2+2)*(d3+2)+(d1+2)*(d3+2)) * sizeof(r8)
      end if
      if (z+1 /= npz+1) then
         nb = nb + 1
        Constructor_with_offset_and_diag%types_send(nb) = up_down_t           
        Constructor_with_offset_and_diag%types_recv(nb) = MPI_DOUBLE_PRECISION
        Constructor_with_offset_and_diag%sizes_send(nb) = 1
        Constructor_with_offset_and_diag%sizes_recv(nb) = (d1+2)*(d2+2)
        Constructor_with_offset_and_diag%start_send(nb) = ((d3 - dim_offset) * (d2+2) * (d1+2)) * sizeof(r8)
        Constructor_with_offset_and_diag%start_recv(nb) = (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)) * sizeof(r8)
      end if
      if (z-1 /= 0) then
         nb = nb + 1
        Constructor_with_offset_and_diag%types_send(nb) = up_down_t           
        Constructor_with_offset_and_diag%types_recv(nb) = MPI_DOUBLE_PRECISION
        Constructor_with_offset_and_diag%sizes_send(nb) = 1
        Constructor_with_offset_and_diag%sizes_recv(nb) = (d1+2)*(d2+2)
        Constructor_with_offset_and_diag%start_send(nb) = ((1 + dim_offset) * (d1+2) * (d2+2)) * sizeof(r8)
        Constructor_with_offset_and_diag%start_recv(nb) = (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+(d1+2)*(d2+2)) * sizeof(r8)
      end if

      if ((x+1 /= npx+1) .and. (y+1 /= npy+1)) then
         nb = nb + 1
         Constructor_with_offset_and_diag%types_send(nb) = z_diag_t       
         Constructor_with_offset_and_diag%types_recv(nb) = MPI_DOUBLE_PRECISION
         Constructor_with_offset_and_diag%sizes_send(nb) = 1
         Constructor_with_offset_and_diag%sizes_recv(nb) = (d3+2)
         Constructor_with_offset_and_diag%start_send(nb) = ((d2 - dim_offset)*(d1+2) + (d1 - dim_offset)) * sizeof(r8)
         Constructor_with_offset_and_diag%start_recv(nb) = (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)) * sizeof(r8)

      end if
      if ((x+1 /= npx+1) .and. (y-1 /= 0)) then
         nb = nb + 1
         Constructor_with_offset_and_diag%types_send(nb) = z_diag_t       
         Constructor_with_offset_and_diag%types_recv(nb) = MPI_DOUBLE_PRECISION
         Constructor_with_offset_and_diag%sizes_send(nb) = 1
         Constructor_with_offset_and_diag%sizes_recv(nb) = (d3+2)
         Constructor_with_offset_and_diag%start_send(nb) = ((1 + dim_offset)*(d1+2) + (d1 - dim_offset)) * sizeof(r8)
         Constructor_with_offset_and_diag%start_recv(nb) = (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+(d3+2)) * sizeof(r8)

      end if
      if ((x-1 /= 0) .and. (y+1 /= npy+1)) then
         nb = nb + 1
         Constructor_with_offset_and_diag%types_send(nb)  = z_diag_t       
         Constructor_with_offset_and_diag%types_recv(nb)  = MPI_DOUBLE_PRECISION
         Constructor_with_offset_and_diag%sizes_send(nb)  = 1
         Constructor_with_offset_and_diag%sizes_recv(nb)  = (d3+2)
         Constructor_with_offset_and_diag%start_send(nb)  = ((d2 - dim_offset) * (d1+2) + (1 + dim_offset)) * sizeof(r8)
         Constructor_with_offset_and_diag%start_recv(nb)  = (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+2*(d3+2)) * sizeof(r8)

      end if
      if ((x-1 /= 0) .and. (y-1 /= 0)) then
         nb = nb + 1
         Constructor_with_offset_and_diag%types_send(nb) = z_diag_t       
         Constructor_with_offset_and_diag%types_recv(nb) = MPI_DOUBLE_PRECISION
         Constructor_with_offset_and_diag%sizes_send(nb) = 1
         Constructor_with_offset_and_diag%sizes_recv(nb) = (d3+2)
         Constructor_with_offset_and_diag%start_send(nb) = ((1 + dim_offset) * (d1+2) + (1 + dim_offset)) * sizeof(r8)
         Constructor_with_offset_and_diag%start_recv(nb) = (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+3*(d3+2)) * sizeof(r8)
      end if
      if ((x+1 /= npx+1) .and. (z+1 /= npz+1)) then
         nb = nb + 1
         Constructor_with_offset_and_diag%types_send(nb) = y_diag_t       
         Constructor_with_offset_and_diag%types_recv(nb) = MPI_DOUBLE_PRECISION
         Constructor_with_offset_and_diag%sizes_send(nb) = 1
         Constructor_with_offset_and_diag%sizes_recv(nb) = (d2+2)
         Constructor_with_offset_and_diag%start_send(nb) = ((d3 - dim_offset) * (d2+2) * (d1+2) + (d1 - dim_offset)) * sizeof(r8)
         Constructor_with_offset_and_diag%start_recv(nb) = (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)) * sizeof(r8)
      end if
      if ((x+1 /= npx+1) .and. (z-1 /= 0)) then
         nb = nb + 1
         Constructor_with_offset_and_diag%types_send(nb) = y_diag_t       
         Constructor_with_offset_and_diag%types_recv(nb) = MPI_DOUBLE_PRECISION
         Constructor_with_offset_and_diag%sizes_send(nb) = 1
         Constructor_with_offset_and_diag%sizes_recv(nb) = (d2+2)
         Constructor_with_offset_and_diag%start_send(nb) = ((1 + dim_offset) * (d2+2) * (d1+2) + (d1 - dim_offset)) * sizeof(r8)
         Constructor_with_offset_and_diag%start_recv(nb) = (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+(d2+2))&
                                                            * sizeof(r8)
      end if
      if ((x-1 /= 0) .and. (z+1 /= npz+1)) then
         nb = nb + 1
         Constructor_with_offset_and_diag%types_send(nb) = y_diag_t       
         Constructor_with_offset_and_diag%types_recv(nb) = MPI_DOUBLE_PRECISION
         Constructor_with_offset_and_diag%sizes_send(nb) = 1
         Constructor_with_offset_and_diag%sizes_recv(nb) = (d2+2)
         Constructor_with_offset_and_diag%start_send(nb) = ((d3 - dim_offset) * (d2+2) * (d1+2) + (1 + dim_offset)) * sizeof(r8)
         Constructor_with_offset_and_diag%start_recv(nb) = (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+2*(d2+2))&
                                                           * sizeof(r8)
      end if
      if ((x-1 /= 0) .and. (z-1 /= 0)) then
         nb = nb + 1
         Constructor_with_offset_and_diag%types_send(nb) = y_diag_t       
         Constructor_with_offset_and_diag%types_recv(nb) = MPI_DOUBLE_PRECISION
         Constructor_with_offset_and_diag%sizes_send(nb) = 1
         Constructor_with_offset_and_diag%sizes_recv(nb) = (d2+2)
         Constructor_with_offset_and_diag%start_send(nb) = ((1 + dim_offset) * (d2+2) * (d1+2) + (1 + dim_offset)) * sizeof(r8)
         Constructor_with_offset_and_diag%start_recv(nb) = (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)&
                                                           +4*(d3+2)+3*(d2+2)) * sizeof(r8)
      end if
      if ((y+1 /= npy+1) .and. (z+1 /= npz+1)) then
         nb = nb + 1
         Constructor_with_offset_and_diag%types_send(nb) = x_diag_t       
         Constructor_with_offset_and_diag%types_recv(nb) = MPI_DOUBLE_PRECISION
         Constructor_with_offset_and_diag%sizes_send(nb) = 1
         Constructor_with_offset_and_diag%sizes_recv(nb) = (d1+2)
         Constructor_with_offset_and_diag%start_send(nb) = ((d3 - dim_offset) * (d2+2) * (d1+2)&
                                                           + (d2 - dim_offset) * (d1+2)) * sizeof(r8)
         Constructor_with_offset_and_diag%start_recv(nb) = (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)&
                                                           +4*(d3+2)+4*(d2+2)) * sizeof(r8)
      end if
      if ((y+1 /= npy+1) .and. (z-1 /= 0)) then
         nb = nb + 1
         Constructor_with_offset_and_diag%types_send(nb) = x_diag_t       
         Constructor_with_offset_and_diag%types_recv(nb) = MPI_DOUBLE_PRECISION
         Constructor_with_offset_and_diag%sizes_send(nb) = 1
         Constructor_with_offset_and_diag%sizes_recv(nb) = (d1+2)
         Constructor_with_offset_and_diag%start_send(nb) = ((1 + dim_offset) * (d2+2) * (d1+2)&
                                                          + (d2 - dim_offset) * (d1+2)) * sizeof(r8)
         Constructor_with_offset_and_diag%start_recv(nb) = (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)&
                                                           +4*(d3+2)+4*(d2+2)+(d1+2)) * sizeof(r8)
      end if
      if ((y-1 /= 0) .and. (z+1 /= npz+1)) then
         nb = nb + 1
         Constructor_with_offset_and_diag%types_send(nb) = x_diag_t       
         Constructor_with_offset_and_diag%types_recv(nb) = MPI_DOUBLE_PRECISION
         Constructor_with_offset_and_diag%sizes_send(nb) = 1
         Constructor_with_offset_and_diag%sizes_recv(nb) = (d1+2)
         Constructor_with_offset_and_diag%start_send(nb) = ((d3 - dim_offset) * (d2+2) * (d1+2) +&
                                                           (1 + dim_offset) * (d1+2)) * sizeof(r8)
         Constructor_with_offset_and_diag%start_recv(nb) = (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+&
                                                           4*(d3+2)+4*(d2+2)+2*(d1+2)) * sizeof(r8)
      end if
      if ((y-1 /= 0) .and. (z-1 /= 0)) then
         nb = nb + 1
         Constructor_with_offset_and_diag%types_send(nb) = x_diag_t       
         Constructor_with_offset_and_diag%types_recv(nb) = MPI_DOUBLE_PRECISION
         Constructor_with_offset_and_diag%sizes_send(nb) = 1
         Constructor_with_offset_and_diag%sizes_recv(nb) = (d1+2)
         Constructor_with_offset_and_diag%start_send(nb) = ((1 + dim_offset) * (d2+2) * (d1+2) +&
                                                           (1 + dim_offset) * (d1+2)) * sizeof(r8)
         Constructor_with_offset_and_diag%start_recv(nb) = (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)&
                                                            +4*(d2+2)+3*(d1+2)) * sizeof(r8)
      end if



      if ((x+1 /= npx+1) .and. (y+1 /= npy+1) .and. (z+1 /= npz+1)) then
         nb = nb + 1
         Constructor_with_offset_and_diag%types_send(nb) = MPI_DOUBLE_PRECISION
         Constructor_with_offset_and_diag%types_recv(nb) = MPI_DOUBLE_PRECISION
         Constructor_with_offset_and_diag%sizes_send(nb) = 1
         Constructor_with_offset_and_diag%sizes_recv(nb) = 1
         Constructor_with_offset_and_diag%start_send(nb) = ((d3 - dim_offset) * (d2+2) * (d1+2) + (d2 - dim_offset)&
                                                           * (d1+2) + (d1 - dim_offset)) * sizeof(r8)
         Constructor_with_offset_and_diag%start_recv(nb) = (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)&
                                                            +4*(d2+2)+4*(d1+2)) * sizeof(r8)
      end if
      if ((x+1 /= npx+1) .and. (y+1 /= npy+1) .and. (z-1 /= 0)) then
         nb = nb + 1
         Constructor_with_offset_and_diag%types_send(nb) = MPI_DOUBLE_PRECISION
         Constructor_with_offset_and_diag%types_recv(nb) = MPI_DOUBLE_PRECISION
         Constructor_with_offset_and_diag%sizes_send(nb) = 1
         Constructor_with_offset_and_diag%sizes_recv(nb) = 1
         Constructor_with_offset_and_diag%start_send(nb) = ((1 + dim_offset) * (d2+2) * (d1+2) + (d2 - dim_offset)&
                                                           * (d1+2) + (d1 - dim_offset)) * sizeof(r8)
         Constructor_with_offset_and_diag%start_recv(nb) = (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+&
                                                            4*(d2+2)+4*(d1+2)+1) * sizeof(r8)
      end if

      if ((x+1 /= npx+1) .and. (y-1 /= 0) .and. (z+1 /= npz+1)) then
         nb = nb + 1
         Constructor_with_offset_and_diag%types_send(nb) = MPI_DOUBLE_PRECISION
         Constructor_with_offset_and_diag%types_recv(nb) = MPI_DOUBLE_PRECISION
         Constructor_with_offset_and_diag%sizes_send(nb) = 1
         Constructor_with_offset_and_diag%sizes_recv(nb) = 1
         Constructor_with_offset_and_diag%start_send(nb) = ((d3 - dim_offset) * (d2+2) * (d1+2) + (1 + dim_offset)&
                                                           * (d1+2) + (d1 - dim_offset)) * sizeof(r8)
         Constructor_with_offset_and_diag%start_recv(nb) = (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+&
                                                           4*(d3+2)+4*(d2+2)+4*(d1+2)+2) * sizeof(r8)
      end if
      if ((x+1 /= npx+1) .and. (y-1 /= 0) .and. (z-1 /= 0)) then
         nb = nb + 1
         Constructor_with_offset_and_diag%types_send(nb) = MPI_DOUBLE_PRECISION
         Constructor_with_offset_and_diag%types_recv(nb) = MPI_DOUBLE_PRECISION
         Constructor_with_offset_and_diag%sizes_send(nb) = 1
         Constructor_with_offset_and_diag%sizes_recv(nb) = 1
         Constructor_with_offset_and_diag%start_send(nb) = ((1 + dim_offset) * (d2+2) * (d1+2) + (1 + dim_offset) &
                                                           * (d1+2) + (d1 - dim_offset)) * sizeof(r8)
         Constructor_with_offset_and_diag%start_recv(nb) = (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)&
                                                            +4*(d2+2)+4*(d1+2)+3) * sizeof(r8)
      end if

      if ((x-1 /= 0) .and. (y+1 /= npy+1) .and. (z+1 /= npz+1)) then
         nb = nb + 1
         Constructor_with_offset_and_diag%types_send(nb) = MPI_DOUBLE_PRECISION
         Constructor_with_offset_and_diag%types_recv(nb) = MPI_DOUBLE_PRECISION
         Constructor_with_offset_and_diag%sizes_send(nb) = 1
         Constructor_with_offset_and_diag%sizes_recv(nb) = 1
         Constructor_with_offset_and_diag%start_send(nb) = ((d3 - dim_offset) * (d2+2) * (d1+2) + (d2 - dim_offset)&
                                                           * (d1+2) + (1 + dim_offset)) * sizeof(r8)
         Constructor_with_offset_and_diag%start_recv(nb) = (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+&
                                                            4*(d2+2)+4*(d1+2)+4) * sizeof(r8)
      end if
      if ((x-1 /= 0) .and. (y+1 /= npy+1) .and. (z-1 /= 0)) then
         nb = nb + 1
         Constructor_with_offset_and_diag%types_send(nb) = MPI_DOUBLE_PRECISION
         Constructor_with_offset_and_diag%types_recv(nb) = MPI_DOUBLE_PRECISION
         Constructor_with_offset_and_diag%sizes_send(nb) = 1
         Constructor_with_offset_and_diag%sizes_recv(nb) = 1
         Constructor_with_offset_and_diag%start_send(nb) = ((1 + dim_offset) * (d2+2) * (d1+2) + (d2 - dim_offset)&
                                                           * (d1+2) + (1 + dim_offset)) * sizeof(r8)
         Constructor_with_offset_and_diag%start_recv(nb) = (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+&
                                                           4*(d2+2)+4*(d1+2)+5) * sizeof(r8)
      end if


      if ((x-1 /= 0) .and. (y-1 /= 0) .and. (z+1 /= npz+1)) then
         nb = nb + 1
         Constructor_with_offset_and_diag%types_send(nb) = MPI_DOUBLE_PRECISION
         Constructor_with_offset_and_diag%types_recv(nb) = MPI_DOUBLE_PRECISION
         Constructor_with_offset_and_diag%sizes_send(nb) = 1
         Constructor_with_offset_and_diag%sizes_recv(nb) = 1
         Constructor_with_offset_and_diag%start_send(nb) = ((d3 - dim_offset) * (d2+2) * (d1+2) + (1 + dim_offset) &
                                                           * (d1+2) + (1 + dim_offset)) * sizeof(r8)
         Constructor_with_offset_and_diag%start_recv(nb) = (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+&
                                                            4*(d3+2)+4*(d2+2)+4*(d1+2)+6) * sizeof(r8)
      end if
      if ((x-1 /= 0) .and. (y-1 /= 0) .and. (z-1 /= 0)) then
         nb = nb + 1
         Constructor_with_offset_and_diag%types_send(nb) = MPI_DOUBLE_PRECISION
         Constructor_with_offset_and_diag%types_recv(nb) = MPI_DOUBLE_PRECISION
         Constructor_with_offset_and_diag%sizes_send(nb) = 1
         Constructor_with_offset_and_diag%sizes_recv(nb) = 1
         Constructor_with_offset_and_diag%start_send(nb) = ((1 + dim_offset) * (d2+2) * (d1+2) + (1 + dim_offset)&
                                                           * (d1+2) + (1 + dim_offset)) * sizeof(r8)
         Constructor_with_offset_and_diag%start_recv(nb) = (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+&
                                                            4*(d2+2)+4*(d1+2)+7) * sizeof(r8)
      end if



  end function


  type(communication_parameters_t) function Constructor_4d_with_offset_and_diag(nmats, d1, d2, d3, nb_num, comm, npx, npy, npz, coords, dim_offset)
     implicit none
     integer, intent(in)                          :: nmats
     integer, intent(in)                          :: d1, d2, d3
     integer, intent(in)                          :: nb_num
     integer, intent(in)                          :: comm
     integer, intent(in)                          :: npx, npy, npz
     integer, dimension(3), intent(in)            :: coords
     integer, intent(in)                          :: dim_offset
     integer, dimension (4) :: global_sizes, recv_global_sizes
     integer :: ierror, ierror2, ierror3
     integer :: rear_front_t, left_right_t, up_down_t, recv_rear_front_t, recv_left_right_t, recv_up_down_t
     integer :: x_diag_t, y_diag_t, z_diag_t, recv_x_diag_t, recv_y_diag_t, recv_z_diag_t
     integer :: mats_t
     integer :: i, j ,k, recv_sum, x,y,z, nb
     integer (kind=MPI_address_kind) :: lb, sizeofreal8
     integer :: d_ndim
     integer, dimension(2)                       :: send_type_left_right, send_type_rear_front, send_type_down_up
     integer, dimension(2)                       :: recv_type_left_right, recv_type_rear_front, recv_type_down_up
     integer, dimension(2)                     :: send_size_left_right, send_size_rear_front, send_size_down_up
     integer, dimension(2)                     :: recv_size_left_right, recv_size_rear_front, recv_size_down_up
     integer (kind=MPI_address_kind), dimension(2) :: send_start_left_right, send_start_rear_front, send_start_down_up
     real(8) :: r8
     integer :: rank
     call MPI_comm_rank(comm, rank, ierror)
      r8 = 0

     Constructor_4d_with_offset_and_diag%dim_offset = dim_offset
     Constructor_4d_with_offset_and_diag%comm = comm
     Constructor_4d_with_offset_and_diag%d1 = d1
     Constructor_4d_with_offset_and_diag%d2 = d2
     Constructor_4d_with_offset_and_diag%d3 = d3

     global_sizes = (/nmats, d1+2, d2+2, d3+2/)
     recv_global_sizes = (/nmats, d1+2, d2+2, d3+2/)

     call MPI_type_create_subarray(4, global_sizes, (/nmats, 1, d2+2, d3+2/), (/0, 0, 0, 0/), &
        MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, left_right_t, ierror)
     call MPI_type_create_subarray(4, global_sizes, (/nmats, d1+2, 1, d3+2/), (/0, 0, 0, 0/), &
        MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, rear_front_t, ierror)
     call MPI_type_create_subarray(4, global_sizes, (/nmats, d1+2, d2+2, 1/), (/0, 0, 0, 0/), &
        MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, up_down_t, ierror)
     call MPI_type_create_subarray(4, global_sizes, (/nmats, 1, 1, d3+2/), (/0, 0, 0, 0/), &
        MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, z_diag_t, ierror)
     call MPI_type_create_subarray(4, global_sizes, (/nmats, 1, d2+2, 1/), (/0, 0, 0, 0/), &
        MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, y_diag_t, ierror)
     call MPI_type_create_subarray(4, global_sizes, (/nmats, d1+2, 1, 1/), (/0, 0, 0, 0/), &
        MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, x_diag_t, ierror)

     call MPI_type_create_subarray(4, global_sizes, (/nmats, 1, 1, 1/), (/0, 0, 0, 0/), &
        MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, mats_t, ierror)

     call MPI_type_commit(left_right_t, ierror)
     call MPI_type_commit(rear_front_t, ierror)
     call MPI_type_commit(up_down_t, ierror)

     call MPI_type_commit(x_diag_t, ierror)
     call MPI_type_commit(y_diag_t, ierror)
     call MPI_type_commit(z_diag_t, ierror)

     call MPI_type_commit(mats_t, ierror)

     allocate(Constructor_4d_with_offset_and_diag%types_send(nb_num))
     allocate(Constructor_4d_with_offset_and_diag%types_recv(nb_num))
     allocate(Constructor_4d_with_offset_and_diag%sizes_send(nb_num))
     allocate(Constructor_4d_with_offset_and_diag%sizes_recv(nb_num))
     allocate(Constructor_4d_with_offset_and_diag%start_send(nb_num))
     allocate(Constructor_4d_with_offset_and_diag%start_recv(nb_num))

     x = coords(1)
     y = coords(2)
     z = coords(3)
     nb = 0


     if (x+1 /= npx+1) then
        nb = nb + 1
        Constructor_4d_with_offset_and_diag%types_send(nb) = left_right_t
        Constructor_4d_with_offset_and_diag%types_recv(nb) = MPI_DOUBLE_PRECISION
        Constructor_4d_with_offset_and_diag%sizes_send(nb) = 1
        Constructor_4d_with_offset_and_diag%sizes_recv(nb) = nmats * (d2+2)*(d3+2)
        Constructor_4d_with_offset_and_diag%start_send(nb) = (d1 - dim_offset) * nmats * sizeof(r8)
        Constructor_4d_with_offset_and_diag%start_recv(nb) = (0) * sizeof(r8)
     end if

     if (x-1 /= 0) then
        nb = nb + 1
        Constructor_4d_with_offset_and_diag%types_send(nb) = left_right_t
        Constructor_4d_with_offset_and_diag%types_recv(nb) = MPI_DOUBLE_PRECISION
        Constructor_4d_with_offset_and_diag%sizes_send(nb) = 1
        Constructor_4d_with_offset_and_diag%sizes_recv(nb) = nmats * (d2+2)*(d3+2)
        Constructor_4d_with_offset_and_diag%start_send(nb) = (1 + dim_offset) * nmats * sizeof(r8)
        Constructor_4d_with_offset_and_diag%start_recv(nb) = ((d2+2)*(d3+2)) * nmats * sizeof(r8)
    end if

      if (y+1 /= npy+1) then
         nb = nb + 1
         Constructor_4d_with_offset_and_diag%types_send(nb) = rear_front_t
         Constructor_4d_with_offset_and_diag%types_recv(nb) = MPI_DOUBLE_PRECISION
         Constructor_4d_with_offset_and_diag%sizes_send(nb) = 1
         Constructor_4d_with_offset_and_diag%sizes_recv(nb) = nmats * (d1+2)*(d3+2)
         Constructor_4d_with_offset_and_diag%start_send(nb) = (d2 - dim_offset) * (d1+2) * nmats * sizeof(r8)
         Constructor_4d_with_offset_and_diag%start_recv(nb) = (2*(d2+2)*(d3+2)) * nmats * sizeof(r8)
      end if

      if (y-1 /= 0) then
         nb = nb + 1
         Constructor_4d_with_offset_and_diag%types_send(nb) = rear_front_t
         Constructor_4d_with_offset_and_diag%types_recv(nb) = MPI_DOUBLE_PRECISION
         Constructor_4d_with_offset_and_diag%sizes_send(nb) = 1
         Constructor_4d_with_offset_and_diag%sizes_recv(nb) = nmats * (d1+2)*(d3+2)
         Constructor_4d_with_offset_and_diag%start_send(nb) = (1 + dim_offset) * (d1+2) * nmats * sizeof(r8)
         Constructor_4d_with_offset_and_diag%start_recv(nb) = (2*(d2+2)*(d3+2)+(d1+2)*(d3+2)) * nmats * sizeof(r8)
      end if

      if (z+1 /= npz+1) then
         nb = nb + 1
        Constructor_4d_with_offset_and_diag%types_send(nb) = up_down_t
        Constructor_4d_with_offset_and_diag%types_recv(nb) = MPI_DOUBLE_PRECISION
        Constructor_4d_with_offset_and_diag%sizes_send(nb) = 1
        Constructor_4d_with_offset_and_diag%sizes_recv(nb) = nmats * (d1+2)*(d2+2)
        Constructor_4d_with_offset_and_diag%start_send(nb) = (d3 - dim_offset) * (d2+2) * (d1+2) * nmats * sizeof(r8)
        Constructor_4d_with_offset_and_diag%start_recv(nb) = (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)) * nmats * sizeof(r8)
      end if

      if (z-1 /= 0) then
         nb = nb + 1
        Constructor_4d_with_offset_and_diag%types_send(nb) = up_down_t
        Constructor_4d_with_offset_and_diag%types_recv(nb) = MPI_DOUBLE_PRECISION
        Constructor_4d_with_offset_and_diag%sizes_send(nb) = 1
        Constructor_4d_with_offset_and_diag%sizes_recv(nb) = nmats * (d1+2)*(d2+2)
        Constructor_4d_with_offset_and_diag%start_send(nb) = (1 + dim_offset) * (d1+2) * (d2+2) * nmats * sizeof(r8)
        Constructor_4d_with_offset_and_diag%start_recv(nb) = (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+(d1+2)*(d2+2)) * nmats * sizeof(r8)
      end if

      if ((x+1 /= npx+1) .and. (y+1 /= npy+1)) then
         nb = nb + 1
         Constructor_4d_with_offset_and_diag%types_send(nb) = z_diag_t
         Constructor_4d_with_offset_and_diag%types_recv(nb) = MPI_DOUBLE_PRECISION
         Constructor_4d_with_offset_and_diag%sizes_send(nb) = 1
         Constructor_4d_with_offset_and_diag%sizes_recv(nb) = nmats * (d3+2)
         Constructor_4d_with_offset_and_diag%start_send(nb) = ((d2 - dim_offset)*(d1+2) + (d1 - dim_offset)) * nmats * sizeof(r8)
         Constructor_4d_with_offset_and_diag%start_recv(nb) = (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)) * nmats * sizeof(r8)

      end if
      if ((x+1 /= npx+1) .and. (y-1 /= 0)) then
         nb = nb + 1
         Constructor_4d_with_offset_and_diag%types_send(nb) = z_diag_t
         Constructor_4d_with_offset_and_diag%types_recv(nb) = MPI_DOUBLE_PRECISION
         Constructor_4d_with_offset_and_diag%sizes_send(nb) = 1
         Constructor_4d_with_offset_and_diag%sizes_recv(nb) = nmats * (d3+2)
         Constructor_4d_with_offset_and_diag%start_send(nb) = ((1 + dim_offset)*(d1+2) + (d1 - dim_offset)) * nmats * sizeof(r8)
         Constructor_4d_with_offset_and_diag%start_recv(nb) = (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+(d3+2))* nmats * sizeof(r8)

      end if
      if ((x-1 /= 0) .and. (y+1 /= npy+1)) then
         nb = nb + 1
         Constructor_4d_with_offset_and_diag%types_send(nb)  = z_diag_t
         Constructor_4d_with_offset_and_diag%types_recv(nb)  = MPI_DOUBLE_PRECISION
         Constructor_4d_with_offset_and_diag%sizes_send(nb)  = 1
         Constructor_4d_with_offset_and_diag%sizes_recv(nb)  = nmats * (d3+2)
         Constructor_4d_with_offset_and_diag%start_send(nb)  = ((d2 - dim_offset) * (d1+2) + (1 + dim_offset)) * nmats * sizeof(r8)
         Constructor_4d_with_offset_and_diag%start_recv(nb)  = (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+2*(d3+2))* nmats * sizeof(r8)

      end if
      if ((x-1 /= 0) .and. (y-1 /= 0)) then
         nb = nb + 1
         Constructor_4d_with_offset_and_diag%types_send(nb) = z_diag_t
         Constructor_4d_with_offset_and_diag%types_recv(nb) = MPI_DOUBLE_PRECISION
         Constructor_4d_with_offset_and_diag%sizes_send(nb) = 1
         Constructor_4d_with_offset_and_diag%sizes_recv(nb) = nmats * (d3+2)
         Constructor_4d_with_offset_and_diag%start_send(nb) = ((1 + dim_offset) * (d1+2) + (1 + dim_offset)) * nmats * sizeof(r8)
         Constructor_4d_with_offset_and_diag%start_recv(nb) = (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+3*(d3+2))* nmats * sizeof(r8)
      end if
      if ((x+1 /= npx+1) .and. (z+1 /= npz+1)) then
         nb = nb + 1
         Constructor_4d_with_offset_and_diag%types_send(nb) = y_diag_t
         Constructor_4d_with_offset_and_diag%types_recv(nb) = MPI_DOUBLE_PRECISION
         Constructor_4d_with_offset_and_diag%sizes_send(nb) = 1
         Constructor_4d_with_offset_and_diag%sizes_recv(nb) = nmats * (d2+2)
         Constructor_4d_with_offset_and_diag%start_send(nb) = ((d3 - dim_offset) * (d2+2) * (d1+2) + (d1 - dim_offset))* nmats * sizeof(r8)
         Constructor_4d_with_offset_and_diag%start_recv(nb) = (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2))* nmats * sizeof(r8)
      end if
      if ((x+1 /= npx+1) .and. (z-1 /= 0)) then
         nb = nb + 1
         Constructor_4d_with_offset_and_diag%types_send(nb) = y_diag_t
         Constructor_4d_with_offset_and_diag%types_recv(nb) = MPI_DOUBLE_PRECISION
         Constructor_4d_with_offset_and_diag%sizes_send(nb) = 1
         Constructor_4d_with_offset_and_diag%sizes_recv(nb) = nmats * (d2+2)
         Constructor_4d_with_offset_and_diag%start_send(nb) = ((1 + dim_offset) * (d2+2) * (d1+2) + (d1 - dim_offset))* nmats * sizeof(r8)
         Constructor_4d_with_offset_and_diag%start_recv(nb) = (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+(d2+2))&
                                                           * nmats * sizeof(r8)
      end if
      if ((x-1 /= 0) .and. (z+1 /= npz+1)) then
         nb = nb + 1
         Constructor_4d_with_offset_and_diag%types_send(nb) = y_diag_t
         Constructor_4d_with_offset_and_diag%types_recv(nb) = MPI_DOUBLE_PRECISION
         Constructor_4d_with_offset_and_diag%sizes_send(nb) = 1
         Constructor_4d_with_offset_and_diag%sizes_recv(nb) = nmats * (d2+2)
         Constructor_4d_with_offset_and_diag%start_send(nb) = ((d3 - dim_offset) * (d2+2) * (d1+2) + (1 + dim_offset))* nmats * sizeof(r8)
         Constructor_4d_with_offset_and_diag%start_recv(nb) = (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+2*(d2+2))&
                                                          * nmats * sizeof(r8)
      end if
      if ((x-1 /= 0) .and. (z-1 /= 0)) then
         nb = nb + 1
         Constructor_4d_with_offset_and_diag%types_send(nb) = y_diag_t
         Constructor_4d_with_offset_and_diag%types_recv(nb) = MPI_DOUBLE_PRECISION
         Constructor_4d_with_offset_and_diag%sizes_send(nb) = 1
         Constructor_4d_with_offset_and_diag%sizes_recv(nb) = nmats * (d2+2)
         Constructor_4d_with_offset_and_diag%start_send(nb) = ((1 + dim_offset) * (d2+2) * (d1+2) + (1 + dim_offset))* nmats * sizeof(r8)
         Constructor_4d_with_offset_and_diag%start_recv(nb) = (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)&
                                                           +4*(d3+2)+3*(d2+2))* nmats * sizeof(r8)
      end if
      if ((y+1 /= npy+1) .and. (z+1 /= npz+1)) then
         nb = nb + 1
         Constructor_4d_with_offset_and_diag%types_send(nb) = x_diag_t
         Constructor_4d_with_offset_and_diag%types_recv(nb) = MPI_DOUBLE_PRECISION
         Constructor_4d_with_offset_and_diag%sizes_send(nb) = 1
         Constructor_4d_with_offset_and_diag%sizes_recv(nb) = nmats * (d1+2)
         Constructor_4d_with_offset_and_diag%start_send(nb) = ((d3 - dim_offset) * (d2+2) * (d1+2)&
                                                           + (d2 - dim_offset) * (d1+2))* nmats * sizeof(r8)
         Constructor_4d_with_offset_and_diag%start_recv(nb) = (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)&
                                                           +4*(d3+2)+4*(d2+2))* nmats * sizeof(r8)
      end if
      if ((y+1 /= npy+1) .and. (z-1 /= 0)) then
         nb = nb + 1
         Constructor_4d_with_offset_and_diag%types_send(nb) = x_diag_t
         Constructor_4d_with_offset_and_diag%types_recv(nb) = MPI_DOUBLE_PRECISION
         Constructor_4d_with_offset_and_diag%sizes_send(nb) = 1
         Constructor_4d_with_offset_and_diag%sizes_recv(nb) = nmats * (d1+2)
         Constructor_4d_with_offset_and_diag%start_send(nb) = ((1 + dim_offset) * (d2+2) * (d1+2)&
                                                          + (d2 - dim_offset) * (d1+2))* nmats * sizeof(r8)
         Constructor_4d_with_offset_and_diag%start_recv(nb) = (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)&
                                                           +4*(d3+2)+4*(d2+2)+(d1+2))* nmats * sizeof(r8)
      end if
      if ((y-1 /= 0) .and. (z+1 /= npz+1)) then
         nb = nb + 1
         Constructor_4d_with_offset_and_diag%types_send(nb) = x_diag_t
         Constructor_4d_with_offset_and_diag%types_recv(nb) = MPI_DOUBLE_PRECISION
         Constructor_4d_with_offset_and_diag%sizes_send(nb) = 1
         Constructor_4d_with_offset_and_diag%sizes_recv(nb) = nmats * (d1+2)
         Constructor_4d_with_offset_and_diag%start_send(nb) = ((d3 - dim_offset) * (d2+2) * (d1+2) +&
                                                           (1 + dim_offset) * (d1+2))* nmats * sizeof(r8)
         Constructor_4d_with_offset_and_diag%start_recv(nb) = (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+&
                                                           4*(d3+2)+4*(d2+2)+2*(d1+2))* nmats * sizeof(r8)
      end if
      if ((y-1 /= 0) .and. (z-1 /= 0)) then
         nb = nb + 1
         Constructor_4d_with_offset_and_diag%types_send(nb) = x_diag_t
         Constructor_4d_with_offset_and_diag%types_recv(nb) = MPI_DOUBLE_PRECISION
         Constructor_4d_with_offset_and_diag%sizes_send(nb) = 1
         Constructor_4d_with_offset_and_diag%sizes_recv(nb) = nmats * (d1+2)
         Constructor_4d_with_offset_and_diag%start_send(nb) = ((1 + dim_offset) * (d2+2) * (d1+2) +&
                                                           (1 + dim_offset) * (d1+2))* nmats * sizeof(r8)
         Constructor_4d_with_offset_and_diag%start_recv(nb) = (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)&
                                                            +4*(d2+2)+3*(d1+2))* nmats * sizeof(r8)
      end if



      if ((x+1 /= npx+1) .and. (y+1 /= npy+1) .and. (z+1 /= npz+1)) then
         nb = nb + 1
         Constructor_4d_with_offset_and_diag%types_send(nb) = mats_t
         Constructor_4d_with_offset_and_diag%types_recv(nb) = MPI_DOUBLE_PRECISION
         Constructor_4d_with_offset_and_diag%sizes_send(nb) = 1
         Constructor_4d_with_offset_and_diag%sizes_recv(nb) = nmats
         Constructor_4d_with_offset_and_diag%start_send(nb) = ((d3 - dim_offset) * (d2+2) * (d1+2) + (d2 - dim_offset)&
                                                           * (d1+2) + (d1 - dim_offset))* nmats * sizeof(r8)
         Constructor_4d_with_offset_and_diag%start_recv(nb) = (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)&
                                                            +4*(d2+2)+4*(d1+2))* nmats * sizeof(r8)
      end if
      if ((x+1 /= npx+1) .and. (y+1 /= npy+1) .and. (z-1 /= 0)) then
         nb = nb + 1
         Constructor_4d_with_offset_and_diag%types_send(nb) = mats_t
         Constructor_4d_with_offset_and_diag%types_recv(nb) = MPI_DOUBLE_PRECISION
         Constructor_4d_with_offset_and_diag%sizes_send(nb) = 1
         Constructor_4d_with_offset_and_diag%sizes_recv(nb) = nmats
         Constructor_4d_with_offset_and_diag%start_send(nb) = ((1 + dim_offset) * (d2+2) * (d1+2) + (d2 - dim_offset)&
                                                           * (d1+2) + (d1 - dim_offset))* nmats * sizeof(r8)
         Constructor_4d_with_offset_and_diag%start_recv(nb) = (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+&
                                                            4*(d2+2)+4*(d1+2)+1)* nmats * sizeof(r8)
      end if

      if ((x+1 /= npx+1) .and. (y-1 /= 0) .and. (z+1 /= npz+1)) then
         nb = nb + 1
         Constructor_4d_with_offset_and_diag%types_send(nb) = mats_t
         Constructor_4d_with_offset_and_diag%types_recv(nb) = MPI_DOUBLE_PRECISION
         Constructor_4d_with_offset_and_diag%sizes_send(nb) = 1
         Constructor_4d_with_offset_and_diag%sizes_recv(nb) = nmats
         Constructor_4d_with_offset_and_diag%start_send(nb) = ((d3 - dim_offset) * (d2+2) * (d1+2) + (1 + dim_offset)&
                                                           * (d1+2) + (d1 - dim_offset))* nmats * sizeof(r8)
         Constructor_4d_with_offset_and_diag%start_recv(nb) = (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+&
                                                           4*(d3+2)+4*(d2+2)+4*(d1+2)+2)* nmats * sizeof(r8)
      end if

      if ((x+1 /= npx+1) .and. (y-1 /= 0) .and. (z-1 /= 0)) then
         nb = nb + 1
         Constructor_4d_with_offset_and_diag%types_send(nb) = mats_t
         Constructor_4d_with_offset_and_diag%types_recv(nb) = MPI_DOUBLE_PRECISION
         Constructor_4d_with_offset_and_diag%sizes_send(nb) = 1
         Constructor_4d_with_offset_and_diag%sizes_recv(nb) = nmats
         Constructor_4d_with_offset_and_diag%start_send(nb) = ((1 + dim_offset) * (d2+2) * (d1+2) + (1 + dim_offset) &
                                                           * (d1+2) + (d1 - dim_offset))* nmats * sizeof(r8)
         Constructor_4d_with_offset_and_diag%start_recv(nb) = (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)&
                                                            +4*(d2+2)+4*(d1+2)+3)* nmats * sizeof(r8)
      end if

      if ((x-1 /= 0) .and. (y+1 /= npy+1) .and. (z+1 /= npz+1)) then
         nb = nb + 1
         Constructor_4d_with_offset_and_diag%types_send(nb) = mats_t
         Constructor_4d_with_offset_and_diag%types_recv(nb) = MPI_DOUBLE_PRECISION
         Constructor_4d_with_offset_and_diag%sizes_send(nb) = 1
         Constructor_4d_with_offset_and_diag%sizes_recv(nb) = nmats
         Constructor_4d_with_offset_and_diag%start_send(nb) = ((d3 - dim_offset) * (d2+2) * (d1+2) + (d2 - dim_offset)&
                                                           * (d1+2) + (1 + dim_offset))* nmats * sizeof(r8)
         Constructor_4d_with_offset_and_diag%start_recv(nb) = (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+&
                                                            4*(d2+2)+4*(d1+2)+4)* nmats * sizeof(r8)
      end if

      if ((x-1 /= 0) .and. (y+1 /= npy+1) .and. (z-1 /= 0)) then
         nb = nb + 1
         Constructor_4d_with_offset_and_diag%types_send(nb) = mats_t
         Constructor_4d_with_offset_and_diag%types_recv(nb) = MPI_DOUBLE_PRECISION
         Constructor_4d_with_offset_and_diag%sizes_send(nb) = 1
         Constructor_4d_with_offset_and_diag%sizes_recv(nb) = nmats
         Constructor_4d_with_offset_and_diag%start_send(nb) = ((1 + dim_offset) * (d2+2) * (d1+2) + (d2 - dim_offset)&
                                                           * (d1+2) + (1 + dim_offset))* nmats * sizeof(r8)
         Constructor_4d_with_offset_and_diag%start_recv(nb) = (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+&
                                                           4*(d2+2)+4*(d1+2)+5)* nmats * sizeof(r8)
      end if

      if ((x-1 /= 0) .and. (y-1 /= 0) .and. (z+1 /= npz+1)) then
         nb = nb + 1
         Constructor_4d_with_offset_and_diag%types_send(nb) = mats_t
         Constructor_4d_with_offset_and_diag%types_recv(nb) = MPI_DOUBLE_PRECISION
         Constructor_4d_with_offset_and_diag%sizes_send(nb) = 1
         Constructor_4d_with_offset_and_diag%sizes_recv(nb) = nmats
         Constructor_4d_with_offset_and_diag%start_send(nb) = ((d3 - dim_offset) * (d2+2) * (d1+2) + (1 + dim_offset) &
                                                           * (d1+2) + (1 + dim_offset))* nmats * sizeof(r8)
         Constructor_4d_with_offset_and_diag%start_recv(nb) = (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+&
                                                            4*(d3+2)+4*(d2+2)+4*(d1+2)+6)* nmats * sizeof(r8)
      end if

      if ((x-1 /= 0) .and. (y-1 /= 0) .and. (z-1 /= 0)) then
         nb = nb + 1
         Constructor_4d_with_offset_and_diag%types_send(nb) = mats_t
         Constructor_4d_with_offset_and_diag%types_recv(nb) = MPI_DOUBLE_PRECISION
         Constructor_4d_with_offset_and_diag%sizes_send(nb) = 1
         Constructor_4d_with_offset_and_diag%sizes_recv(nb) = nmats
         Constructor_4d_with_offset_and_diag%start_send(nb) = ((1 + dim_offset) * (d2+2) * (d1+2) + (1 + dim_offset)&
                                                           * (d1+2) + (1 + dim_offset))* nmats * sizeof(r8)
         Constructor_4d_with_offset_and_diag%start_recv(nb) = (2*(d2+2)*(d3+2)+2*(d1+2)*(d3+2)+2*(d1+2)*(d2+2)+4*(d3+2)+&
                                                            4*(d2+2)+4*(d1+2)+7)* nmats * sizeof(r8)
      end if



  end function



  type(communication_parameters_t) function Constructor(d1, d2, d3, comm, ndim, npx, npy, npz)
     implicit none
     integer, intent(in)                          :: d1, d2, d3
     integer, intent(in)                          :: comm
     integer, intent(in)                          :: ndim
     integer, intent(in)                          :: npx, npy, npz
     integer, dimension (3) :: global_sizes, recv_global_sizes
     integer :: ierror, ierror2, ierror3
     integer :: rear_front_t, left_right_t, up_down_t, recv_rear_front_t, recv_left_right_t, recv_up_down_t
     integer :: i, j ,k, recv_sum
     integer (kind=MPI_address_kind) :: lb, sizeofreal8
     integer :: d_ndim
     integer, dimension(2)                       :: send_type_left_right, send_type_rear_front, send_type_down_up
     integer, dimension(2)                       :: recv_type_left_right, recv_type_rear_front, recv_type_down_up
     integer, dimension(2)                     :: send_size_left_right, send_size_rear_front, send_size_down_up
     integer, dimension(2)                     :: recv_size_left_right, recv_size_rear_front, recv_size_down_up
     integer (kind=MPI_address_kind), dimension(2) :: send_start_left_right, send_start_rear_front, send_start_down_up
     real(8) :: r8

     Constructor%comm = comm
     Constructor%d1 = d1
     Constructor%d2 = d2
     Constructor%d3 = d3

     global_sizes = (/d1, d2, d3/)
     recv_global_sizes = (/d1+2, d2+2, d3+2/)

     call MPI_type_create_subarray(3, global_sizes, (/1, d2, d3/), (/0, 0, 0/), &
        MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, left_right_t, ierror)
     call MPI_type_create_subarray(3, global_sizes, (/d1, 1, d3/), (/0, 0, 0/), &
        MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, rear_front_t, ierror)
     call MPI_type_create_subarray(3, global_sizes, (/d1, d2, 1/), (/0, 0, 0/), &
        MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, up_down_t, ierror)

     call MPI_type_create_subarray(3, recv_global_sizes, (/1, d2, d3/), (/0, 0, 0/), &
        MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, recv_left_right_t, ierror)
     call MPI_type_create_subarray(3, recv_global_sizes, (/d1, 1, d3/), (/0, 0, 0/), &
        MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, recv_rear_front_t, ierror)
     call MPI_type_create_subarray(3, recv_global_sizes, (/d1, d2, 1/), (/0, 0, 0/), &
        MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, recv_up_down_t, ierror)

     call MPI_type_commit(left_right_t, ierror)
     call MPI_type_commit(rear_front_t, ierror)
     call MPI_type_commit(up_down_t, ierror)

     call MPI_type_commit(recv_left_right_t, ierror)
     call MPI_type_commit(recv_rear_front_t, ierror)
     call MPI_type_commit(recv_up_down_t, ierror)
     d_ndim = 2 * ndim
     allocate(Constructor%types_send(d_ndim))
     allocate(Constructor%types_recv(d_ndim))
     allocate(Constructor%sizes_send(d_ndim))
     allocate(Constructor%sizes_recv(d_ndim))
     allocate(Constructor%start_send(d_ndim))
     allocate(Constructor%start_recv(d_ndim))
     i = 0
     if (npx .gt. 1) then
        i = i + 2
        Constructor%types_send(i-1) = left_right_t       
        Constructor%types_send(i)   = left_right_t       
        Constructor%types_recv(i-1) = recv_left_right_t
        Constructor%types_recv(i)   = recv_left_right_t
        Constructor%sizes_send(i-1) = 1
        Constructor%sizes_send(i)   = 1
        Constructor%sizes_recv(i-1) = 1
        Constructor%sizes_recv(i)   = 1
        Constructor%start_send(i-1) = 0 * sizeof(r8)
        Constructor%start_send(i)   = (d1-1) * sizeof(r8)
        Constructor%start_recv(i-1) = (1 * (d2 + 2) * (d1 + 2) + 1 * (d1 + 2)) * sizeof(r8)
        Constructor%start_recv(i)   = (1 * (d2 + 2) * (d1 + 2) + 1 * (d1 + 2) + (d1-1 + 2)) * sizeof(r8)
     end if
     if (npy .gt. 1) then
     i = i + 2
        Constructor%types_send(i-1) = rear_front_t        
        Constructor%types_send(i)   = rear_front_t        
        Constructor%types_recv(i-1) = recv_rear_front_t
        Constructor%types_recv(i)   = recv_rear_front_t
        Constructor%sizes_send(i-1) = 1
        Constructor%sizes_send(i)   = 1
        Constructor%sizes_recv(i-1) = 1
        Constructor%sizes_recv(i)   = 1
        Constructor%start_send(i-1) = 0 * sizeof(r8)
        Constructor%start_send(i)   = (d2-1) * d1  * sizeof(r8)
        Constructor%start_recv(i-1) = (1 * (d2 + 2) * (d1 + 2) + 1) * sizeof(r8)
        Constructor%start_recv(i)   = (1 * (d2 + 2) * (d1 + 2) + 1 + (d2-1 + 2) * (d1 + 2)) * sizeof(r8)
     end if
     if (npz .gt. 1) then
        i = i + 2
        Constructor%types_send(i-1) = up_down_t           
        Constructor%types_send(i)   = up_down_t           
        Constructor%types_recv(i-1) = recv_up_down_t
        Constructor%types_recv(i)   = recv_up_down_t
        Constructor%sizes_send(i-1) = 1
        Constructor%sizes_send(i)   = 1
        Constructor%sizes_recv(i-1) = 1
        Constructor%sizes_recv(i)   = 1
        Constructor%start_send(i-1) = 0 * sizeof(r8)
        Constructor%start_send(i)   = (d3-1) * d2 * d1 * sizeof(r8)
        Constructor%start_recv(i-1) = (1 * (d1 + 2) + 1) * sizeof(r8)
        Constructor%start_recv(i)   = (1 * (d1 + 2) + 1 + (d3-1 + 2) * (d2 + 2) * (d1 + 2)) * sizeof(r8)
     end if

  end function

end module communication_parameters_module
