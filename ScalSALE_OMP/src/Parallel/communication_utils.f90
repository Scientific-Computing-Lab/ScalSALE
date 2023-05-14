module communication_utils_module
   implicit none

contains

   subroutine MPI_User_fn_array_by_min_val(invec, inoutvec, len1, datatype)
      real(8), intent(inout) :: invec(len1), inoutvec(len1)
      integer, intent(inout) :: len1, datatype

      integer :: i

      if (invec(1) < inoutvec(1)) then
         do i=1, len1
            inoutvec(i) = invec(i)
         end do
      end if
   end subroutine MPI_User_fn_array_by_min_val

   subroutine MPI_User_fn_array_by_max_val(invec, inoutvec, len1, datatype)
      real(8),intent(inout) :: invec(len1), inoutvec(len1)
      integer,intent(inout) :: len1, datatype
      integer :: i
      if (invec(1) > inoutvec(1)) then
         do i=1, len1
            inoutvec(i) = invec(i)
         end do
      end if
   end subroutine MPI_User_fn_array_by_max_val


   subroutine Create_diag_topology(npx, npy, npz, coords, comm_dist_graph, num_neighbors)
      use mpi
      implicit none
      integer, intent(out) :: comm_dist_graph, num_neighbors
      integer                              :: my_rank
      integer                              :: npx, npy, npz
      integer                              :: x, y, z   
      integer                              :: nb
      integer :: ierr

      integer, dimension(:), allocatable   :: destinations, weights
      integer, dimension(8)    :: destinations1, weights1
      integer, dimension(1)   :: sources, degrees
      integer, dimension(1)   :: sources1, degrees1
      integer, dimension(3)   :: coords

      integer, dimension(:), allocatable :: neighbors, srcw, dest, destw
      integer :: nnodes, nnedges
      integer :: p,q
      logical :: flag
      p = npx
      q = npy
      call Count_neighbors(npx, npy, npz, coords, num_neighbors)
      allocate(destinations(num_neighbors))
      allocate(weights(num_neighbors))
      call MPI_comm_rank(MPI_COMM_WORLD, my_rank, ierr)
      x = coords(1) - 1
      y = coords(2) - 1
      z = coords(3) - 1

      nb = 0


      if (x+1 /= npx) then
         nb = nb + 1
         destinations(nb) = z * (npx * npy) + npx * y + (x+1)
       weights(nb) = 1
      end if
      if (x-1 /= -1) then
         nb = nb + 1
         destinations(nb) = z * (npx * npy) + npx * y + (x-1)
       weights(nb) = 1
      end if


      if (y+1 /= npy) then
         nb = nb + 1
         destinations(nb) = z * (npx * npy) + npx * (y+1) + x
       weights(nb) = 1
      end if
      if (y-1 /= -1) then
         nb = nb + 1
         destinations(nb) = z * (npx * npy) + npx * (y-1) + x
       weights(nb) = 1
      end if

      if (z+1 /= npz) then
         nb = nb + 1
         destinations(nb) = (z+1) * (npx * npy) + npx * y + x
       weights(nb) = 1
      end if
      if (z-1 /= -1) then
         nb = nb + 1
         destinations(nb) = (z-1) * (npx * npy) + npx * y + x
       weights(nb) = 1
      end if

      if (x+1 /= npx .and. y+1 /= npy) then
         nb = nb + 1
         destinations(nb) = z * (npx * npy) + npx * (y+1) + (x+1)
       weights(nb) = 2
      end if
      if (x+1 /= npx .and. y-1 /= -1) then
         nb = nb + 1
         destinations(nb) = z * (npx * npy) + npx * (y-1) + (x+1)
       weights(nb) = 2
      end if
      if (x-1 /= -1 .and. y+1 /= npy) then
         nb = nb + 1
         destinations(nb) = z * (npx * npy) + npx * (y+1) + (x-1)
       weights(nb) = 2
      end if
      if (x-1 /= -1 .and. y-1 /= -1) then
         nb = nb + 1
         destinations(nb) = z * (npx * npy) + npx * (y-1) + (x-1)
       weights(nb) = 2
      end if
      if (x+1 /= npx .and. z+1 /= npz) then
         nb = nb + 1
         destinations(nb) = (z+1) * (npx * npy) + npx * y + (x+1)
       weights(nb) = 2
      end if
      if (x+1 /= npx .and. z-1 /= -1) then
         nb = nb + 1
         destinations(nb) = (z-1) * (npx * npy) + npx * y + (x+1)
       weights(nb) = 2
      end if

      if (x-1 /= -1 .and. z+1 /= npz) then
         nb = nb + 1
         destinations(nb) = (z+1) * (npx * npy) + npx * y + (x-1)
       weights(nb) = 2
      end if
      if (x-1 /= -1 .and. z-1 /= -1) then
         nb = nb + 1
         destinations(nb) = (z-1) * (npx * npy) + npx * y + (x-1)
       weights(nb) = 2
      end if
      if (y+1 /= npy .and. z+1 /= npz) then
         nb = nb + 1
         destinations(nb) = (z+1) * (npx * npy) + npx * (y+1) + x
       weights(nb) = 2
      end if
      if (y+1 /= npy .and. z-1 /= -1) then
         nb = nb + 1
         destinations(nb) = (z-1) * (npx * npy) + npx * (y+1) + x
       weights(nb) = 2
      end if
      if (y-1 /= -1 .and. z+1 /= npz) then
         nb = nb + 1
         destinations(nb) = (z+1) * (npx * npy) + npx * (y-1) + x
       weights(nb) = 2
      end if
      if (y-1 /= -1 .and. z-1 /= -1) then
         nb = nb + 1
         destinations(nb) = (z-1) * (npx * npy) + npx * (y-1) + x
       weights(nb) = 2
      end if
      if (x+1 /= npx .and. y+1 /= npy .and. z+1 /= npz) then
         nb = nb + 1
         destinations(nb) = (z+1) * (npx * npy) + npx * (y+1) + (x+1)
       weights(nb) = 3
      end if
      if (x+1 /= npx .and. y+1 /= npy .and. z-1 /= -1) then
         nb = nb + 1
         destinations(nb) = (z-1) * (npx * npy) + npx * (y+1) + (x+1)
       weights(nb) = 3
      end if

      if (x+1 /= npx .and. y-1 /= -1 .and. z+1 /= npz) then
         nb = nb + 1
         destinations(nb) = (z+1) * (npx * npy) + npx * (y-1) + (x+1)
       weights(nb) = 3
      end if
      if (x+1 /= npx .and. y-1 /= -1 .and. z-1 /= -1) then
         nb = nb + 1
         destinations(nb) = (z-1) * (npx * npy) + npx * (y-1) + (x+1)
       weights(nb) = 3
      end if

      if (x-1 /= -1 .and. y+1 /= npy .and. z+1 /= npz) then
         nb = nb + 1
         destinations(nb) = (z+1) * (npx * npy) + npx * (y+1) + (x-1)
       weights(nb) = 3
      end if
      if (x-1 /= -1 .and. y+1 /= npy .and. z-1 /= -1) then
         nb = nb + 1
         destinations(nb) = (z-1) * (npx * npy) + npx * (y+1) + (x-1)
       weights(nb) = 3
      end if
      if (x-1 /= -1 .and. y-1 /= -1 .and. z+1 /= npz) then
         nb = nb + 1
         destinations(nb) = (z+1) * (npx * npy) + npx * (y-1) + (x-1)
       weights(nb) = 3
      end if
      if (x-1 /= -1 .and. y-1 /= -1 .and. z-1 /= -1) then
         nb = nb + 1
         destinations(nb) = (z-1) * (npx * npy) + npx * (y-1) + (x-1)
       weights(nb) = 3
      end if

      call MPI_dist_graph_create_adjacent(MPI_COMM_WORLD, nb, destinations, weights, &
                                          nb, destinations, weights, &
                                          MPI_info_null, .true., comm_dist_graph, ierr)
      call MPI_comm_rank(comm_dist_graph, my_rank, ierr)
deallocate(destinations)
      deallocate(weights)



   end subroutine Create_diag_topology

   subroutine Count_neighbors(npx, npy, npz, coords, num_neighbors)
      implicit none
      integer, intent(in)                  :: npx, npy, npz
      integer, dimension(3), intent(in)    :: coords
      integer, intent(out)                 :: num_neighbors
      integer                              :: x, y, z   

      x = coords(1)
      y = coords(2)
      z = coords(3)

      num_neighbors = 0


      if (x+1 /= npx+1) then
         num_neighbors = num_neighbors + 1
      end if
      if (x-1 /= 0) then
         num_neighbors = num_neighbors + 1
      end if

      if (y+1 /= npy+1) then
         num_neighbors = num_neighbors + 1
      end if
      if (y-1 /= 0) then
         num_neighbors = num_neighbors + 1
      end if

      if (z+1 /= npz+1) then
         num_neighbors = num_neighbors + 1
      end if
      if (z-1 /= 0) then
         num_neighbors = num_neighbors + 1
      end if

      if (x+1 /= npx+1 .and. y+1 /= npy+1) then
         num_neighbors = num_neighbors + 1
      end if
      if (x+1 /= npx+1 .and. y-1 /= 0) then
         num_neighbors = num_neighbors + 1
      end if
      if (x-1 /= 0 .and. y+1 /= npy+1) then
         num_neighbors = num_neighbors + 1
      end if
      if (x-1 /= 0 .and. y-1 /= 0) then
         num_neighbors = num_neighbors + 1
      end if
      if (x+1 /= npx+1 .and. z+1 /= npz+1) then
         num_neighbors = num_neighbors + 1
      end if
      if (x+1 /= npx+1 .and. z-1 /= 0) then
         num_neighbors = num_neighbors + 1
      end if
      if (x-1 /= 0 .and. z+1 /= npz+1) then
         num_neighbors = num_neighbors + 1
      end if
      if (x-1 /= 0 .and. z-1 /= 0) then
         num_neighbors = num_neighbors + 1
      end if
      if (y+1 /= npy+1 .and. z+1 /= npz+1) then
         num_neighbors = num_neighbors + 1
      end if
      if (y+1 /= npy+1 .and. z-1 /= 0) then
         num_neighbors = num_neighbors + 1
      end if
      if (y-1 /= 0 .and. z+1 /= npz+1) then
         num_neighbors = num_neighbors + 1
      end if
      if (y-1 /= 0 .and. z-1 /= 0) then
         num_neighbors = num_neighbors + 1
      end if
      if (x+1 /= npx+1 .and. y+1 /= npy+1 .and. z+1 /= npz+1) then
         num_neighbors = num_neighbors + 1
      end if
      if (x+1 /= npx+1 .and. y+1 /= npy+1 .and. z-1 /= 0) then
         num_neighbors = num_neighbors + 1
      end if
      if (x+1 /= npx+1 .and. y-1 /= 0 .and. z+1 /= npz+1) then
         num_neighbors = num_neighbors + 1
      end if
      if (x+1 /= npx+1 .and. y-1 /= 0 .and. z-1 /= 0) then
         num_neighbors = num_neighbors + 1
      end if

      if (x-1 /= 0 .and. y+1 /= npy+1 .and. z+1 /= npz+1) then
         num_neighbors = num_neighbors + 1
      end if
      if (x-1 /= 0 .and. y+1 /= npy+1 .and. z-1 /= 0) then
         num_neighbors = num_neighbors + 1
      end if
      if (x-1 /= 0 .and. y-1 /= 0 .and. z+1 /= npz+1) then
         num_neighbors = num_neighbors + 1
      end if
      if (x-1 /= 0 .and. y-1 /= 0 .and. z-1 /= 0) then
         num_neighbors = num_neighbors + 1
      end if

   end subroutine Count_neighbors


end module communication_utils_module
