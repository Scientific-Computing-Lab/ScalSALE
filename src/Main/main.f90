program main
    use problem_module, only : problem_t
    use datafile_module, only : datafile_t
    use mpi
    use omp_lib
    implicit none

    type(problem_t) , allocatable    :: prob
    type(datafile_t) :: df_obj

    integer :: num_args, ix
    integer :: ierr, rank

    character(:), allocatable :: arg
    integer :: arglen, stat
    call get_command_argument(number=1, length=arglen)  ! Assume for simplicity success
    if (arglen == 0) then
        df_obj = datafile_t("../Datafiles/datafile.json")
        else
            allocate (character(arglen) :: arg)
            call get_command_argument(number=1, value=arg, status=stat)
            df_obj = datafile_t(arg)
    end if

    call MPI_init(ierr)
    call mpi_comm_rank(MPI_COMM_WORLD, rank, ierr)

    allocate(prob)



    prob = problem_t(df_obj)
    call prob%Start_calculation()
    call MPI_FINALIZE(ierr)

end program main
