module boundary_parameters_module
    use data_module                  , only : data_t
    use parallel_parameters_module   , only : parallel_parameters_t
    use communication_module, only : communication_t
    use communication_parameters_module, only : communication_parameters_t
    implicit none
    private

    type, public :: boundary_parameters_t
        private
        integer, public :: nxp,nyp,nzp,dimension

        integer, dimension(:), pointer,public         :: edge_num       

        real(8)      , dimension(:), pointer         :: angles         

        real(8), dimension(3), public :: normal_i1     
        real(8), dimension(3), public :: normal_inxp   
        real(8), dimension(3), public :: normal_j1     
        real(8), dimension(3), public :: normal_jnyp   
        real(8), dimension(3), public :: normal_k1     
        real(8), dimension(3), public :: normal_knzp   

        real(8), dimension(:), pointer :: edge_vector_x   
        real(8), dimension(:), pointer :: edge_vector_y   
        real(8), dimension(:), pointer :: edge_vector_z   

        real(8), public :: plane_const_i1     
        real(8), public :: plane_const_inxp   
        real(8), public :: plane_const_j1     
        real(8), public :: plane_const_jnyp   
        real(8), public :: plane_const_k1     
        real(8), public :: plane_const_knzp   

        logical, public               :: is_wall_x_top    
        logical, public               :: is_wall_x_bot    
        logical, public               :: is_wall_y_top    
        logical, public               :: is_wall_y_bot    
        logical, public               :: is_wall_z_top    
        logical, public               :: is_wall_z_bot    
        logical, public               :: is_parallel
        integer, public               :: mesh_type

        integer, dimension(:), pointer, public :: boundary_type
        type (communication_t),pointer,public         :: mpi_comm       

        type (parallel_parameters_t), pointer, public ::   parallel_params

    contains
        procedure, public :: Calculate_boundary_normal_3d
        procedure, public :: Calculate_edge_vector_3d
        procedure, public :: Set_communication
        procedure, public :: Point_to_edge_vectors
        procedure, public :: Point_to_edges
        procedure, public :: Build_boundary

        procedure, public :: Write_boundary_parameters
        generic :: write(unformatted) => Write_boundary_parameters

        procedure, public :: Read_boundary_parameters
        generic :: read(unformatted) => Read_boundary_parameters

    end type

    interface boundary_parameters_t
        module procedure Constructor

    end interface boundary_parameters_t
contains

    type(boundary_parameters_t) function Constructor(nxp, nyp, nzp, dimension, parallel_params, mesh_type, boundary_conditions)
        implicit none
        integer                  , intent(in)                   :: nxp           
        integer                  , intent(in)                   :: nyp          
        integer                  , intent(in)                   :: nzp          
        integer                  , intent(in)                   :: dimension          
        type(parallel_parameters_t), pointer, intent(inout)     :: parallel_params
        integer                  , intent(in)                   :: mesh_type
        integer, dimension(:), allocatable, intent(in)          :: boundary_conditions
        Constructor%nxp = nxp
        Constructor%nyp = nyp
        Constructor%nzp = nzp
        Constructor%dimension = dimension
        Constructor%mesh_type = mesh_type
        Constructor%parallel_params => parallel_params
        Constructor%is_parallel = parallel_params%is_parallel
        allocate (Constructor%angles(10))
        Constructor%angles = 90d0

        if (dimension == 3) then
            Constructor%normal_i1 = 0
            Constructor%normal_inxp = 0
            Constructor%normal_j1 = 0
            Constructor%normal_jnyp = 0
            Constructor%normal_k1 = 0
            Constructor%normal_knzp = 0
            Constructor%plane_const_i1 = 0
            Constructor%plane_const_inxp = 0
            Constructor%plane_const_j1 = 0
            Constructor%plane_const_jnyp = 0
            Constructor%plane_const_k1 = 0
            Constructor%plane_const_knzp = 0
            allocate (Constructor%edge_vector_x (12))
            allocate (Constructor%edge_vector_y (12))
            allocate (Constructor%edge_vector_z (12))
        else
        end if
        Constructor%is_wall_x_top = parallel_params%is_wall_x_top
        Constructor%is_wall_x_bot = parallel_params%is_wall_x_bot
        Constructor%is_wall_y_top = parallel_params%is_wall_y_top
        Constructor%is_wall_y_bot = parallel_params%is_wall_y_bot
        Constructor%is_wall_z_top = parallel_params%is_wall_z_top
        Constructor%is_wall_z_bot = parallel_params%is_wall_z_bot
        call Constructor%Build_boundary(boundary_conditions)

    end function

    subroutine Calculate_edge_vector_3d(this, coordinates_values)
        use geometry_module, only: Plane_normal
        implicit none
        class(boundary_parameters_t), intent(inout)        :: this
        type(data_t), dimension(:), pointer,intent(inout)  :: coordinates_values


        real(8), dimension(:, :, :), pointer :: x      
        real(8), dimension(:, :, :), pointer :: y      
        real(8), dimension(:, :, :), pointer :: z      

        real(8), dimension(:),allocatable :: pt_1, pt_2, pt_3  
        real(8) :: dx, dy, dz 
        real(8) :: length 
        real(8) :: sender_id
        integer :: my_x_dim, my_y_dim, my_z_dim, np, npz, npx ,npy,nxp,nyp,nzp

        my_x_dim = this%parallel_params%my_coords(1)
        my_y_dim = this%parallel_params%my_coords(2)
        my_z_dim = this%parallel_params%my_coords(3)

        np = this%parallel_params%np
        npx = this%parallel_params%npx
        npy = this%parallel_params%npy
        npz = this%parallel_params%npz
        call coordinates_values(1)%Point_to_data (x)
        call coordinates_values(2)%Point_to_data (y)
        call coordinates_values(3)%Point_to_data (z)

        nxp = this%nxp
        nyp = this%nyp
        nzp = this%nzp
        allocate(pt_1(3))
        allocate(pt_2(3))
        allocate(pt_3(3))

        sender_id = -1
        if ( (my_x_dim == npx .and. my_y_dim == 1 .and. my_z_dim == 1) .or. np == 1) then
            pt_1(1) = x(nxp, 1, 1)
            pt_1(2) = y(nxp, 1, 1)
            pt_1(3) = z(nxp, 1, 1)
            sender_id = this%parallel_params%my_rank
        end if
        call this%mpi_comm%Send_get_array_by_max_val(sender_id, pt_1)
        sender_id = -1
        if ( (my_x_dim == 1 .and. my_y_dim == 1 .and. my_z_dim == 1) .or. np == 1) then

            pt_2(1) = x(1, 1, 1)
            pt_2(2) = y(1, 1, 1)
            pt_2(3) = z(1, 1, 1)
            sender_id = this%parallel_params%my_rank
        end if
        call this%mpi_comm%Send_get_array_by_max_val(sender_id, pt_2)

        dx = pt_1(1) - pt_2(1)
        dy = pt_1(2) - pt_2(2)
        dz = pt_1(3) - pt_2(3)
        length = sqrt(dx**2 + dy**2 + dz**2)
        this%edge_vector_x(1) = dx/length
        this%edge_vector_y(1) = dy/length
        this%edge_vector_z(1) = dz/length
        sender_id = -1
        if ( (my_x_dim == npx .and. my_y_dim == npy .and. my_z_dim == 1) .or. np == 1) then
            pt_1(1) = x(nxp, nyp, 1)
            pt_1(2) = y(nxp, nyp, 1)
            pt_1(3) = z(nxp, nyp, 1)
            sender_id = this%parallel_params%my_rank
        end if
        call this%mpi_comm%Send_get_array_by_max_val(sender_id, pt_1)
        sender_id = -1
        if ( (my_x_dim == npx .and. my_y_dim == 1 .and. my_z_dim == 1) .or. np == 1) then
            pt_2(1) = x(nxp, 1, 1)
            pt_2(2) = y(nxp, 1, 1)
            pt_2(3) = z(nxp, 1, 1)
            sender_id = this%parallel_params%my_rank
        end if
        call this%mpi_comm%Send_get_array_by_max_val(sender_id, pt_2)

        dx = pt_1(1) - pt_2(1)
        dy = pt_1(2) - pt_2(2)
        dz = pt_1(3) - pt_2(3)
        length = sqrt(dx**2 + dy**2 + dz**2)
        this%edge_vector_x(2) = dx/length
        this%edge_vector_y(2) = dy/length
        this%edge_vector_z(2) = dz/length
        sender_id = -1
        if ( (my_x_dim == npx .and. my_y_dim == npy .and. my_z_dim == 1) .or. np == 1) then
            pt_1(1)=X(nxp,nyp,1)
            pt_1(2)=Y(nxp,nyp,1)
            pt_1(3)=Z(nxp,nyp,1)
            sender_id = this%parallel_params%my_rank
        end if
        call this%mpi_comm%Send_get_array_by_max_val(sender_id, pt_1)
        sender_id = -1
        if ( (my_x_dim == npx .and. my_y_dim == npy .and. my_z_dim == 1) .or. np == 1) then
            pt_2(1) = x(1, nyp, 1)
            pt_2(2) = y(1, nyp, 1)
            pt_2(3) = z(1, nyp, 1)
            sender_id = this%parallel_params%my_rank
        end if
        call this%mpi_comm%Send_get_array_by_max_val(sender_id, pt_2)

        dx = pt_1(1) - pt_2(1)
        dy = pt_1(2) - pt_2(2)
        dz = pt_1(3) - pt_2(3)
        length = sqrt(dx**2 + dy**2 + dz**2)
        this%edge_vector_x(3) = dx/length
        this%edge_vector_y(3) = dy/length
        this%edge_vector_z(3) = dz/length
        sender_id = -1
        if ( (my_x_dim == 1 .and. my_y_dim == npy .and. my_z_dim == 1) .or. np == 1) then
            pt_1(1) = x(1, nyp, 1)
            pt_1(2) = y(1, nyp, 1)
            pt_1(3) = z(1, nyp, 1)
            sender_id = this%parallel_params%my_rank
        end if
        call this%mpi_comm%Send_get_array_by_max_val(sender_id, pt_1)
        sender_id = -1
        if ( (my_x_dim == 1 .and. my_y_dim == 1 .and. my_z_dim == 1) .or. np == 1) then
            pt_2(1) = x(1, 1, 1)
            pt_2(2) = y(1, 1, 1)
            pt_2(3) = z(1, 1, 1)
            sender_id = this%parallel_params%my_rank
        end if
        call this%mpi_comm%Send_get_array_by_max_val(sender_id, pt_2)

        dx = pt_1(1) - pt_2(1)
        dy = pt_1(2) - pt_2(2)
        dz = pt_1(3) - pt_2(3)
        length = sqrt(dx**2 + dy**2 + dz**2)
        this%edge_vector_x(4) = dx/length
        this%edge_vector_y(4) = dy/length
        this%edge_vector_z(4) = dz/length
        sender_id = -1
        if ( (my_x_dim == 1 .and. my_y_dim == 1 .and. my_z_dim == npz) .or. np == 1) then
            pt_1(1) = x(1, 1, nzp)
            pt_1(2) = y(1, 1, nzp)
            pt_1(3) = z(1, 1, nzp)
            sender_id = this%parallel_params%my_rank
        end if
        call this%mpi_comm%Send_get_array_by_max_val(sender_id, pt_1)
        sender_id = -1
        if ( (my_x_dim == 1 .and. my_y_dim == 1 .and. my_z_dim == 1) .or. np == 1) then
            pt_2(1) = x(1, 1, 1)
            pt_2(2) = y(1, 1, 1)
            pt_2(3) = z(1, 1, 1)
            sender_id = this%parallel_params%my_rank
        end if
        call this%mpi_comm%Send_get_array_by_max_val(sender_id, pt_2)

        dx = pt_1(1) - pt_2(1)
        dy = pt_1(2) - pt_2(2)
        dz = pt_1(3) - pt_2(3)
        length = sqrt(dx**2 + dy**2 + dz**2)
        this%edge_vector_x(5) = dx/length
        this%edge_vector_y(5) = dy/length
        this%edge_vector_z(5) = dz/length
        sender_id = -1
        if ( (my_x_dim == npx .and. my_y_dim == 1 .and. my_z_dim == npz) .or. np == 1) then
            pt_1(1) = x(nxp, 1, nzp)
            pt_1(2) = y(nxp, 1, nzp)
            pt_1(3) = z(nxp, 1, nzp)
            sender_id = this%parallel_params%my_rank
        end if
        call this%mpi_comm%Send_get_array_by_max_val(sender_id, pt_1)
        sender_id = -1
        if ( (my_x_dim == npx .and. my_y_dim == 1 .and. my_z_dim == 1) .or. np == 1) then
            pt_2(1) = x(nxp, 1, 1)
            pt_2(2) = y(nxp, 1, 1)
            pt_2(3) = z(nxp, 1, 1)
            sender_id = this%parallel_params%my_rank
        end if
        call this%mpi_comm%Send_get_array_by_max_val(sender_id, pt_2)

        dx = pt_1(1) - pt_2(1)
        dy = pt_1(2) - pt_2(2)
        dz = pt_1(3) - pt_2(3)
        length = sqrt(dx**2 + dy**2 + dz**2)
        this%edge_vector_x(6) = dx/length
        this%edge_vector_y(6) = dy/length
        this%edge_vector_z(6) = dz/length
        sender_id = -1
        if ( (my_x_dim == npx .and. my_y_dim == npy .and. my_z_dim == npz) .or. np == 1) then
            pt_1(1) = x(nxp, nyp, nzp)
            pt_1(2) = y(nxp, nyp, nzp)
            pt_1(3) = z(nxp, nyp, nzp)
            sender_id = this%parallel_params%my_rank
        end if
        call this%mpi_comm%Send_get_array_by_max_val(sender_id, pt_1)
        sender_id = -1
        if ( (my_x_dim == npx .and. my_y_dim == npy .and. my_z_dim == 1) .or. np == 1) then
            pt_2(1) = x(nxp,nyp,1)
            pt_2(2) = y(nxp,nyp,1)
            pt_2(3) = z(nxp,nyp,1)
            sender_id = this%parallel_params%my_rank
        end if
        call this%mpi_comm%Send_get_array_by_max_val(sender_id, pt_2)

        dx = pt_1(1) - pt_2(1)
        dy = pt_1(2) - pt_2(2)
        dz = pt_1(3) - pt_2(3)
        length = sqrt(dx**2 + dy**2 + dz**2)
        this%edge_vector_x(7) = dx/length
        this%edge_vector_y(7) = dy/length
        this%edge_vector_z(7) = dz/length
        sender_id = -1
        if ( (my_x_dim == 1 .and. my_y_dim == npy .and. my_z_dim == npz) .or. np == 1) then
            pt_1(1) = x(1, nyp, nzp)
            pt_1(2) = y(1, nyp, nzp)
            pt_1(3) = z(1, nyp, nzp)
            sender_id = this%parallel_params%my_rank
        end if
        call this%mpi_comm%Send_get_array_by_max_val(sender_id, pt_1)
        sender_id = -1
        if ( (my_x_dim == 1 .and. my_y_dim == npy .and. my_z_dim == 1) .or. np == 1) then
            pt_2(1) = x(1, nyp, 1)
            pt_2(2) = y(1, nyp, 1)
            pt_2(3) = z(1, nyp, 1)
            sender_id = this%parallel_params%my_rank
        end if
        call this%mpi_comm%Send_get_array_by_max_val(sender_id, pt_2)

        dx = pt_1(1) - pt_2(1)
        dy = pt_1(2) - pt_2(2)
        dz = pt_1(3) - pt_2(3)
        length = sqrt(dx**2 + dy**2 + dz**2)
        this%edge_vector_x(8) = dx/length
        this%edge_vector_y(8) = dy/length
        this%edge_vector_z(8) = dz/length
        sender_id = -1
        if ( (my_x_dim == npx .and. my_y_dim == 1 .and. my_z_dim == npz) .or. np == 1) then
            pt_1(1) = x(nxp, 1, nzp)
            pt_1(2) = y(nxp, 1, nzp)
            pt_1(3) = z(nxp, 1, nzp)
            sender_id = this%parallel_params%my_rank
        end if
        call this%mpi_comm%Send_get_array_by_max_val(sender_id, pt_1)
        sender_id = -1
        if ( (my_x_dim == 1 .and. my_y_dim == 1 .and. my_z_dim == npz) .or. np == 1) then
            pt_2(1) = x(1, 1, nzp)
            pt_2(2) = y(1, 1, nzp)
            pt_2(3) = z(1, 1, nzp)
            sender_id = this%parallel_params%my_rank
        end if
        call this%mpi_comm%Send_get_array_by_max_val(sender_id, pt_2)

        dx = pt_1(1) - pt_2(1)
        dy = pt_1(2) - pt_2(2)
        dz = pt_1(3) - pt_2(3)
        length = sqrt(dx**2 + dy**2 + dz**2)
        this%edge_vector_x(9) = dx/length
        this%edge_vector_y(9) = dy/length
        this%edge_vector_z(9) = dz/length
        sender_id = -1
        if ( (my_x_dim == npx .and. my_y_dim == npy .and. my_z_dim == npz) .or. np == 1) then
            pt_1(1) = x(nxp, nyp, nzp)
            pt_1(2) = y(nxp, nyp, nzp)
            pt_1(3) = z(nxp, nyp, nzp)
            sender_id = this%parallel_params%my_rank
        end if
        call this%mpi_comm%Send_get_array_by_max_val(sender_id, pt_1)
        sender_id = -1
        if ( (my_x_dim == npx .and. my_y_dim == 1 .and. my_z_dim == npz) .or. np == 1) then
            pt_2(1) = x(nxp, 1, nzp)
            pt_2(2) = y(nxp, 1, nzp)
            pt_2(3) = z(nxp, 1, nzp)
            sender_id = this%parallel_params%my_rank
        end if
        call this%mpi_comm%Send_get_array_by_max_val(sender_id, pt_2)

        dx = pt_1(1) - pt_2(1)
        dy = pt_1(2) - pt_2(2)
        dz = pt_1(3) - pt_2(3)
        length = sqrt(dx**2 + dy**2 + dz**2)
        this%edge_vector_x(10) = dx/length
        this%edge_vector_y(10) = dy/length
        this%edge_vector_z(10) = dz/length
        sender_id = -1
        if ( (my_x_dim == npx .and. my_y_dim == npy .and. my_z_dim == npz) .or. np == 1) then
            pt_1(1) = x(nxp, nyp, nzp)
            pt_1(2) = y(nxp, nyp, nzp)
            pt_1(3) = z(nxp, nyp, nzp)
            sender_id = this%parallel_params%my_rank
        end if
        call this%mpi_comm%Send_get_array_by_max_val(sender_id, pt_1)
        sender_id = -1
        if ( (my_x_dim == 1 .and. my_y_dim == npy .and. my_z_dim == npz) .or. np == 1) then
            pt_2(1) = x(1, nyp, nzp)
            pt_2(2) = y(1, nyp, nzp)
            pt_2(3) = z(1, nyp, nzp)
            sender_id = this%parallel_params%my_rank
        end if
        call this%mpi_comm%Send_get_array_by_max_val(sender_id, pt_2)

        dx = pt_1(1) - pt_2(1)
        dy = pt_1(2) - pt_2(2)
        dz = pt_1(3) - pt_2(3)
        length = sqrt(dx**2 + dy**2 + dz**2)
        this%edge_vector_x(11) = dx/length
        this%edge_vector_y(11) = dy/length
        this%edge_vector_z(11) = dz/length
        sender_id = -1
        if ( (my_x_dim == 1 .and. my_y_dim == npy .and. my_z_dim == npz) .or. np == 1) then
            pt_1(1) = x(1, nyp, nzp)
            pt_1(2) = y(1, nyp, nzp)
            pt_1(3) = z(1, nyp, nzp)
            sender_id = this%parallel_params%my_rank
        end if
        call this%mpi_comm%Send_get_array_by_max_val(sender_id, pt_1)
        sender_id = -1
        if ( (my_x_dim == 1 .and. my_y_dim == 1 .and. my_z_dim == npz) .or. np == 1) then
            pt_2(1) = x(1, 1, nzp)
            pt_2(2) = y(1, 1, nzp)
            pt_2(3) = z(1, 1, nzp)
            sender_id = this%parallel_params%my_rank
        end if
        call this%mpi_comm%Send_get_array_by_max_val(sender_id, pt_2)

        dx = pt_1(1) - pt_2(1)
        dy = pt_1(2) - pt_2(2)
        dz = pt_1(3) - pt_2(3)
        length = sqrt(dx**2 + dy**2 + dz**2)
        this%edge_vector_x(12) = dx/length
        this%edge_vector_y(12) = dy/length
        this%edge_vector_z(12) = dz/length
        deallocate(pt_1)
        deallocate(pt_2)
        deallocate(pt_3)
    end subroutine Calculate_edge_vector_3d


    subroutine Calculate_boundary_normal_3d(this, coordinates_values)
        use geometry_module, only: Plane_normal
        implicit none
        class(boundary_parameters_t), intent(inout)        :: this
        type(data_t), dimension(:), pointer,intent(inout)  :: coordinates_values


        real(8), dimension(:,:,:), pointer :: x,y,z
        real(8) :: sender_id
        integer :: my_x_dim, my_y_dim, my_z_dim, np, npz, npx ,npy

        real(8), dimension(:),allocatable :: pt_1, pt_2, pt_3  
        call coordinates_values(1)%Point_to_data (x)
        call coordinates_values(2)%Point_to_data (y)
        call coordinates_values(3)%Point_to_data (z)

        my_x_dim = this%parallel_params%my_coords(1)
        my_y_dim = this%parallel_params%my_coords(2)
        my_z_dim = this%parallel_params%my_coords(3)
        np = this%parallel_params%np
        npx = this%parallel_params%npx
        npy = this%parallel_params%npy
        npz = this%parallel_params%npz
        allocate(pt_1(3))
        allocate(pt_2(3))
        allocate(pt_3(3))
        sender_id = -1
        if ( (my_x_dim == 1 .and. my_y_dim == 1 .and. my_z_dim == 1) .or. np == 1) then
            pt_1(1) = x(1, 1, 1)
            pt_1(2) = y(1, 1, 1)
            pt_1(3) = z(1, 1, 1)
            sender_id = this%parallel_params%my_rank
        end if


        call this%mpi_comm%Send_get_array_by_max_val(sender_id, pt_1)
        sender_id = -1
        if ( (my_x_dim == 1 .and. my_y_dim == 1 .and. my_z_dim == npz) .or. np == 1) then
            pt_2(1) = x(1, 1, this%nzp)
            pt_2(2) = y(1, 1, this%nzp)
            pt_2(3) = z(1, 1, this%nzp)
            sender_id = this%parallel_params%my_rank
        end if
        call this%mpi_comm%Send_get_array_by_max_val(sender_id, pt_2)
        sender_id = - 1
        if ( (my_x_dim == 1 .and. my_y_dim == npy .and. my_z_dim == 1) .or. np == 1) then

            pt_3(1) = x(1, this%nyp, 1)
            pt_3(2) = y(1, this%nyp, 1)
            pt_3(3) = z(1, this%nyp, 1)
            sender_id = this%parallel_params%my_rank
        end if
        call this%mpi_comm%Send_get_array_by_max_val(sender_id, pt_3)
        call Plane_normal(pt_1(1), pt_1(2), pt_1(3), pt_2(1), pt_2(2), pt_2(3), pt_3(1), pt_3(2), pt_3(3), &
            this%normal_i1(1), this%normal_i1(2), this%normal_i1(3))
        this%plane_const_i1 = this%normal_i1(1) * pt_1(1) + this%normal_i1(2) * pt_1(2) + this%normal_i1(3) * pt_1(3)

        sender_id = - 1
        if ( (my_x_dim == npx .and. my_y_dim == 1 .and. my_z_dim == 1) .or. np == 1) then
            pt_1(1) = x(this%nxp, 1, 1)
            pt_1(2) = y(this%nxp, 1, 1)
            pt_1(3) = z(this%nxp, 1, 1)
            sender_id = this%parallel_params%my_rank
        end if
        call this%mpi_comm%Send_get_array_by_max_val(sender_id, pt_1)

        sender_id = - 1
        if ( (my_x_dim == npx .and. my_y_dim == 1 .and. my_z_dim == npz) .or. np == 1) then
            pt_2(1) = x(this%nxp, 1, this%nzp)
            pt_2(2) = y(this%nxp, 1, this%nzp)
            pt_2(3) = z(this%nxp, 1, this%nzp)
            sender_id = this%parallel_params%my_rank
        end if
        call this%mpi_comm%Send_get_array_by_max_val(sender_id, pt_2)

        sender_id = - 1
        if ( (my_x_dim == npx .and. my_y_dim == npy .and. my_z_dim == 1) .or. np == 1) then
            pt_3(1) = x(this%nxp, this%nyp, 1)
            pt_3(2) = y(this%nxp, this%nyp, 1)
            pt_3(3) = z(this%nxp, this%nyp, 1)
            sender_id = this%parallel_params%my_rank
        end if
        call this%mpi_comm%Send_get_array_by_max_val(sender_id, pt_3)
        call Plane_normal(pt_1(1), pt_1(2), pt_1(3), pt_2(1), pt_2(2), pt_2(3), pt_3(1), pt_3(2), pt_3(3), &
            this%normal_inxp(1), this%normal_inxp(2), this%normal_inxp(3))
        this%plane_const_inxp = this%normal_inxp(1) * pt_1(1) + this%normal_inxp(2) * pt_1(2) + this%normal_inxp(3) * pt_1(3)
        sender_id = - 1
        if ( (my_x_dim == 1 .and. my_y_dim == 1 .and. my_z_dim == 1) .or. np == 1) then
            pt_1(1) = x(1, 1, 1)
            pt_1(2) = y(1, 1, 1)
            pt_1(3) = z(1, 1, 1)
            sender_id = this%parallel_params%my_rank
        end if
        call this%mpi_comm%Send_get_array_by_max_val(sender_id, pt_1)
        sender_id = - 1
        if ( (my_x_dim == 1 .and. my_y_dim == 1 .and. my_z_dim == npz) .or. np == 1) then
            pt_2(1) = x(1, 1, this%nzp)
            pt_2(2) = y(1, 1, this%nzp)
            pt_2(3) = z(1, 1, this%nzp)
            sender_id = this%parallel_params%my_rank
        end if
        call this%mpi_comm%Send_get_array_by_max_val(sender_id, pt_2)

        sender_id = - 1
        if ( (my_x_dim == npx .and. my_y_dim == 1 .and. my_z_dim == 1) .or. np == 1) then
            pt_3(1) = x(this%nxp, 1, 1)
            pt_3(2) = y(this%nxp, 1, 1)
            pt_3(3) = z(this%nxp, 1, 1)
            sender_id = this%parallel_params%my_rank
        end if
        call this%mpi_comm%Send_get_array_by_max_val(sender_id, pt_3)
        call Plane_normal(pt_1(1), pt_1(2), pt_1(3), pt_2(1), pt_2(2), pt_2(3), pt_3(1), pt_3(2), pt_3(3), &
            this%normal_j1(1), this%normal_j1(2), this%normal_j1(3))
        this%plane_const_j1 = this%normal_j1(1) * pt_1(1) + this%normal_j1(2) * pt_1(2) + this%normal_j1(3) * pt_1(3)
        sender_id = -1
        if ( (my_x_dim == 1 .and. my_y_dim == npy .and. my_z_dim == 1) .or. np == 1) then
            pt_1(1) = x(1, this%nyp, 1)
            pt_1(2) = y(1, this%nyp, 1)
            pt_1(3) = z(1, this%nyp, 1)
            sender_id = this%parallel_params%my_rank
        end if
        call this%mpi_comm%Send_get_array_by_max_val(sender_id, pt_1)
        sender_id = -1
        if ( (my_x_dim == 1 .and. my_y_dim == npy .and. my_z_dim == npz) .or. np == 1) then
            pt_2(1) = x(1, this%nyp, this%nzp)
            pt_2(2) = y(1, this%nyp, this%nzp)
            pt_2(3) = z(1, this%nyp, this%nzp)
            sender_id = this%parallel_params%my_rank
        end if
        call this%mpi_comm%Send_get_array_by_max_val(sender_id, pt_2)

        sender_id = -1
        if ( (my_x_dim == npx .and. my_y_dim == npy .and. my_z_dim == 1) .or. np == 1) then
            pt_3(1) = x(this%nxp, this%nyp, 1)
            pt_3(2) = y(this%nxp, this%nyp, 1)
            pt_3(3) = z(this%nxp, this%nyp, 1)
            sender_id = this%parallel_params%my_rank
        end if
        call this%mpi_comm%Send_get_array_by_max_val(sender_id, pt_3)
        call Plane_normal(pt_1(1), pt_1(2), pt_1(3), pt_2(1), pt_2(2), pt_2(3), pt_3(1), pt_3(2), pt_3(3), &
            this%normal_jnyp(1), this%normal_jnyp(2), this%normal_jnyp(3))
        this%plane_const_jnyp = this%normal_jnyp(1) * pt_1(1) + this%normal_jnyp(2) * pt_1(2) + this%normal_jnyp(3) * pt_1(3)

        sender_id = -1
        if ( (my_x_dim == 1 .and. my_y_dim == 1 .and. my_z_dim == 1) .or. np == 1) then
            pt_1(1) = x(1, 1, 1)
            pt_1(2) = y(1, 1, 1)
            pt_1(3) = z(1, 1, 1)
            sender_id = this%parallel_params%my_rank
        end if
        call this%mpi_comm%Send_get_array_by_max_val(sender_id, pt_1)
        sender_id = -1
        if ( (my_x_dim == 1 .and. my_y_dim == npy .and. my_z_dim == 1) .or. np == 1) then
            pt_2(1) = x(1, this%nyp, 1)
            pt_2(2) = y(1, this%nyp, 1)
            pt_2(3) = z(1, this%nyp, 1)
            sender_id = this%parallel_params%my_rank
        end if
        call this%mpi_comm%Send_get_array_by_max_val(sender_id, pt_2)
        sender_id = -1
        if ( (my_x_dim == npx .and. my_y_dim == 1 .and. my_z_dim == 1) .or. np == 1) then
            pt_3(1) = x(this%nxp, 1, 1)
            pt_3(2) = y(this%nxp, 1, 1)
            pt_3(3) = z(this%nxp, 1, 1)
            sender_id = this%parallel_params%my_rank
        end if
        call this%mpi_comm%Send_get_array_by_max_val(sender_id, pt_3)
        call Plane_normal(pt_1(1), pt_1(2), pt_1(3), pt_2(1), pt_2(2), pt_2(3), pt_3(1), pt_3(2), pt_3(3), &
            this%normal_k1(1), this%normal_k1(2), this%normal_k1(3))
        this%plane_const_k1 = this%normal_k1(1) * pt_1(1) + this%normal_k1(2) * pt_1(2) + this%normal_k1(3) * pt_1(3)
        sender_id = - 1
        if ( (my_x_dim == 1 .and. my_y_dim == 1 .and. my_z_dim == npz) .or. np == 1) then
            pt_1(1) = x(1, 1, this%nzp)
            pt_1(2) = y(1, 1, this%nzp)
            pt_1(3) = z(1, 1, this%nzp)
            sender_id = this%parallel_params%my_rank
        end if
        call this%mpi_comm%Send_get_array_by_max_val(sender_id, pt_1)
        sender_id = - 1
        if ( (my_x_dim == 1 .and. my_y_dim == npy .and. my_z_dim == npz) .or. np == 1) then
            pt_2(1) = x(1, this%nyp, this%nzp)
            pt_2(2) = y(1, this%nyp, this%nzp)
            pt_2(3) = z(1, this%nyp, this%nzp)
            sender_id = this%parallel_params%my_rank
        end if
        call this%mpi_comm%Send_get_array_by_max_val(sender_id, pt_2)
        sender_id = - 1
        if ( (my_x_dim == npx .and. my_y_dim == npy .and. my_z_dim == npz) .or. np == 1) then
            pt_3(1) = x(this%nxp, 1, this%nzp)
            pt_3(2) = y(this%nxp, 1, this%nzp)
            pt_3(3) = z(this%nxp, 1, this%nzp)
            sender_id = this%parallel_params%my_rank
        end if
        call this%mpi_comm%Send_get_array_by_max_val(sender_id, pt_3)
        call Plane_normal(pt_1(1), pt_1(2), pt_1(3), pt_2(1), pt_2(2), pt_2(3), pt_3(1), pt_3(2), pt_3(3), &
            this%normal_knzp(1), this%normal_knzp(2), this%normal_knzp(3))
        this%plane_const_knzp = this%normal_knzp(1) * pt_1(1) + this%normal_knzp(2) * pt_1(2) + this%normal_knzp(3) * pt_1(3)

        deallocate(pt_1)
        deallocate(pt_2)
        deallocate(pt_3)
    end subroutine Calculate_boundary_normal_3d


    subroutine Build_boundary(this, bc)
        implicit none
        class(boundary_parameters_t), intent(inout)     :: this
        integer, dimension(:), allocatable, intent(in)  :: bc

        integer :: number_real_boundary = 0
        integer :: counter = 1

        if (this%is_wall_x_bot) number_real_boundary = number_real_boundary + 1
        if (this%is_wall_x_top) number_real_boundary = number_real_boundary + 1
        if (this%is_wall_y_bot) number_real_boundary = number_real_boundary + 1
        if (this%is_wall_y_top) number_real_boundary = number_real_boundary + 1
        if (this%is_wall_z_bot) number_real_boundary = number_real_boundary + 1
        if (this%is_wall_z_top) number_real_boundary = number_real_boundary + 1
        allocate(this%edge_num(number_real_boundary))
        if (this%dimension == 2) allocate(this%boundary_type(4))
        if (this%dimension == 3) allocate(this%boundary_type(6))



        if (this%is_wall_z_bot) then
            this%edge_num(counter) = 5
            counter = counter + 1
        end if

        if (this%is_wall_z_top) then
            this%edge_num(counter) = 6
            counter = counter + 1
        end if

        if (this%is_wall_y_bot) then
            this%edge_num(counter) = 3
            counter = counter + 1
        end if
        if (this%is_wall_y_top) then
            this%edge_num(counter) = 4
            counter = counter + 1
        end if

        if (this%is_wall_x_bot) then
            this%edge_num(counter) = 1
            counter = counter + 1
        end if

        if (this%is_wall_x_top) then
            this%edge_num(counter) = 2
            counter = counter + 1
        end if

        this%boundary_type = -1
        if (this%is_wall_x_bot) then
            if (this%dimension == 2) then
                if (this%mesh_type == 2)  this%boundary_type(1) = bc(1)
            else
                if (this%mesh_type == 2)  this%boundary_type(1) = bc(1)
                if (this%mesh_type == 1)  this%boundary_type(1) = bc(1)
            end if
        end if

        if (this%is_wall_x_top) then
            if (this%dimension == 2) then
                if (this%mesh_type == 2)  this%boundary_type(2) = bc(2)
            else
                if (this%mesh_type == 2)  this%boundary_type(2) = bc(2)
                if (this%mesh_type == 1)  this%boundary_type(2) = bc(2)
            end if
        end if

        if (this%is_wall_y_bot) then
            if (this%dimension == 2) then
                if (this%mesh_type == 2)  this%boundary_type(3) = bc(3)
            else
                if (this%mesh_type == 2)  this%boundary_type(3) = bc(3)
                if (this%mesh_type == 1)  this%boundary_type(3) = bc(3)
            end if
        end if

        if (this%is_wall_y_top) then
            if (this%dimension == 2) then
                if (this%mesh_type == 2)  this%boundary_type(4) = bc(4)
            else
                if (this%mesh_type == 2)  this%boundary_type(4) = bc(4)
                if (this%mesh_type == 1)  this%boundary_type(4) = bc(4)
            end if
        end if

        if (this%is_wall_z_bot) then
            if (this%dimension == 2) then
                if (this%mesh_type == 2)  this%boundary_type(5) = bc(5)
            else
                if (this%mesh_type == 2)  this%boundary_type(5) = bc(5)
                if (this%mesh_type == 1)  this%boundary_type(5) = bc(5)
            end if
        end if

        if (this%is_wall_z_top) then
            if (this%dimension == 2) then
                if (this%mesh_type == 2)  this%boundary_type(6) = bc(6)
            else
                if (this%mesh_type == 2)  this%boundary_type(6) = bc(6)
                if (this%mesh_type == 1)  this%boundary_type(6) = bc(6)
            end if
        end if
    end subroutine Build_boundary





    subroutine Set_communication(this, comm, comm_params)
        class (boundary_parameters_t)            :: this 
        type(communication_t), pointer            :: comm
        type(communication_parameters_t), pointer :: comm_params

        this%mpi_comm => comm
    end subroutine Set_communication


    subroutine Point_to_edge_vectors(this, ptr_x, ptr_y, ptr_z)
        implicit none
        class (boundary_parameters_t)                , intent(in)    :: this   
        real(8), dimension(:), pointer, intent(out)   :: ptr_x  
        real(8), dimension(:), pointer, intent(out)   :: ptr_y  
        real(8), dimension(:), pointer, intent(out)   :: ptr_z  

        ptr_x => this%edge_vector_x
        ptr_y => this%edge_vector_y
        ptr_z => this%edge_vector_z
    end subroutine Point_to_edge_vectors


    subroutine Point_to_edges(this, ptr)
        implicit none
        class (boundary_parameters_t)               , intent(in)      :: this  
        integer, dimension(:), pointer    , intent(out)     :: ptr   

        ptr => this%edge_num
    end subroutine Point_to_edges



    subroutine Point_to_angle(this, ptr)
        implicit none
        class (boundary_parameters_t)               , intent(in)    :: this   
        real(8), dimension(:), pointer    , intent(out)   :: ptr    

        ptr => this% angles
    end subroutine Point_to_angle


    subroutine Write_boundary_parameters(this, unit, iostat, iomsg)
        class (boundary_parameters_t), intent(in) :: this  
        integer,      intent(in)     :: unit
        integer,      intent(out)    :: iostat
        character(*), intent(in out)  :: iomsg

#ifdef DEBUG
        write(*,*) '@@@ in Write_boundary_parameters @@@'
#endif

        write(unit, iostat=iostat, iomsg=iomsg) &
            size(this%edge_num), &
            size(this%angles), &
            size(this%edge_vector_x), &
            size(this%edge_vector_y), &
            size(this%edge_vector_z)

        write(unit, iostat=iostat, iomsg=iomsg) &
            this%nxp, &
            this%nyp, &
            this%nzp, &
            this%dimension, &
            this%edge_num, &
            this%angles, &
            this%normal_i1, &
            this%normal_inxp, &
            this%normal_j1, &
            this%normal_jnyp, &
            this%normal_k1, &
            this%normal_knzp, &
            this%edge_vector_x, &
            this%edge_vector_y, &
            this%edge_vector_z, &
            this%plane_const_i1, &
            this%plane_const_inxp, &
            this%plane_const_j1, &
            this%plane_const_jnyp, &
            this%plane_const_k1, &
            this%plane_const_knzp, &
            this%is_wall_x_top, &
            this%is_wall_x_bot, &
            this%is_wall_y_top, &
            this%is_wall_y_bot, &
            this%is_wall_z_top, &
            this%is_wall_z_bot, &
            this%is_parallel

#ifdef DEBUG
        write(*,*) &
            size(this%edge_num), &
            size(this%angles), &
            size(this%edge_vector_x), &
            size(this%edge_vector_y), &
            size(this%edge_vector_z)

        write(*,*) &
            this%nxp, &
            this%nyp, &
            this%nzp, &
            this%dimension, &
            this%edge_num, &
            this%angles, &
            this%normal_i1, &
            this%normal_inxp, &
            this%normal_j1, &
            this%normal_jnyp, &
            this%normal_k1, &
            this%normal_knzp, &
            this%edge_vector_x, &
            this%edge_vector_y, &
            this%edge_vector_z, &
            this%plane_const_i1, &
            this%plane_const_inxp, &
            this%plane_const_j1, &
            this%plane_const_jnyp, &
            this%plane_const_k1, &
            this%plane_const_knzp, &
            this%is_wall_x_top, &
            this%is_wall_x_bot, &
            this%is_wall_y_top, &
            this%is_wall_y_bot, &
            this%is_wall_z_top, &
            this%is_wall_z_bot, &
            this%is_parallel, &
            '###'

        write(*,*) '@@@ end Write_boundary_parameters @@@'
#endif

    end subroutine Write_boundary_parameters

    subroutine Read_boundary_parameters(this, unit, iostat, iomsg)
        class (boundary_parameters_t), intent(in out) :: this  
        integer,      intent(in)     :: unit
        integer,      intent(out)    :: iostat
        character(*), intent(in out)  :: iomsg
        integer :: size_edge_num
        integer, dimension(2) :: shape_boundary_indices_map
        integer :: size_angles
        integer :: size_edge_vector_x
        integer :: size_edge_vector_y
        integer :: size_edge_vector_z

#ifdef DEBUG
        write(*,*) "@@@ in Read_boundary_parameters @@@"
#endif

        read(unit, iostat=iostat, iomsg=iomsg) &
            size_edge_num, &
            size_angles, &
            size_edge_vector_x, &
            size_edge_vector_y, &
            size_edge_vector_z

        deallocate(this%edge_num)
        allocate(this%edge_num(size_edge_num))
        deallocate(this%angles)
        allocate(this%angles(size_angles))
        deallocate(this%edge_vector_x)
        allocate(this%edge_vector_x(size_edge_vector_x))
        deallocate(this%edge_vector_y)
        allocate(this%edge_vector_y(size_edge_vector_y))
        deallocate(this%edge_vector_z)
        allocate(this%edge_vector_z(size_edge_vector_z))

        read(unit, iostat=iostat, iomsg=iomsg) &
            this%nxp, &
            this%nyp, &
            this%nzp, &
            this%dimension, &
            this%edge_num, &
            this%angles, &
            this%normal_i1, &
            this%normal_inxp, &
            this%normal_j1, &
            this%normal_jnyp, &
            this%normal_k1, &
            this%normal_knzp, &
            this%edge_vector_x, &
            this%edge_vector_y, &
            this%edge_vector_z, &
            this%plane_const_i1, &
            this%plane_const_inxp, &
            this%plane_const_j1, &
            this%plane_const_jnyp, &
            this%plane_const_k1, &
            this%plane_const_knzp, &
            this%is_wall_x_top, &
            this%is_wall_x_bot, &
            this%is_wall_y_top, &
            this%is_wall_y_bot, &
            this%is_wall_z_top, &
            this%is_wall_z_bot, &
            this%is_parallel


#ifdef DEBUG
        write(*,*) &
            size(this%edge_num), &
            size(this%angles), &
            size(this%edge_vector_x), &
            size(this%edge_vector_y), &
            size(this%edge_vector_z)
        write(*,*) &
            this%nxp, &
            this%nyp, &
            this%nzp, &
            this%dimension, &
            this%edge_num, &
            this%angles, &
            this%normal_i1, &
            this%normal_inxp, &
            this%normal_j1, &
            this%normal_jnyp, &
            this%normal_k1, &
            this%normal_knzp, &
            this%edge_vector_x, &
            this%edge_vector_y, &
            this%edge_vector_z, &
            this%plane_const_i1, &
            this%plane_const_inxp, &
            this%plane_const_j1, &
            this%plane_const_jnyp, &
            this%plane_const_k1, &
            this%plane_const_knzp, &
            this%is_wall_x_top, &
            this%is_wall_x_bot, &
            this%is_wall_y_top, &
            this%is_wall_y_bot, &
            this%is_wall_z_top, &
            this%is_wall_z_bot, &
            this%is_parallel, &
            '###'

        write(*,*) "@@@ end Read_boundary_parameters @@@"
#endif

    end subroutine Read_boundary_parameters


end module boundary_parameters_module
