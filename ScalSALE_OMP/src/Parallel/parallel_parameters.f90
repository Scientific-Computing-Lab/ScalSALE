
module parallel_parameters_module

    implicit none
    private
    public :: parallel_parameters_t

    type :: parallel_parameters_t
        private
        integer, public               :: virt_nxp   
        integer, public               :: virt_nyp   
        integer, public               :: virt_nzp   

        integer, public               :: virt_nx    
        integer, public               :: virt_ny    
        integer, public               :: virt_nz    

        integer, public               :: nxp        
        integer, public               :: nyp        
        integer, public               :: nzp        

        integer, public               :: nx         
        integer, public               :: ny         
        integer, public               :: nz         

        integer, public               :: np         
        integer, public               :: npx        
        integer, public               :: npy        
        integer, public               :: npz        

        integer, dimension(3), public :: my_coords  

        integer, dimension(:),pointer, public :: i_virt  
        integer, dimension(:),pointer, public :: j_virt  
        integer, dimension(:),pointer, public :: k_virt  

        integer, public               :: my_id      
        integer, public               :: my_rank      

        logical, public               :: is_wall_x_top    
        logical, public               :: is_wall_x_bot    
        logical, public               :: is_wall_y_top    
        logical, public               :: is_wall_y_bot    
        logical, public               :: is_wall_z_top    
        logical, public               :: is_wall_z_bot    

        logical, public               :: is_parallel
    contains
        procedure, public :: get_virtual_array
        procedure, public :: get_physical_array
        procedure, public :: get_virtual_index
        procedure, public :: get_virtual_index_i
        procedure, public :: get_virtual_index_j
        procedure, public :: get_virtual_index_k
        procedure, public :: get_physical_index
        procedure, public :: get_virtual_from_to
        procedure, public :: Point_to_virtual_array
        procedure, public :: Inside_physical_world
    end type parallel_parameters_t

    interface parallel_parameters_t
        module procedure Constructor
    end interface parallel_parameters_t

contains

    type(parallel_parameters_t) function Constructor(my_id, my_coords, np, npx, npy, npz, virt_nxp, virt_nyp, virt_nzp&
        , nxp, nyp, nzp, nx&
        , ny, nz)
        implicit none
        integer, intent(in)                 :: my_id
        integer, intent(in)                 :: np
        integer, intent(in)                 :: npx, npy, npz
        integer, intent(in)                 :: virt_nxp, virt_nyp, virt_nzp
        integer, intent(in)                 :: nxp, nyp, nzp
        integer, intent(in)                 :: nx, ny, nz
        integer, dimension(3), intent(in)   :: my_coords

        Constructor%my_id = my_id
        Constructor%my_rank = my_id

        Constructor%np = np
        Constructor%npx = npx
        Constructor%npy = npy
        Constructor%npz = npz

        Constructor%virt_nxp = virt_nxp
        Constructor%virt_nyp = virt_nyp
        Constructor%virt_nzp = virt_nzp

        Constructor%virt_nx = virt_nxp - 1
        Constructor%virt_ny = virt_nyp - 1
        Constructor%virt_nz = virt_nzp - 1
        if (Constructor%virt_nzp == 1) Constructor%virt_nz = 1

        Constructor%nxp = nxp
        Constructor%nyp = nyp
        Constructor%nzp = nzp

        Constructor%nx = nx
        Constructor%ny = ny
        Constructor%nz = nz

        Constructor%my_coords = my_coords

        Constructor%is_wall_x_top = .false.
        Constructor%is_wall_x_bot = .false.
        Constructor%is_wall_y_top = .false.
        Constructor%is_wall_y_bot = .false.
        Constructor%is_wall_z_top = .false.
        Constructor%is_wall_z_bot = .false.
        if (my_coords(1) == npx) then
            Constructor%is_wall_x_top = .true.
        end if
        if (my_coords(1) == 1) then
            Constructor%is_wall_x_bot = .true.
        end if

        if (my_coords(2) == npy) then
            Constructor%is_wall_y_top = .true.
        end if
        if (my_coords(2) == 1) then
            Constructor%is_wall_y_bot = .true.
        end if

        if (my_coords(3) == npz) then
            Constructor%is_wall_z_top = .true.
        end if
        if (my_coords(3) == 1) then
            Constructor%is_wall_z_bot = .true.
        end if

        if (virt_nzp == 1) then
            Constructor%is_wall_z_bot = .false.
            Constructor%is_wall_z_top = .false.
        end if

        Constructor%is_parallel = .true.
        if (virt_nxp == nxp .and. virt_nyp == nyp .and. virt_nzp == nzp) then
            Constructor%is_parallel = .false.
        end if

        allocate(Constructor%i_virt(0:Constructor%nxp + 1))
        allocate(Constructor%j_virt(0:Constructor%nyp + 1))
        allocate(Constructor%k_virt(0:Constructor%nzp + 1))
        call Constructor%get_virtual_array(1,Constructor%i_virt)
        call Constructor%get_virtual_array(2,Constructor%j_virt)
        call Constructor%get_virtual_array(3,Constructor%k_virt)
    end function

    subroutine get_virtual_array(this, ind, virt_array)
        class(parallel_parameters_t), intent(inout)      :: this             
        integer, dimension(:), pointer, intent(out)      :: virt_array
        integer,  intent(in)                             :: ind

        integer :: i, stride, to

        if (ind == 1) then
            stride = this%nx
            to = this%nxp + 1
        elseif (ind == 2) then
            stride = this%ny
            to = this%nyp + 1
        else
            stride = this%nz
            to = this%nzp + 1
        end if

        allocate(virt_array(0:to))

        do i = 0, to
            virt_array(i) = (this%my_coords(ind) - 1) * stride + i
        end do
    end subroutine get_virtual_array

    function get_virtual_index(this, i ,ind) result (virt_index)
        class(parallel_parameters_t), intent(inout)      :: this             
        integer,  intent(in)                             :: ind, i
        integer :: virt_index
        integer :: stride

        if (ind == 1) then
            stride = this%nx
        elseif (ind == 2) then
            stride = this%ny
        else
            stride = this%nz
        end if

        virt_index = (this%my_coords(ind) - 1) * stride + i
    end function get_virtual_index

    function get_physical_index(this, i ,ind) result (phys_index)
        class(parallel_parameters_t), intent(inout)      :: this             
        integer,  intent(in)                             :: ind, i
        integer :: phys_index
        integer :: stride, plus

        if (ind == 1) then
            stride = this%nx
            plus = this%nxp
        elseif (ind == 2) then
            stride = this%ny
            plus = this%nyp
        else
            stride = this%nz
            plus = this%nzp
        end if

        if ( (i > ( (this%my_coords(ind) - 1) * stride)).and.(i <= (this%my_coords(ind) * stride))) then
            phys_index = mod(i, stride)
            if (phys_index == 0) phys_index = stride
        elseif (i == ((this%my_coords(ind) - 1) * stride)) then
            phys_index = 0
        elseif (i == ((this%my_coords(ind)) * stride + 1)) then
            phys_index = plus
        else
            phys_index = -1
        end if

    end function get_physical_index

    subroutine get_physical_array(this, ind, phys_array)
        implicit none
        class(parallel_parameters_t), intent(inout)      :: this             
        integer, dimension(:), pointer, intent(out)      :: phys_array
        integer,  intent(in)                             :: ind

        integer :: i, stride, to
        integer :: plus

        if (ind == 1) then
            stride = this%nx
            to = this%virt_nxp
            plus = this%nxp
        elseif (ind == 2) then
            stride = this%ny
            to = this%virt_nyp
            plus = this%nyp
        else
            stride = this%nz
            to = this%virt_nzp
            plus = this%nzp
        end if

        do i = 0, to
            if (i > ( (this%my_coords(ind) - 1) * stride) .and. i <= (this%my_coords(ind) * stride)) then
                phys_array(i) = mod(i, stride)
                if (phys_array(i) == 0) phys_array(i) = stride
            elseif (i == ((this%my_coords(ind) - 1) * stride)) then
                phys_array(i) = 0
            elseif (i == ((this%my_coords(ind)) * stride + 1)) then
                phys_array(i) = plus
            else
                phys_array(i) = -1
            end if
        end do
    end subroutine get_physical_array


    subroutine get_virtual_from_to(this, from, to, ind)
        implicit none
        class(parallel_parameters_t), intent(inout)      :: this             
        integer, intent(out) :: from, to
        integer, intent(in)  :: ind

        integer :: stride, to2
        if (ind == 1) then
            stride = this%nx
            to2 = this%nxp + 1
        elseif (ind == 2) then
            stride = this%ny
            to2 = this%nyp + 1
        else
            stride = this%nz
            to2 = this%nzp + 1
        end if

        from = (this%my_coords(ind) - 1) * stride
        to = (this%my_coords(ind) - 1) * stride + to2

    end subroutine get_virtual_from_to

    function get_virtual_index_i(this, i) result (virt_index)
        class(parallel_parameters_t), intent(inout)      :: this             
        integer,  intent(in)                             ::  i
        integer :: virt_index
        integer :: stride
        virt_index = (this%my_coords(1) - 1) * this%nx + i
    end function get_virtual_index_i

    function get_virtual_index_j(this, i) result (virt_index)
        class(parallel_parameters_t), intent(inout)      :: this             
        integer,  intent(in)                             ::  i
        integer :: virt_index

        virt_index = (this%my_coords(2) - 1) * this%ny + i

    end function get_virtual_index_j

    function get_virtual_index_k(this, i) result (virt_index)
        class(parallel_parameters_t), intent(inout)      :: this             
        integer,  intent(in)                             ::  i
        integer :: virt_index
        integer :: stride

        virt_index = (this%my_coords(3) - 1) * this%nz + i
    end function get_virtual_index_k

    subroutine Point_to_virtual_array(this, ptr_x, ptr_y, ptr_z)
        class(parallel_parameters_t), intent(inout)      :: this             
        integer       , dimension(:), pointer, intent(out)    :: ptr_x   
        integer       , dimension(:), pointer, intent(out)    :: ptr_y   
        integer       , dimension(:), pointer, intent(out)    :: ptr_z   

        ptr_x => this%i_virt
        ptr_y => this%j_virt
        ptr_z => this%k_virt

    end subroutine Point_to_virtual_array

    function Inside_physical_world(this, ind, axis) result(inside)
        implicit none
        class(parallel_parameters_t), intent(inout)      :: this             
        integer, intent(in)                              :: ind 
        integer, intent(in)                              :: axis 
        integer :: inside

        integer :: start, end

        call this%get_virtual_from_to(start, end, axis)

        if (ind >= start .and. ind <= end) then
            inside = 0
        end if

        if (ind < start) then
            inside = -1
        end if
        if (ind > start .and. ind > end ) then
            inside = 1
        end if
    end function Inside_physical_world
end module parallel_parameters_module
