
module materials_in_cells_module
    use cell_quantity_module, only : cell_quantity_t
    use cell_boundary_condition_module, only : cell_boundary_condition_t, cell_bc_wrapper_t
    use boundary_parameters_module      , only : boundary_parameters_t

    implicit none
    private
    public :: materials_in_cells_t

    type, extends(cell_quantity_t) :: materials_in_cells_t
        private


    contains

        procedure, public :: Clean_materials_in_cells

        procedure, public :: Write_materials_in_cells
        procedure, public :: Write_quantity_abstract => Write_materials_in_cells

        procedure, public :: Read_materials_in_cells
        procedure, public :: Read_quantity_abstract => Read_materials_in_cells

    end type


    interface materials_in_cells_t

        module procedure Constructor_2d
        module procedure Constructor_3d_sedov_taylor
        module procedure Constructor_3d_pyramid
        module procedure Constructor_3d_sod

    end interface materials_in_cells_t

contains

    type(materials_in_cells_t) function Constructor_2d(d1, d2, bc, bc_params, nc, number_layers_i, number_layers_j, &
        number_cells_i, number_cells_j, mat_index)
        integer                                           , intent(in) :: d1           
        integer                                           , intent(in) :: d2           
        type(cell_bc_wrapper_t), dimension(:), pointer    , intent(in) :: bc           
        integer                                           , intent(in) :: nc
        integer                                           , intent(in) :: number_layers_i  
        integer                                           , intent(in) :: number_layers_j  
        integer, dimension(:), allocatable                , intent(in) :: number_cells_i 
        integer, dimension(:), allocatable                , intent(in) :: number_cells_j 
        type(boundary_parameters_t), pointer, intent(in) :: bc_params

        integer, dimension(:), allocatable                , intent(in) :: mat_index

        real(8), dimension (d1, d2, 1)                                  :: init_values  
        integer :: i_curr, j_curr, num_mat, lay_j, lay_i, i, j, ele_i, ele_j
        j_curr = 1
        i_curr = 1
        num_mat = 1
        do lay_j = 1, number_layers_j
            i_curr = 1
            do lay_i = 1, number_layers_i
                ele_j = number_cells_j(lay_j) + j_curr 
                ele_i = number_cells_i(lay_i) + i_curr 
                do j = j_curr, ele_j
                    do i = i_curr, ele_i
                        init_values(i, j, 1) = mat_index(num_mat)
                    end do
                end do
                i_curr = i_curr + number_cells_i(lay_i) 
                num_mat = num_mat + 1
            end do
            j_curr = j_curr + number_cells_j(lay_j)
        end do
        call Constructor_2d%Init_cell_quantity_init_arr (init_values, d1, d2, 1, bc, bc_params)

    end function


    type(materials_in_cells_t) function Constructor_3d_sod(d1, d2, d3, bc,bc_params, nc, mat_index, parallel_params)
        use parallel_parameters_module ,only : parallel_parameters_t
        implicit none
        type(parallel_parameters_t), pointer              , intent(inout) :: parallel_params
        integer                                           , intent(in) :: d1           
        integer                                           , intent(in) :: d2           
        integer                                           , intent(in) :: d3           
        type(cell_bc_wrapper_t), dimension(:), pointer    , intent(in) :: bc           
        integer                                           , intent(in) :: nc
        integer, dimension(:), allocatable                , intent(in) :: mat_index
        type(boundary_parameters_t), pointer              , intent(in) :: bc_params

        real(8), allocatable, dimension (:,:,:)                        :: init_values  
        integer :: i_curr, j_curr, k_curr, num_mat, lay_k, lay_j, lay_i, i, j, k, ele_i, ele_j, ele_k
        integer :: i_from, i_to, j_from, j_to,k_from, k_to

        allocate(init_values(d1,d2,d3))
        init_values = 0d0

        k_curr = 1
        j_curr = 1
        i_curr = 1
        num_mat = 1
write(*,*) " making the sod"
        call parallel_params%get_virtual_from_to(i_from, i_to, 1)
        call parallel_params%get_virtual_from_to(j_from, j_to, 2)
        call parallel_params%get_virtual_from_to(k_from, k_to, 3)
        i_from = i_from + 1
        j_from = j_from + 1
        k_from = k_from + 1

        i_to = i_to - 1
        j_to = j_to - 1
        k_to = k_to - 1
        DO K=1, d3 - 1
            DO J=1, d2 - 1
                DO I=1, d1 - 1
                    if (i+j+k < (d3-1)+2+(d3-1)/2) then
                        init_values(i, j, k) = 1
                    else if (i+j+k > (d3-1)+2+(d3-1)/2) then
                        init_values(i, j, k) = 2
                    else
                        init_values(i, j, k) = 3
                    end if
                END DO
            END DO
        END DO

        call Constructor_3d_sod%Init_cell_quantity_init_arr (init_values, d1, d2, d3, bc, bc_params)
        deallocate(init_values)
    end function

    type(materials_in_cells_t) function Constructor_3d_pyramid(d1, d2, d3, bc,bc_params&
                                                    , number_layers, number_cells, mat_index, parallel_params)
        use parallel_parameters_module , only : parallel_parameters_t
        use geometry_module            , only : Create_start_index_layer
        implicit none
        type(parallel_parameters_t), pointer              , intent(inout) :: parallel_params
        integer                                           , intent(in) :: d1           
        integer                                           , intent(in) :: d2           
        integer                                           , intent(in) :: d3           
        type(cell_bc_wrapper_t), dimension(:), pointer    , intent(in) :: bc           
        integer                                           , intent(in) :: number_layers 
        integer, dimension(:), allocatable                , intent(in) :: number_cells
        type(boundary_parameters_t), pointer              , intent(in) :: bc_params
        integer, dimension(:), allocatable                , intent(in) :: mat_index

        real(8), allocatable, dimension (:,:,:)                        :: init_values  
        integer :: i,j,k,layer_index, counter, k_from, k_to
        integer, dimension(:, :), allocatable :: start_index

        allocate(init_values(d1,d2,d3))
        init_values = 0d0

        call parallel_params%get_virtual_from_to(k_from, k_to, 3)

        k_from = k_from + 1
        k_to = k_to - 1

        start_index = Create_start_index_layer(number_cells, number_layers, 1)
        do layer_index = 1, size(start_index(:, 1)) - 1 
            if (start_index(layer_index + 1,1) > k_to .or. start_index(layer_index,1) < k_from) cycle 
            do k = start_index(layer_index,1), start_index(layer_index + 1,1)
                do j = 1, d2
                    do i = 1, d1
                        init_values(i, j, k) = mat_index(layer_index)
                    end do
                end do
            end do
        end do

        call Constructor_3d_pyramid%Init_cell_quantity_init_arr (init_values, d1, d2, d3, bc, bc_params)
        deallocate(init_values)
        deallocate(start_index)
    end function

    type(materials_in_cells_t) function Constructor_3d_sedov_taylor(d1, d2, d3, bc,bc_params, parallel_params)
        use parallel_parameters_module ,only : parallel_parameters_t
        implicit none
        type(parallel_parameters_t), pointer              , intent(inout) :: parallel_params
        integer                                           , intent(in) :: d1           
        integer                                           , intent(in) :: d2           
        integer                                           , intent(in) :: d3           
        type(cell_bc_wrapper_t), dimension(:), pointer    , intent(in) :: bc           
        type(boundary_parameters_t)          , pointer    , intent(in) :: bc_params


        real(8), allocatable, dimension (:,:,:)                        :: init_values  
        integer :: i_curr, j_curr, k_curr, num_mat, lay_k, lay_j, lay_i, i, j, k, ele_i, ele_j, ele_k
        integer :: i_from, i_to, j_from, j_to,k_from, k_to

        allocate(init_values(d1,d2,d3))
        init_values = 0d0

        k_curr = 1
        j_curr = 1
        i_curr = 1
        num_mat = 1

        call parallel_params%get_virtual_from_to(i_from, i_to, 1)
        call parallel_params%get_virtual_from_to(j_from, j_to, 2)
        call parallel_params%get_virtual_from_to(k_from, k_to, 3)
        i_from = i_from + 1
        j_from = j_from + 1
        k_from = k_from + 1

        i_to = i_to - 1
        j_to = j_to - 1
        k_to = k_to - 1

        write(*,*) " making sedov-taylor"



        DO K=1, d3 - 1
            DO J=1, d2 - 1
                DO I=1, d1 - 1
                        if (parallel_params%i_virt(i) == 1 .and. &
                            parallel_params%j_virt(j) == 1 .and. &
                            parallel_params%k_virt(k) == 1) then
                            !write(*,*), "in section 1",i,j,k
                            init_values(i, j, k) = 1
                        else
                            init_values(i, j, k) = 2
                        end if
                    !write(*,*),"init_values: ",i,j,k, init_values(i,j,k)
                END DO
            END DO
        END DO


        call Constructor_3d_sedov_taylor%Init_cell_quantity_init_arr (init_values, d1, d2, d3, bc, bc_params)
        deallocate(init_values)
    end function


    subroutine Clean_materials_in_cells (this)
        class (materials_in_cells_t), intent(in out) :: this 

        call this%Clean_cell_quantity ()
    end subroutine Clean_materials_in_cells

    subroutine Write_materials_in_cells(this, unit, iostat, iomsg)
        class (materials_in_cells_t), intent(in) :: this  
        integer,      intent(in)     :: unit
        integer,      intent(out)    :: iostat
        character(*), intent(in out)  :: iomsg

#ifdef DEBUG
        write(*,*) '@@@ in Write_materials_in_cells @@@'
#endif

        call this%Write_cell_quantity(unit, iostat=iostat, iomsg=iomsg)

#ifdef DEBUG
        write(*,*) '@@@ end Write_materials_in_cells @@@'
#endif

    end subroutine Write_materials_in_cells

    subroutine Read_materials_in_cells(this, unit, iostat, iomsg)
        class (materials_in_cells_t), intent(in out) :: this  
        integer,      intent(in)     :: unit
        integer,      intent(out)    :: iostat
        character(*), intent(in out)  :: iomsg

#ifdef DEBUG
        write(*,*) '@@@ in Read_materials_in_cells @@@'
#endif

        call this%Read_cell_quantity(unit, iostat=iostat, iomsg=iomsg)

#ifdef DEBUG
        write(*,*) '@@@ end Read_materials_in_cells @@@'
#endif

    end subroutine Read_materials_in_cells

end module materials_in_cells_module

