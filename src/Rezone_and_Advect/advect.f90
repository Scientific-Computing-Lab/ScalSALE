
module advect_module
    use vertex_boundary_condition_module, only : vertex_boundary_condition_t, vertex_bc_wrapper_t
    use cell_boundary_condition_module  , only : cell_boundary_condition_t  , cell_bc_wrapper_t
    use num_materials_in_cells_module   , only : num_materials_in_cells_t
    use materials_in_cells_module       , only : materials_in_cells_t
    use material_advect_module          , only : material_advect_t
    use vertex_mass_module              , only : vertex_mass_t
    use cell_mass_module                , only : cell_mass_t
    use velocity_module                 , only : velocity_t
    use momentum_module                 , only : momentum_t
    use material_module                 , only : material_t
    use density_module                  , only : density_t
    use energy_module                   , only : energy_t
    use rezone_module                   , only : rezone_t
    use volume_module                   , only : volume_t
    use vof_module                      , only : vof_t
    use mesh_base_module                , only : mesh_base_t
    use data_module                     , only : data_t
    use parallel_parameters_module      , only : parallel_parameters_t
    use communication_parameters_module , only : communication_parameters_t
    use communication_module            , only : communication_t
    use material_quantity_module        , only : material_quantity_t
    implicit none
    private
    public :: advect_t

    type :: advect_t
        private


        type (rezone_t), pointer :: rezone

        type (material_t)       , pointer :: materials
        type (material_advect_t), pointer :: adv_mats

        type (data_t)                                   :: area_top_out
        type (data_t)                                   :: area_top_in

        type (material_quantity_t), pointer                          :: a
        type (material_quantity_t), pointer                          :: b
        type (material_quantity_t), pointer                          :: c
        type (material_quantity_t), pointer                          :: side

        type (vertex_mass_t)                            :: adv_vertex_mass
        type (cell_mass_t)                              :: adv_cell_mass
        type (energy_t)                                 :: adv_sie
        type (vof_t)                                    :: vof_advect
        type (parallel_parameters_t), pointer           :: parallel_params

        type (num_materials_in_cells_t)       , pointer :: num_mat_cells
        type (materials_in_cells_t)           , pointer :: mat_cells
        type (vertex_mass_t)                  , pointer :: vertex_mass
        type (cell_mass_t)                    , pointer :: total_cell_mass
        type (density_t)                      , pointer :: total_density
        type (energy_t)                       , pointer :: total_sie
        type (vof_t)                          , pointer :: total_vof

        type (velocity_t)                     , pointer :: velocity
        type (volume_t)                       , pointer :: volume
        class(mesh_base_t)                    , pointer :: mesh

        type (momentum_t)             , pointer, public :: momentum
        type (momentum_t)             , pointer, public :: adv_momentum

        type (data_t)                                   :: vel_adv_weight

        integer                                         :: advect_init_layer_mat
        integer                                         :: a0
        integer                                         :: b0
        integer                                         :: nx
        integer                                         :: ny
        integer                                         :: nz
        integer                                         :: nxp
        integer                                         :: nyp
        integer                                         :: nzp
        integer                                         :: n_materials
        integer                                         :: cyc_delete
        logical                                         :: shorter_advect
        logical                                         :: fix_overflow
        logical                                         :: old_volmat
        logical                                         :: old_line_calc
        integer                                         :: wilkins_scheme
        real(8)                                         :: total_out
        real(8)                                         :: total_in
        real(8)                                         :: total_left_out
        real(8)                                         :: total_left_in
        real(8)                                         :: total_right_out
        real(8)                                         :: total_right_in
        real(8)                                         :: total_bottom_out
        real(8)                                         :: total_bottom_in
        real(8)                                         :: total_top_out
        real(8)                                         :: total_top_in
        real(8)                                         :: emf
        real(8)                                         :: emfm
        real(8)                                         :: emf1

    contains


        procedure, public :: Clean_advect

        procedure, public :: Calculate_advect_2d

        procedure, public :: Calculate_advect_3d

        procedure, public :: Calculate_areas_for_advect

        procedure, public :: Calc_area_for_each_side

        procedure, public :: Volume_material_new

        procedure, public :: Volume_material_old

        procedure, public :: Volume_material_3d

        procedure, public :: Find_donnor_acceptor

        procedure, public :: Put_quantities_into_cell

        procedure, public :: Calculate_cell_quantities_in_advect

        procedure, public :: Advect_velocities

        procedure, public :: Fix_overflow_error

        procedure, public :: Fix_overflow_for_passed_cell

        procedure, public :: Calculate_vof_line_new

        procedure, private :: Line_calc_3d

        procedure, public :: Point_to_material

        procedure, public :: Point_to_adv_material

        procedure, public :: Point_to_adv_vof_data

        procedure, public :: Point_to_adv_cell_mass_data

        procedure, public :: Point_to_adv_sie_data

        procedure, public :: Point_to_momentum_data

        procedure, public :: Point_to_adv_momentum_data

        procedure, public :: Point_to_adv_vertex_mass_data

        procedure, public :: Set_advect_init_layer_mat

        procedure, public :: Get_all_areas

        procedure, public :: Set_communication
      
        procedure, public :: Write_advect
        generic :: write(unformatted) => Write_advect


        procedure, public :: Read_advect
        generic :: read(unformatted) => Read_advect
    end type


    interface advect_t
        module procedure Constructor
    end interface advect_t

contains

    type(advect_t) function Constructor(nxp, nyp, nzp, nmats, mat_ids, bcc, bcv, bc_params,&
        line_calc,shorter_advect, fix_overflow, rezone, mesh, &
        materials, mat_cells, num_mat_cells, total_sie, total_vof, total_density, total_cell_mass, &
        vertex_mass, volume, velocity, emf, emfm, wilkins_scheme, parallel_params)
        use boundary_parameters_module      , only : boundary_parameters_t

        implicit none
        integer                                 , intent(in) :: nxp
        integer                                 , intent(in) :: nyp
        integer                                 , intent(in) :: nzp
        logical                                 , intent(in) :: shorter_advect
        logical                                 , intent(in) :: fix_overflow
        logical                                 , intent(in) :: line_calc
        type (rezone_t)                , pointer, intent(in) :: rezone
        type (energy_t)                , pointer, intent(in) :: total_sie
        type (density_t)               , pointer, intent(in) :: total_density
        type (vof_t)                   , pointer, intent(in) :: total_vof
        type (cell_mass_t)             , pointer, intent(in) :: total_cell_mass
        type (vertex_mass_t)           , pointer, intent(in) :: vertex_mass
        type (volume_t)                , pointer, intent(in) :: volume
        type (velocity_t)              , pointer, intent(in) :: velocity
        class(mesh_base_t)             , pointer, intent(in) :: mesh
        type (materials_in_cells_t)    , pointer, intent(in) :: mat_cells
        type (num_materials_in_cells_t), pointer, intent(in) :: num_mat_cells
        type (parallel_parameters_t)   , pointer, intent(in) :: parallel_params
        type (material_t),  pointer, intent(inout) :: materials
        type(boundary_parameters_t)    , pointer, intent(in) :: bc_params
        integer                               , intent(in)           :: nmats
        integer,dimension(:), allocatable         , intent(in)           :: mat_ids
        type (cell_bc_wrapper_t  ), dimension(:), pointer, intent(in) :: bcc
        type (vertex_bc_wrapper_t), dimension(:), pointer, intent(in) :: bcv
        real(8)                                          , intent(in) :: emf
        real(8)                                          , intent(in) :: emfm

        integer, intent(in) :: wilkins_scheme

        integer :: i

        allocate(Constructor%momentum)
        allocate(Constructor%adv_momentum)
        allocate(Constructor%a)
        allocate(Constructor%b)
        allocate(Constructor%c)
        allocate(Constructor%side)

        Constructor%shorter_advect = shorter_advect
        Constructor%fix_overflow   = fix_overflow
        Constructor%old_line_calc    = line_calc
        Constructor%wilkins_scheme = wilkins_scheme
        Constructor%old_volmat     = .true.
        Constructor%a0             = 1d0
        Constructor%b0             = 0d0
        Constructor%adv_cell_mass = cell_mass_t  (0d0, nxp, nyp, nzp, bcc,bc_params)
        if (nzp == 1) then
            Constructor%adv_vertex_mass  = vertex_mass_t(0d0, nxp+1, nyp+1, 1, bcv,bc_params)
            Constructor%momentum     = momentum_t(0d0, nxp+1, nyp+1, 1, bcv, bc_params)
            Constructor%adv_momentum = momentum_t(0d0, nxp+1, nyp+1, 1, bcv, bc_params)
            Constructor%vel_adv_weight = data_t(0d0, nxp+1, nyp+1, 1)

        else
            Constructor%vel_adv_weight = data_t(0d0, nxp+1, nyp+1, nzp+1)
            Constructor%adv_vertex_mass  = vertex_mass_t(0d0, nxp+1, nyp+1, nzp+1, bcv,bc_params)
            Constructor%momentum     = momentum_t(0d0, nxp+1, nyp+1, nzp+1, bcv, bc_params)
            Constructor%adv_momentum = momentum_t(0d0, nxp+1, nyp+1, nzp+1, bcv, bc_params)
        end if

        Constructor%vof_advect  = vof_t   (0d0, nxp, nyp, nzp, bcc,bc_params)
        Constructor%adv_sie     = energy_t(0d0, nxp, nyp, nzp, bcc,bc_params)

        Constructor%n_materials = nmats
        Constructor%nx  = nxp - 1
        Constructor%ny  = nyp - 1
        if (nzp == 1) then
            Constructor%nz  = 1
        else
            Constructor%nz  = nzp - 1
            Constructor%a = material_quantity_t(0d0, nxp, nyp, nzp, 1)
            Constructor%b = material_quantity_t(0d0, nxp, nyp, nzp, 1)
            Constructor%c = material_quantity_t(0d0, nxp, nyp, nzp, 1)
            Constructor%side = material_quantity_t(0d0, nxp, nyp, nzp,1 )
        end if

        Constructor%nxp = nxp
        Constructor%nyp = nyp
        Constructor%nzp = nzp
        Constructor%emf  = emf
        Constructor%emf1  = 1 - emf
        Constructor%emfm = emfm
        Constructor%cyc_delete = 0
        allocate(material_advect_t :: Constructor%adv_mats)


        Constructor%mesh => mesh
        Constructor%total_sie => total_sie
        Constructor%total_density => total_density
        Constructor%total_vof => total_vof
        Constructor%total_cell_mass => total_cell_mass
        Constructor%vertex_mass=> vertex_mass
        Constructor%volume => volume
        Constructor%velocity => velocity
        Constructor%rezone => rezone
        Constructor%mat_cells => mat_cells
        Constructor%num_mat_cells => num_mat_cells
        Constructor%materials => materials
        Constructor%parallel_params => parallel_params

        Constructor%adv_mats = material_advect_t(nxp, nyp, nzp, nmats, mat_ids, bcc, bc_params)



    end function

    subroutine Calculate_advect_2d(this)
        use vof_module, only : vof_t

        class (advect_t)    , intent(in out) :: this


        !        real(8), dimension(:, :, :), pointer :: area_top_out
        !        real(8), dimension(:, :, :), pointer :: area_top_in
        real(8), dimension(:, :, :), pointer :: vof
        real(8), dimension(:, :, :), pointer :: cell_mass_adv
        real(8), dimension(:, :, :), pointer :: vof_adv
        real(8), dimension(:, :, :), pointer :: sie_adv
        real(8), dimension(:, :, :, :), pointer :: mat_vof_adv
        real(8), dimension(:, :, :), pointer :: mat_vof_delete
        real(8), dimension(:, :, :), pointer :: sie_vof_delete
        real(8), dimension(:,:,:,:), pointer :: area_top_in,area_top_out
        integer :: i, j, i1, j1, nx, ny, tmp_mat


        this%cyc_delete=this%cyc_delete+1
        call this%total_vof       %Point_to_data(vof)
        call this%adv_cell_mass%Point_to_data(cell_mass_adv)
        call this%vof_advect      %Point_to_data(vof_adv)
        call this%adv_sie         %Point_to_data(sie_adv)
        call this%adv_mats%vof%Point_to_data(mat_vof_adv)
        call this%adv_mats%area_top_in%Point_to_data(area_top_in)
        call this%adv_mats%area_top_out%Point_to_data(area_top_out)
        nx = this%nx
        ny = this%ny


        sie_adv = 0d0
        cell_mass_adv = 0d0
        vof_adv = 0d0

        area_top_in = 0d0
        area_top_out = 0d0

        this%adv_mats%area_right_out = 0d0
        this%adv_mats%area_right_in = 0d0
        this%adv_mats%area_left_out = 0d0
        this%adv_mats%area_left_in = 0d0
        this%adv_mats%area_bottom_out = 0d0
        this%adv_mats%area_bottom_in = 0d0
        mat_vof_adv = 0d0

        !        do tmp_mat = 1, this%n_materials
        !            call this%adv_mats(tmp_mat)%Point_to_areas(area_top_in, area_top_out)
        !            call this%adv_mats(tmp_mat)%vof%Point_to_data(mat_vof_adv)
        !            this%adv_mats(tmp_mat)%area_right_out  = 0d0
        !            this%adv_mats(tmp_mat)%area_right_in   = 0d0
        !            this%adv_mats(tmp_mat)%area_left_out   = 0d0
        !            this%adv_mats(tmp_mat)%area_left_in    = 0d0
        !            this%adv_mats(tmp_mat)%area_bottom_out = 0d0
        !            this%adv_mats(tmp_mat)%area_bottom_in  = 0d0
        !            area_top_in  = 0d0
        !            area_top_out = 0d0
        !            mat_vof_adv  = 0d0
        !        end do



        do j = 1, ny
            do i = 1, nx
                if ((vof(i, j    , 1) >= this%emf) .or. (vof(i + 1, j    , 1) >= this%emf) .or. (vof(i - 1, j, 1) >= this%emf) .or.&
                    (vof(i, j + 1, 1) >= this%emf) .or. (vof(i    , j - 1, 1) >= this%emf)) then
                    call this%Calculate_areas_for_advect(i, j)
                    call this%Calculate_cell_quantities_in_advect(i, j, 1)
                    call this%Fix_overflow_error(i, j, 1)

                end if
            end do
        end do


        call this%total_density%Exchange_virtual_space_blocking()
        call this%vertex_mass%Calculate_vertex_mass_2d(this%mesh%coordinates, this%total_density, this%total_cell_mass , &
            this%wilkins_scheme, this%mesh%cyl)
        call this%vertex_mass%Exchange_virtual_space_blocking()
        call this%adv_vertex_mass%Calculate_vertex_mass_2d  (this%mesh%coordinates, this%total_density, this%adv_cell_mass, &
            this%wilkins_scheme, this%mesh%cyl)
        call this%adv_vertex_mass%Exchange_virtual_space_blocking()

        do j = 1, ny
            do i = 1, nx
                call this%Put_quantities_into_cell(i, j, 1)
            end do
        end do
        call this%total_cell_mass%Point_to_data(mat_vof_delete)




        call this%Advect_velocities()



        return
    end subroutine Calculate_advect_2d

    subroutine Calculate_areas_for_advect(this, i, j)
        implicit none

        class (advect_t)                     :: this

        integer, intent(in) :: i, j
        real(8), dimension(:, :, :), pointer :: x
        real(8), dimension(:, :, :), pointer :: y
        real(8), dimension(:, :, :), pointer :: x_material
        real(8), dimension(:, :, :), pointer :: y_material
        real(8), dimension(:, :, :), pointer :: vol
        real(8), dimension(:, :, :, :), pointer :: area_top_out
        real(8), dimension(:, :, :, :), pointer :: area_top_in

        real(8), dimension(:), allocatable :: area_out_tmp
        real(8), dimension(:), allocatable :: area_in_tmp

        real(8) :: total_area_out_tmp
        real(8) :: total_area_in_tmp

        real(8) :: x1    , x2    , x3    , x4
        real(8) :: y1    , y2    , y3    , y4
        real(8) :: x_lag1, x_lag2, x_lag3, x_lag4
        real(8) :: y_lag1, y_lag2, y_lag3, y_lag4

        real(8) :: vol_right, vol_left
        real(8) :: vol_top, vol_bottom
        real(8) :: vol_center

        integer :: tmp_mat
        call this%adv_mats%area_top_in%Point_to_data(area_top_in)
        call this%adv_mats%area_top_out%Point_to_data(area_top_out)

        call this%rezone%Point_to_coordinates_2d(x_material, y_material)
        call this%mesh  %Point_to_data (x, y)
        call this%volume%Point_to_data (vol)
        allocate(area_out_tmp(this%n_materials))
        allocate(area_in_tmp(this%n_materials))

        area_out_tmp = 0d0
        area_in_tmp  = 0d0

        x1 = x(i+1, j  , 1)
        x2 = x(i+1, j+1, 1)
        x3 = x(i  , j+1, 1)
        x4 = x(i  , j  , 1)

        y1 = y(i+1, j  , 1)
        y2 = y(i+1, j+1, 1)
        y3 = y(i  , j+1, 1)
        y4 = y(i  , j  , 1)

        x_lag1 = x_material(i+1, j  , 1)
        x_lag2 = x_material(i+1, j+1, 1)
        x_lag3 = x_material(i  , j+1, 1)
        x_lag4 = x_material(i  , j  , 1)

        y_lag1 = y_material(i+1, j  , 1)
        y_lag2 = y_material(i+1, j+1, 1)
        y_lag3 = y_material(i  , j+1, 1)
        y_lag4 = y_material(i  , j  , 1)

        vol_left   = vol(i-1, j  , 1)
        vol_bottom = vol(i  , j-1, 1)
        vol_right  = vol(i+1, j  , 1)
        vol_top    = vol(i  , j+1, 1)
        vol_center = vol(i  , j  , 1)


        if (i == 1) then
            vol_left  = vol_center
        else if (i == this%nx) then
            vol_right = vol_center
        end if
        if (j == 1) then
            vol_bottom  = vol_center
        else if (j == this%ny) then
            vol_top     = vol_center
        end if


        if (this%shorter_advect .and. i /= 1) then
            this%total_left_out = 0d0
            this%total_left_in  = 0d0
            do tmp_mat = 1, this%n_materials
                this%adv_mats%area_left_out(tmp_mat) = -this%adv_mats%area_right_in(tmp_mat)
                this%adv_mats%area_left_in(tmp_mat)  = -this%adv_mats%area_right_out(tmp_mat)
                this%total_left_out = this%total_left_out + this%adv_mats%area_left_out(tmp_mat)
                this%total_left_in  = this%total_left_in  + this%adv_mats%area_left_in(tmp_mat)
            end do
        else
            call this%Calc_area_for_each_side(x3, y3, x4, y4, x_lag3, y_lag3, x_lag4, y_lag4, &
                vol_left, vol_center, .true. , -1, i, j       , &
                area_out_tmp, area_in_tmp, total_area_out_tmp , total_area_in_tmp)
            do tmp_mat = 1, this%n_materials
                this%adv_mats%area_left_out(tmp_mat) = area_out_tmp(tmp_mat)
                this%adv_mats%area_left_in(tmp_mat)  = area_in_tmp (tmp_mat)
            end do
            this%total_left_out = total_area_out_tmp
            this%total_left_in  = total_area_in_tmp
        end if


        if (this%shorter_advect .and. j /= 1) then
            this%total_bottom_out = 0d0
            this%total_bottom_in  = 0d0
            do tmp_mat = 1, this%n_materials
                !                call this%adv_mats%Point_to_areas(area_top_in, area_top_out)
                this%adv_mats%area_bottom_out(tmp_mat) = -area_top_in (tmp_mat, i, 1, 1)
                this%adv_mats%area_bottom_in(tmp_mat)  = -area_top_out(tmp_mat, i, 1, 1)
                this%total_bottom_out = this%total_bottom_out + this%adv_mats%area_bottom_out(tmp_mat)
                this%total_bottom_in  = this%total_bottom_in  + this%adv_mats%area_bottom_in(tmp_mat)
            end do
        else
            call this%Calc_area_for_each_side(x4, y4, x1, y1, x_lag4, y_lag4, x_lag1, y_lag1, &
                vol_bottom, vol_center, .false. , -1, i, j    , &
                area_out_tmp, area_in_tmp, total_area_out_tmp , total_area_in_tmp)
            do tmp_mat = 1, this%n_materials
                this%adv_mats%area_bottom_out(tmp_mat) = area_out_tmp(tmp_mat)
                this%adv_mats%area_bottom_in(tmp_mat)  = area_in_tmp (tmp_mat)
            end do
            this%total_bottom_out = total_area_out_tmp
            this%total_bottom_in  = total_area_in_tmp
        end if


        call this%Calc_area_for_each_side(x1, y1, x2, y2, x_lag1, y_lag1, x_lag2, y_lag2, &
            vol_right, vol_center, .true. , 1, i, j       , &
            area_out_tmp, area_in_tmp, total_area_out_tmp , total_area_in_tmp)
        do tmp_mat = 1, this%n_materials
            this%adv_mats%area_right_out(tmp_mat) = area_out_tmp(tmp_mat)
            this%adv_mats%area_right_in(tmp_mat)  = area_in_tmp (tmp_mat)
        end do

        this%total_right_out = total_area_out_tmp
        this%total_right_in  = total_area_in_tmp

        call this%Calc_area_for_each_side(x2, y2, x3, y3, x_lag2, y_lag2, x_lag3, y_lag3, &
            vol_top, vol_center, .false. , 1, i, j        , &
            area_out_tmp, area_in_tmp, total_area_out_tmp , total_area_in_tmp)
        do tmp_mat = 1, this%n_materials
            !            call this%adv_mats%Point_to_areas(area_top_in, area_top_out)
            area_top_out(tmp_mat, i, 1, 1) = area_out_tmp(tmp_mat)
            area_top_in (tmp_mat, i, 1, 1) = area_in_tmp (tmp_mat)
        end do
        this%total_top_out = total_area_out_tmp
        this%total_top_in  = total_area_in_tmp

        do tmp_mat = 1, this%n_materials
            !            call this%adv_mats(tmp_mat)%Point_to_area_top_out(area_top_out)
            this%adv_mats%total_out(tmp_mat) = this%adv_mats%area_right_out(tmp_mat) + area_top_out(tmp_mat, i, 1, 1) + &
                this%adv_mats%area_left_out(tmp_mat)  + this%adv_mats%area_bottom_out(tmp_mat)
        end do
        this%total_out = this%total_right_out + this%total_top_out + this%total_left_out + this%total_bottom_out

        deallocate(area_out_tmp, area_in_tmp)
    end subroutine Calculate_areas_for_advect

    subroutine Calc_area_for_each_side(this, &
        x1    , y1    , x2    , y2    , &
        x1_lag, y1_lag, x2_lag, y2_lag, &
        vol_neighbor, vol,i_not_j, up_not_down, i, j, &
        area_out, area_in, total_area_out, total_area_in)
        use geometry_module , only : Cut_line, Quad_volume, Triangle_volume
        implicit none

        class(advect_t), intent(in out) :: this

        integer, intent(in) :: i, j

        logical, intent(in) :: i_not_j
        integer, intent(in) :: up_not_down

        real(8), dimension(:), allocatable, intent(in out) :: area_out
        real(8), dimension(:), allocatable, intent(in out) :: area_in
        real(8)                           , intent(out) :: total_area_out
        real(8)                           , intent(out) :: total_area_in

        real(8), intent(in) :: vol
        real(8), intent(in) :: vol_neighbor

        real(8), intent(in) :: x1    , y1    , x2    , y2
        real(8), intent(in) :: x1_lag, y1_lag, x2_lag, y2_lag

        real(8) :: x_cut, y_cut
        integer :: cut_inside

        real(8) :: quad_vol_half
        real(8) :: tri_vol_half
        real(8) :: weight

        integer :: tmp_mat

        real(8), dimension(:), allocatable :: area_out_1, area_out_2
        real(8), dimension(:), allocatable :: area_in_1 , area_in_2
        real(8)                            :: total_area_out_1, total_area_out_2
        real(8)                            :: total_area_in_1 , total_area_in_2


        allocate(area_out_1(this%n_materials), area_in_1(this%n_materials), &
            area_out_2(this%n_materials), area_in_2(this%n_materials))
        call Cut_line(x1, y1, x2, y2, x1_lag, y1_lag, x2_lag, y2_lag, cut_inside, x_cut, y_cut, 0)

        if (cut_inside == 0) then

            quad_vol_half = 0.5d0 * Quad_volume(x2_lag, y2_lag, x1_lag, y1_lag, x1, y1, x2, y2, this%mesh%cyl)
            weight = this%a0 * sign(1d0, quad_vol_half) + this%b0 * 4d0 * quad_vol_half / (vol + vol_neighbor)


            call this%Volume_material_old(i_not_j, up_not_down, i, j, quad_vol_half, weight, &
                x1_lag, y1_lag, x2_lag, y2_lag                   , &
                x1    , y1    , x2    , y2                       , &
                area_out, area_in, total_area_out, total_area_in)
        else

            tri_vol_half =  0.5d0 * Triangle_volume(x1, y1, x_cut, y_cut, x1_lag, y1_lag, this%mesh%cyl)
            weight = this%a0 * sign(1d0, tri_vol_half) + this%b0 * 4d0 * tri_vol_half / (vol + vol_neighbor)

            call this%Volume_material_old(i_not_j, up_not_down, i, j, tri_vol_half, weight, &
                x1_lag, y1_lag, x_cut, y_cut                    , &
                x1    , y1    , x_cut, y_cut                    , &
                area_out_1, area_in_1, total_area_out_1, total_area_in_1)

            tri_vol_half = -0.5d0 * Triangle_volume(x2, y2, x_cut, y_cut, x2_lag, y2_lag, this%mesh%cyl)
            weight = this%a0 * sign(1d0, tri_vol_half) + this%b0 * 4d0 * tri_vol_half / (vol + vol_neighbor)
            call this%Volume_material_old(i_not_j, up_not_down, i, j, tri_vol_half, weight, &
                x_cut, y_cut, x2_lag, y2_lag                    , &
                x_cut, y_cut, x2    , y2                        , &
                area_out_2, area_in_2, total_area_out_2, total_area_in_2)

            total_area_in  = total_area_in_1  + total_area_in_2
            total_area_out = total_area_out_1 + total_area_out_2
            do tmp_mat = 1, this%n_materials
                area_in (tmp_mat) = area_in_1 (tmp_mat) + area_in_2 (tmp_mat)
                area_out(tmp_mat) = area_out_1(tmp_mat) + area_out_2(tmp_mat)
            end do
        end if

        deallocate(area_out_1, area_in_1, &
            area_out_2, area_in_2)
        return
    end subroutine Calc_area_for_each_side


    subroutine Calculate_vof_line_new(this, x, y, cyl, ind_mat, i, j, k, nx, ny, nz, a, b, c, &
        x_line1, y_line1, x_line2, y_line2, max_num_iter)
        use geometry_module , only : Quad_volume, Triangle_volume
        class(advect_t), intent(in out) :: this

        real(8)                             , intent(in)  :: cyl
        integer                             , intent(in)  :: ind_mat, i, j, k, nx, ny, nz
        real(8), dimension(:, :, :), pointer, intent(in)  :: x, y
        real(8)                             , intent(out) :: a, b, c, x_line1, y_line1, x_line2, y_line2
        integer                             , intent(in out) :: max_num_iter

        real(8), dimension(:, :, :,:), pointer :: mat_vof
        real(8), dimension(:, :, :), pointer :: vof
        real(8), dimension(:, :, :), pointer :: vol

        real(8), dimension(4) :: c_temp
        real(8), dimension(3) :: cn
        integer, dimension(4) :: ind
        integer, dimension(4) :: i_vert
        integer, dimension(4) :: j_vert


        integer :: im, ip, jm, jp
        integer :: ii, jj, jj_1
        integer :: ind_temp
        integer :: case
        integer :: i_next, i_prev
        integer :: is_next_i_3
        integer :: iter_count

        real(8) :: omcyl
        real(8) :: vol_temp
        real(8) :: f1, f2, f3, f4
        real(8) :: df_dx, df_dy
        real(8) :: cc
        real(8) :: vof_temp
        real(8) :: f
        real(8) :: x1 , x2 , x3 , x4 , y1 , y2 , y3 , y4
        real(8) :: xi1, xi2, xi3, xi4, yi1, yi2, yi3, yi4
        real(8) :: a13, b13, c13, x13, y13, det13,a12, b12, c12, x12, y12, det12, a24, b24, c24, x24, y24, det24, &
            a34, b34, c34,  det34
        real(8) :: fx1, fy1, gx1, gy1, fx2, fy2, gx2, gy2
        real(8) :: vof_exchange
        real(8) :: dv_dc, dv
        real(8) :: vol_new
        real(8) :: c_old
        real(8) :: c_new
        real(8) :: c_low, c_high

        call this%volume%Point_to_data(vol)

        im = max(i - 1, 1 )
        ip = min(i + 1, nx)


        jm = max(j - 1, 1 )
        jp = min(j + 1, ny)

        if (ind_mat /= 0) then
            call this%materials%vof%Point_to_data(mat_vof)
            f1 = 0.25d0 * (mat_vof(ind_mat,ip, j , 1) * vol(ip, j , 1) + mat_vof(ind_mat,ip, jm, 1) * vol(ip, jm, 1) + &
                mat_vof(ind_mat,i , jm, 1) * vol(i , jm, 1) + mat_vof(ind_mat,i , j , 1) * vol(i , j , 1)) / &
                (vol(ip, j , 1) + vol(ip, jm, 1) + vol(i , jm, 1) + vol(i , j , 1))

            f2 = 0.25d0 * (mat_vof(ind_mat,ip, jp, 1) * vol(ip, jp, 1) + mat_vof(ind_mat,ip, j , 1) * vol(ip, j , 1) + &
                mat_vof(ind_mat,i , j , 1) * vol(i , j , 1) + mat_vof(ind_mat,i , jp, 1) * vol(i , jp, 1)) / &
                (vol(ip, jp, 1) + vol(ip, j , 1) + vol(i , j , 1) + vol(i , jp, 1))

            f3 = 0.25d0 * (mat_vof(ind_mat,i , jp, 1) * vol(i , jp, 1) + mat_vof(ind_mat,i , j , 1) * vol(i , j , 1) + &
                mat_vof(ind_mat,im, j , 1) * vol(im, j , 1) + mat_vof(ind_mat,im, jp, 1) * vol(im, jp, 1)) / &
                (vol(i , jp, 1) + vol(i , j , 1) + vol(im, j , 1) + vol(im, jp, 1))

            f4 = 0.25d0 * (mat_vof(ind_mat,i , j , 1) * vol(i , j , 1) + mat_vof(ind_mat,i , jm, 1) * vol(i , jm, 1) + &
                mat_vof(ind_mat,im, jm, 1) * vol(im, jm, 1) + mat_vof(ind_mat,im, j , 1) * vol(im, j , 1)) / &
                (vol(i , j , 1) + vol(i , jm, 1) + vol(im, jm, 1) + vol(im, j , 1))
            vof_temp = mat_vof(ind_mat, i, j, 1)

        else
            call this%total_vof%Point_to_data(vof)
            f1 = 0.25d0 * (vof(ip, j , 1) * vol(ip, j , 1) + vof(ip, jm, 1) * vol(ip, jm, 1) + &
                vof(i , jm, 1) * vol(i , jm, 1) + vof(i , j , 1) * vol(i , j , 1)) / &
                (vol(ip, j , 1) + vol(ip, jm, 1) + vol(i , jm, 1) + vol(i , j , 1))

            f2 = 0.25d0 * (vof(ip, jp, 1) * vol(ip, jp, 1) + vof(ip, j , 1) * vol(ip, j , 1) + &
                vof(i , j , 1) * vol(i , j , 1) + vof(i , jp, 1) * vol(i , jp, 1)) / &
                (vol(ip, jp, 1) + vol(ip, j , 1) + vol(i , j , 1) + vol(i , jp, 1))

            f3 = 0.25d0 * (vof(i , jp, 1) * vol(i , jp, 1) + vof(i , j , 1) * vol(i , j , 1) + &
                vof(im, j , 1) * vol(im, j , 1) + vof(im, jp, 1) * vol(im, jp, 1)) / &
                (vol(i , jp, 1) + vol(i , j , 1) + vol(im, j , 1) + vol(im, jp, 1))

            f4 = 0.25d0 * (vof(i , j , 1) * vol(i , j , 1) + vof(i , jm, 1) * vol(i , jm, 1) + &
                vof(im, jm, 1) * vol(im, jm, 1) + vof(im, j , 1) * vol(im, j , 1)) / &
                (vol(i , j , 1) + vol(i , jm, 1) + vol(im, jm, 1) + vol(im, j , 1))
            vof_temp = vof(i, j, 1)

        end if






        df_dx = (f2 - f4) * (y(i, j + 1, 1) - y(i + 1, j, 1)) - (f1 - f3) * (y(i    , j    , 1) - y(i + 1, j + 1, 1))

        df_dy = (f4 - f2) * (x(i, j + 1, 1) - x(i + 1, j, 1)) - (f1 - f3) * (x(i + 1, j + 1, 1) - x(i    , j    , 1))

        a = df_dx
        b = df_dy


        if (abs(a) < 1d-20 .and. abs(b) < 1d-20) then
            a = 1d0
            b = 0d0
        end if



        c_temp(1) = - a * x(i    , j    , 1) - b * y(i    , j    , 1)
        c_temp(2) = - a * x(i + 1, j    , 1) - b * y(i + 1, j    , 1)
        c_temp(3) = - a * x(i + 1, j + 1, 1) - b * y(i + 1, j + 1, 1)
        c_temp(4) = - a * x(i    , j + 1, 1) - b * y(i    , j + 1, 1)

        i_vert(1) = i
        i_vert(2) = i + 1
        i_vert(3) = i + 1
        i_vert(4) = i

        j_vert(1) = j
        j_vert(2) = j
        j_vert(3) = j + 1
        j_vert(4) = j + 1

        ind(1) = 1
        ind(2) = 2
        ind(3) = 3
        ind(4) = 4



        do ii = 1, 3
            do jj = 1, 4 - ii
                jj_1 = jj + 1
                if(c_temp(jj_1) < c_temp(jj)) then
                    cc           = c_temp(jj)
                    ind_temp     = ind(jj)
                    c_temp(jj)   = c_temp(jj_1)

                    ind(jj)      = ind(jj_1)
                    c_temp(jj_1) = cc
                    ind(jj_1)    = ind_temp
                end if
            end do
        end do

        x1 = x(i_vert(ind(1)), j_vert(ind(1)), 1)
        x2 = x(i_vert(ind(2)), j_vert(ind(2)), 1)

        y1 = y(i_vert(ind(1)), j_vert(ind(1)), 1)
        y2 = y(i_vert(ind(2)), j_vert(ind(2)), 1)

        i_next = ind(1) + 1
        if (i_next == 5) i_next = 1
        i_prev = ind(1) - 1
        if (i_prev == 0) i_prev = 4

        if (ind(3) == i_next .or. ind(3) == i_prev) then
            is_next_i_3 = 1
            x3 = x(i_vert(ind(3)), j_vert(ind(3)), 1)
            x4 = x(i_vert(ind(4)), j_vert(ind(4)), 1)
            y3 = y(i_vert(ind(3)), j_vert(ind(3)), 1)
            y4 = y(i_vert(ind(4)), j_vert(ind(4)), 1)

        else
            is_next_i_3 = 0
            x3 = x(i_vert(ind(4)), j_vert(ind(4)), 1)
            x4 = x(i_vert(ind(3)), j_vert(ind(3)), 1)
            y3 = y(i_vert(ind(4)), j_vert(ind(4)), 1)
            y4 = y(i_vert(ind(3)), j_vert(ind(3)), 1)
        end if

        vol_temp = abs(Quad_volume(x(i    , j    , 1), y(i    , j    , 1), x(i + 1, j    , 1), y(i + 1, j    , 1), &
            x(i + 1, j + 1, 1), y(i + 1, j + 1, 1), x(i    , j + 1, 1), y(i    , j + 1, 1), cyl))
        a13   = y3 - y1
        b13   = x1 - x3
        c13   = (y1 - y3) * x1 + (x3 - x1) * y1
        det13 = a13 * b - a * b13
        x13   = (-c13 * b + c_temp(2) * b13) / det13
        y13   = ( c13 * a - c_temp(2) * a13) / det13

        f = abs(Triangle_volume(x1, y1, x2, y2, x13, y13, cyl)) / abs(vol_temp)
        if (f >= vof_temp) then
            case = 1
            vof_exchange  = 0d0
            a12   = y2 - y1
            b12   = x1 - x2
            c12   = (y1 - y2) * x1 + (x2 - x1) * y1
            det12 = a12 * b - a * b12

            fx1 =  b13     / det13
            gx1 = -c13 * b / det13
            fy1 = -a13     / det13
            gy1 =  c13 * a / det13
            fx2 =  b12     / (det12 + 1d-30)
            gx2 = -c12 * b / (det12 + 1d-30)
            fy2 = -a12     / (det12 + 1d-30)
            gy2 =  c12 * a / (det12 + 1d-30)
            xi1 =  x1
            yi1 =  y1
            xi2 =  x1
            yi2 =  y1

            if (Quad_volume(xi1, yi1, xi2, yi2, x3, y3, x2, y2, cyl) < 0d0) then
                fx2 =  b13     / det13
                gx2 = -c13 * b / det13
                fy2 = -a13     / det13
                gy2 =  c13 * a / det13
                fx1 =  b12     / (det12 + 1d-30)
                gx1 = -c12 * b / (det12 + 1d-30)
                fy1 = -a12     / (det12 + 1d-30)
                gy1 =  c12 * a / (det12 + 1d-30)
            end if

        else

            a24   = y4 - y2
            b24   = x2 - x4
            c24   = (y2 - y4) * x2 + (x4 - x2) * y2
            det24 = a24 * b - a * b24
            if (is_next_i_3 == 1) then
                x24  = (-c24 * b + c_temp(3) * b24) / det24
                y24  = ( c24 * a - c_temp(3) * a24) / det24
                vof_exchange = f
                f   = f + abs(Quad_volume(x2, y2, x24, y24, x3, y3, x13, y13, cyl)) / abs(vol_temp)
            else
                x24  = (-c13 * b + c_temp(3) * b13) / det13
                y24  = ( c13 * a - c_temp(3) * a13) / det13
                vof_exchange = f
                f   = f + abs(Quad_volume(x2, y2, x4, y4, x24, y24, x13, y13, cyl)) / abs(vol_temp)
            end if

            if (f >= vof_temp) then
                case = 2
                fx1 =  b13     / det13
                gx1 = -c13 * b / det13
                fy1 = -a13     / det13
                gy1 =  c13 * a / det13
                fx2 =  b24     / det24
                gx2 = -c24 * b / det24
                fy2 = -a24     / det24
                gy2 =  c24 * a / det24
                xi1 =  x2
                yi1 =  y2
                xi2 =  x13
                yi2 =  y13
                if ((is_next_i_3 == 1 .and. Quad_volume(xi1, yi1, xi2, yi2, x3 ,  y3, x24, y24, cyl) < 0d0) .or.  &
                    (is_next_i_3 == 0 .and. Quad_volume(xi1, yi1, xi2, yi2, x24, y24,  x4,  y4, cyl) < 0d0)) then
                    fx2 =  b13     / det13
                    gx2 = -c13 * b / det13
                    fy2 = -a13     / det13
                    gy2 =  c13 * a / det13
                    fx1 =  b24     / det24
                    gx1 = -c24 * b / det24
                    fy1 = -a24     / det24
                    gy1 =  c24 * a / det24
                    xi2 =  x2
                    yi2 =  y2
                    xi1 =  x13
                    yi1 =  y13
                end if

            else
                case = 3
                vof_exchange  = f
                a34   = y4 - y3
                b34   = x3 - x4
                c34   = (y3 - y4) * x3 + (x4 - x3) * y3
                det34 = a34 * b - a * b34
                if (is_next_i_3 == 1) then
                    fx1 =  b24     / det24
                    gx1 = -c24 * b / det24
                    fy1 = -a24     / det24
                    gy1 =  c24 * a / det24
                else
                    fx1 =  b13     / det13
                    gx1 = -c13 * b / det13
                    fy1 = -a13     / det13
                    gy1 =  c13 * a / det13
                end if
                fx2 =  b34     / det34
                gx2 = -c34 * b / det34
                fy2 = -a34     / det34
                gy2 =  c34 * a / det34
                if (is_next_i_3 == 1) then
                    xi1 = x3
                    yi1 = y3
                else
                    xi1 = x4
                    yi1 = y4
                end if
                xi2 = x24
                yi2 = y24
                if ((is_next_i_3 == 1 .and. Quad_volume(xi1, yi1, xi2, yi2, x4, y4, x4, y4, cyl) < 0d0) .or.  &
                    (is_next_i_3 == 0 .and. Quad_volume(xi1, yi1, xi2, yi2, x3, y3, x3, y3, cyl) < 0d0)) then
                    if (is_next_i_3 == 1) then
                        fx2 =  b24     / det24
                        gx2 = -c24 * b / det24
                        fy2 = -a24     / det24
                        gy2 =  c24 * a / det24
                    else
                        fx2 =  b13     / det13
                        gx2 = -c13 * b / det13
                        fy2 = -a13     / det13
                        gy2 =  c13 * a / det13
                    end if
                    fx1 =  b34     / det34
                    gx1 = -c34 * b / det34
                    fy1 = -a34     / det34
                    gy1 =  c34 * a / det34
                    if (is_next_i_3 == 1) then
                        xi2 = x3
                        yi2 = y3
                    else
                        xi2 = x4
                        yi2 = y4
                    end if
                    xi1 = x24
                    yi1 = y24
                end if
            end if
        end if

        dv     = (vof_temp - vof_exchange) * vol_temp
        c_low  = c_temp(case)
        c_high = c_temp(case + 1)
        c_new = (c_low + c_high) / 2d0
        iter_count = 0


        omcyl = 1d0 - cyl
        c_old = c_new
        vol_new = 100

        do while (iter_count == 0 .or. ((abs((c_old - c_new) / (c_old + c_new + 1d-38)) > 1d-6 .or. &
            abs((vol_new - dv) / (vol_new + dv + 1d-38)) > 1d-3) .and. iter_count < 99))
            c_old = c_new
            iter_count = iter_count + 1
            xi3 = gx1 + fx1 * c_new
            xi4 = gx2 + fx2 * c_new
            yi3 = gy1 + fy1 * c_new
            yi4 = gy2 + fy2 * c_new
            dv_dc = (fx1 * cyl * ((xi3 - xi2) * (yi1 - yi2) - (xi1 - xi2) * (yi3 - yi2)) + &
                ((xi1 + xi2 + xi3) * cyl + 3d0 * omcyl) * &
                (fx1 * (yi1 - yi2) - (xi1 - xi2) * fy1) + &
                (fx1 + fx2) * cyl * ((xi1 - xi4) * (yi3 - yi4) - (xi3 - xi4) * (yi1 - yi4)) + &
                ((xi3 + xi4 + xi1) * cyl + 3d0 * omcyl) * &
                (-fx2 * (yi3 - yi4) + (xi1 - xi4) * (fy1 - fy2) - &
                (fx1 - fx2) * (yi1 - yi4) - (xi3 - xi4) * (-fy2))) / 6d0
            vol_new = Quad_volume(xi1, yi1, xi2, yi2, xi3, yi3, xi4, yi4, cyl)
            if (dv_dc /= 0d0) then
                c_new = c_old - (vol_new - dv) / dv_dc
            else
            end if
            if (c_new > c_high) c_new = 2d0 * c_high - c_new
            if (c_new < c_low)  c_new = 2d0 * c_low  - c_new
        end do



        max_num_iter = MAX(max_num_iter, iter_count)
        c = c_new
        x_line1 = gx1 + fx1 * c
        x_line2 = gx2 + fx2 * c
        y_line1 = gy1 + fy1 * c
        y_line2 = gy2 + fy2 * c

        return
    end subroutine Calculate_vof_line_new


    subroutine Volume_material_new(this, i_not_j, up_not_down, i, j, f, a1, xp1, yp1, xp2, yp2, &
        x1, y1, x2, y2, area_out, area_in, total_out, total_in)
        use geometry_module, only : Volume_fraction2
        implicit none

        class(advect_t), intent(in out) :: this

        integer              , intent(in)    :: i, j, up_not_down
        logical              , intent(in)    :: i_not_j
        real(8)              , intent(in)    :: xp1, yp1, xp2, yp2, x1, y1, x2, y2, f, a1
        real(8)              , intent(out)   :: total_out, total_in
        real(8), dimension(:), intent(in out) :: area_out
        real(8), dimension(:), intent(in out) :: area_in

        real(8), dimension(:, :, :), pointer :: x_coor
        real(8), dimension(:, :, :), pointer :: y_coor
        real(8), dimension(:, :, :), pointer :: vof
        real(8), dimension(:, :, :, :), pointer :: mat_vof
        real(8), dimension(:, :, :), pointer :: vol_lag
        real(8), dimension(:, :, :), pointer :: mat_id
        real(8), dimension(:, :, :), pointer :: n_materials_in_cell

        real(8)               :: f_tot, fph, sg, f_tot_with_rem, df, df1, norm_f, &
            a, b, c, side, xx1, yy1, xx2 ,yy2, vof_max, mat_in_cell
        real(8)               :: emf1
        integer               :: tmp_mat, ii_max, numb, id, ia, ja, jd, n_mats, mat_index
        integer :: max_num_iter
        real(8), dimension(:), allocatable :: fm, fre_main


        allocate(fm(this%n_materials), fre_main(this%n_materials))
        emf1 = 1 - this%emf
        n_mats = this%n_materials

        call this%total_vof%Point_to_data(vof)
        call this%rezone   %Point_to_volume(vol_lag)
        call this%rezone   %Point_to_coordinates_2d(x_coor, y_coor)
        call this%mat_cells%Point_to_data(mat_id)
        call this%num_mat_cells%Point_to_data(n_materials_in_cell)
        call this%Find_donnor_acceptor(f, i_not_j, up_not_down, i, j, ia, ja, id, jd)
        call this%materials%vof%Point_to_data(mat_vof)

        area_out = 0d0
        area_in = 0d0
        fm = 0d0
        fre_main = 0d0

        if (vof(id, jd, 1) < this%emf) then
            fph = 0.d0
        else
            if (vof(id, jd, 1) < emf1) then

                call this%Calculate_vof_line_new(x_coor, y_coor, this%mesh%cyl, 0, &
                    id, jd, 1, this%nx, this%ny, this%nz, &
                    a, b, c, xx1, yy1, xx2, yy2, max_num_iter)

                fph = Volume_fraction2(a, b, c, f, xp2, yp2, xp1, yp1, x1, y1, x2, y2, this%mesh%cyl)


            else
                fph = f
            end if
            f_tot = 0d0

            if (n_materials_in_cell(id, jd, 1) > 1) then
                do tmp_mat = 1, n_mats
                    if (mat_vof(tmp_mat, id, jd, 1) < this%emf) cycle
                    call this%Calculate_vof_line_new(x_coor, y_coor, this%mesh%cyl, tmp_mat, &
                        id, jd, 1, this%nx, this%ny, this%nz,   &
                        a, b, c, xx1, yy1, xx2, yy2, max_num_iter)
                    fm(tmp_mat) = Volume_fraction2(a, b, c, f, xp2, yp2, xp1, yp1, x1, y1, x2, y2, this%mesh%cyl)
                    sg = sign(1.0d0, fm(tmp_mat))
                    mat_in_cell = 0.5d0 * mat_vof(tmp_mat,id, jd, 1) * vol_lag(id, jd, 1)
                    fm(tmp_mat) = sg * min(dabs(fm(tmp_mat)), mat_in_cell)
                    fre_main(tmp_mat) = mat_in_cell - dabs(fm(tmp_mat))

                    if (fre_main(tmp_mat) > 0d0) then
                        f_tot_with_rem = f_tot_with_rem + fm(tmp_mat)
                    end if
                    f_tot = f_tot + fm(tmp_mat)
                end do
                df = fph - f_tot
                sg = sign(1d0, fph)
                if (sg * df > 0) then
                    do while(df /= 0d0)
                        norm_f = 1 / f_tot_with_rem
                        f_tot_with_rem = 0d0
                        df1 = dabs(df)
                        do tmp_mat = 1, n_mats
                            if (fre_main(tmp_mat) == 0d0) cycle
                            df1 = min(df1, fre_main(tmp_mat) / (fm(tmp_mat) * norm_f))
                        end do

                        df1 = df1 * sign(1d0, df)
                        do tmp_mat = 1, n_mats
                            if (fre_main(tmp_mat) == 0d0) cycle

                            fm(tmp_mat) = fm(tmp_mat) + df1 * fm(tmp_mat) * norm_f
                            fre_main(tmp_mat) = fre_main(tmp_mat) - dabs(df1 * fm(tmp_mat) * norm_f)
                            if (fre_main(tmp_mat) > 0d0) then
                                f_tot_with_rem = f_tot_with_rem + fm(tmp_mat)
                            end if
                        end do
                        df = df - df1

                    end do
                else if (sg * df < 0) then
                    norm_f = 1 / f_tot
                    do tmp_mat = 1, n_mats
                        if (mat_vof(tmp_mat,id, jd, 1) < this%emf) cycle
                        fm(tmp_mat) = fm(tmp_mat) + df * fm(tmp_mat) * norm_f
                    end do
                end if
            else if (mat_id(id, jd, 1) > 0) then
                fm(int(mat_id(id, jd, 1))) = fph
            else
                do tmp_mat = 1, n_mats
                    if (mat_vof(tmp_mat,id, jd, 1) >= this%emf) then
                        fm(tmp_mat) = fph
                        exit
                    end if
                end do
            end if
        end if
        if ((vof(ia, ja, 1) > emf1) .and. (vof(id, jd, 1) > emf1) .and.  &
            (mat_id(ia, ja, 1) == mat_id(id, jd, 1)) .and. (n_materials_in_cell(id, jd, 1) == 1)) then
            mat_index = int(mat_id(id, jd, 1))
            area_out(mat_index) = f * (1d0 - a1)
            area_in (mat_index) = f * (1d0 + a1)
            total_out = area_out(mat_index)
            total_in  = area_in (mat_index)
        else
            total_out = 0d0
            total_in  = 0d0
            do tmp_mat = 1, n_mats
                sg = sign(1d0, fm(tmp_mat))
                area_out(tmp_mat) = fm(tmp_mat) * (1d0 - sg)
                area_in (tmp_mat) = fm(tmp_mat) * (1d0 + sg)
                total_out = total_out + area_out(tmp_mat)
                total_in  = total_in  + area_in(tmp_mat)
            end do
        end if
        deallocate(fm, fre_main)
        return
    end subroutine Volume_material_new


    subroutine Volume_material_old(this, i_not_j, up_not_down, i, j, f ,a1 ,xp1, yp1, xp2, yp2, &
        x1, y1, x2, y2, area_out, area_in, total_out, total_in)
        use geometry_module, only : Volume_fraction2
        implicit none

        class(advect_t), intent(inout) :: this

        integer              , intent(in)    :: i, j, up_not_down
        logical              , intent(in)    :: i_not_j
        real(8)              , intent(in)    :: xp1, yp1, xp2, yp2, x1, y1, x2, y2, f, a1
        real(8)              , intent(out)   :: total_out, total_in
        real(8), dimension(:), intent(in out) :: area_out
        real(8), dimension(:), intent(in out) :: area_in

        real(8), dimension(:, :, :), pointer :: x_coor
        real(8), dimension(:, :, :), pointer :: y_coor
        real(8), dimension(:, :, :), pointer :: vof
        real(8), dimension(:, :, :), pointer :: mat_id
        real(8), dimension(:, :, :), pointer :: n_materials_in_cell

        real(8), dimension(:, :, :, :), pointer :: mat_vof
        real(8), dimension(:, :, :, :), pointer :: density_vof
        real(8), dimension(:, :, :, :), pointer :: cell_mass_vof

        real(8) :: f_tot, fph, sg, fre_main_tot, df, norm_f, emf1
        real(8) :: a, b, c, side, xx1, yy1, xx2 ,yy2, vof_max, mat_in_cell, fmmt
        integer :: tmp_mat, indx_max, num_mat_cell_tmp, id, ia, ja, jd, n_mats, ndi, ndj, mat_index
        integer :: max_num_iter
        real(8), dimension(:), allocatable :: fm

        allocate(fm(this%n_materials))
        emf1 = 1 - this%emf
        n_mats = this%n_materials

        call this%rezone   %Point_to_coordinates_2d(x_coor, y_coor)
        call this%total_vof%Point_to_data(vof)
        call this%mat_cells%Point_to_data(mat_id)
        call this%num_mat_cells%Point_to_data(n_materials_in_cell)
        call this%materials%vof%Point_to_data(mat_vof)
        call this%materials%cell_mass%Point_to_data(cell_mass_vof)
        call this%materials%density  %Point_to_data(density_vof)
        call this%Find_donnor_acceptor(f, i_not_j, up_not_down, i, j, ia, ja, id, jd)

        area_out = 0d0
        area_in = 0d0
        fm = 0d0

        if (vof(id, jd, 1) < this%emf) then
            fph = 0d0
        else
            if (vof(id, jd, 1) < emf1) then
                call this%Calculate_vof_line_new(x_coor, y_coor, this%mesh%cyl, 0, &
                    id, jd, 1, this%nx, this%ny, this%nz, &
                    a, b, c, xx1, yy1, xx2, yy2, max_num_iter)
                fph = Volume_fraction2(a, b, c, f, xp2, yp2, xp1, yp1, x1, y1, x2, y2, this%mesh%cyl)
            else
                fph = f
            end if

            f_tot = 0d0
            num_mat_cell_tmp = 1
            vof_max = 0d0
            indx_max = n_mats

            do tmp_mat = 1, n_mats
                if (mat_vof(tmp_mat,id, jd, 1) < this%emf) cycle

                if (num_mat_cell_tmp < n_materials_in_cell(id, jd, 1)) then

                    num_mat_cell_tmp = num_mat_cell_tmp + 1



                    call this%Calculate_vof_line_new (x_coor, y_coor, this%mesh%cyl, tmp_mat, &
                        id, jd, 1, this%nx, this%ny, this%nz,   &
                        a, b, c, xx1, yy1, xx2, yy2, max_num_iter)
                    fm(tmp_mat) = Volume_fraction2 (a, b, c, f, xp2, yp2, xp1, yp1, x1, y1, x2, y2, this%mesh%cyl)
                    if (abs(fm(tmp_mat)) > vof_max) then
                        vof_max  = abs(fm(tmp_mat))
                        indx_max = tmp_mat
                    end if
                    sg = sign(1d0, fm(tmp_mat))
                    if (sg == sign(1d0, f)) then
                        ndi = id
                        ndj = jd
                    else
                        ndi = ia
                        ndj = ja
                    end if
                    fm(tmp_mat) = sg * min(abs(fm(tmp_mat)), 0.5d0 * cell_mass_vof(tmp_mat,ndi, ndj, 1) / (density_vof(tmp_mat,ndi, ndj, 1) + 2d-38))
                    f_tot  = f_tot + fm(tmp_mat)
                else
                    fm(tmp_mat) = fph - f_tot
                    sg = sign(1d0, fm(tmp_mat))
                    if (sg == sign(1d0, f)) then
                        ndi = id
                        ndj = jd
                    else
                        ndi = ia
                        ndj = ja
                    end if
                    fmmt   = sg * min(abs(fm(tmp_mat)), 0.5d0 * cell_mass_vof(tmp_mat,ndi, ndj, 1) / (density_vof(tmp_mat,ndi, ndj, 1) + 2d-38))
                    fm(indx_max) = fm(indx_max) + fm(tmp_mat) - fmmt
                    fm(tmp_mat) = fmmt
                end if
            end do
        end if

        if ((vof(ia, ja, 1) > emf1) .and. (vof(id, jd, 1) > emf1) .and.  &
            (mat_id(ia, ja, 1) == mat_id(id, jd, 1)) .and. (n_materials_in_cell(id, jd, 1) == 1)) then
            mat_index = int(mat_id(id, jd, 1))
            area_out(mat_index) = f * (1d0 - a1)
            area_in (mat_index) = f * (1d0 + a1)
            total_out = area_out(mat_index)
            total_in  = area_in (mat_index)
        else
            total_out = 0d0
            total_in  = 0d0
            do tmp_mat = 1, n_mats
                sg = sign(1d0, fm(tmp_mat))
                area_out(tmp_mat) = fm(tmp_mat) * (1d0 - sg)
                area_in (tmp_mat) = fm(tmp_mat) * (1d0 + sg)
                total_out = total_out + area_out(tmp_mat)
                total_in  = total_in  + area_in(tmp_mat)
            end do
        end if
        deallocate(fm)
        return
    end subroutine Volume_material_old


    subroutine Volume_material_3d(this, i, j, k, itd, jtd, k_donner, ita, jta, k_accept, fxt, axt, &
        xl1, yl1, zl1, xl2, yl2, zl2, xl3, yl3, zl3, xl4, yl4, zl4, x1, y1, &
        z1, x2, y2, z2, x3, y3, z3, x4, y4, z4, vol1, vol2, vol3, vol4, iprint)
        use geometry_module, only : Volume_fraction_3d
        implicit none

        class(advect_t), intent(inout)   :: this
        integer, intent(in)              :: i, j, k, itd, jtd, k_donner, ita, jta, k_accept, iprint
        real(8), intent(in)              :: fxt, axt, xl1, yl1, zl1, xl2, yl2, zl2, xl3, yl3, &
            zl3, xl4, yl4, zl4, x1, y1, z1, x2, y2, z2, x3, y3, &
            z3, x4, y4, z4, vol1, vol2, vol3, vol4
        real(8)                          :: fxtmt, xc1, yc1, zc1, xc2, yc2, zc2, xlc1, ylc1, zlc1, xlc2, &
            ylc2, zlc2, ftph, vofmax, ftot, sg, a, b, c, side

        real(8), dimension(:, :, :), pointer :: tot_vof, num_mat_cell

        real(8), dimension(:, :, :), pointer :: vof
        real(8), dimension(:, :, :), pointer :: density
        real(8), dimension(:, :, :), pointer :: cell_mass

        real(8), dimension(:, :, :), pointer :: mat_id
        real(8), dimension(:, :, :, :), pointer :: mat_vof
        real(8), dimension(:, :, :, :), pointer :: density_vof
        real(8), dimension(:, :, :, :), pointer :: a_adv
        real(8), dimension(:, :, :, :), pointer :: b_adv
        real(8), dimension(:, :, :, :), pointer :: c_adv
        real(8), dimension(:, :, :, :), pointer :: side_adv
        real(8), dimension(:, :, :, :), pointer :: a_mat, b_mat, c_mat, side_mat
        real(8), dimension(:, :, :, :), pointer :: cell_mass_vof

        integer :: mat_max
        integer                              :: ii, numb,  id, jd, kd



        call this%total_vof      %Point_to_data(tot_vof)
        call this%mat_cells      %Point_to_data(mat_id)
        call this%num_mat_cells  %Point_to_data(num_mat_cell)
        call this%total_density  %Point_to_data(density)
        call this%total_cell_mass%Point_to_data(cell_mass)
        call this%materials%vof%Point_to_data(mat_vof)
        call this%materials%density%Point_to_data(density_vof)
        call this%materials%cell_mass%Point_to_data(cell_mass_vof)
        call this%adv_mats%a%Point_to_data(a_mat)
        call this%adv_mats%b%Point_to_data(b_mat)
        call this%adv_mats%c%Point_to_data(c_mat)
        call this%adv_mats%side%Point_to_data(side_mat)
        call this%a%Point_to_data(a_adv)
        call this%b%Point_to_data(b_adv)
        call this%c%Point_to_data(c_adv)
        call this%side%Point_to_data(side_adv)

        do ii = 1, this%n_materials
            this%adv_mats%sxt1(ii) = 0d0
            this%adv_mats%sxt2(ii) = 0d0
            this%adv_mats%fxtm(ii) = 0d0
        end do
        xc1 = (x1 + x3) / 2d0
        yc1 = (y1 + y3) / 2d0
        zc1 = (z1 + z3) / 2d0
        xc2 = (x2 + x4) / 2d0
        yc2 = (y2 + y4) / 2d0
        zc2 = (z2 + z4) / 2d0

        xlc1 = (xl1 + xl3) / 2d0
        ylc1 = (yl1 + yl3) / 2d0
        zlc1 = (zl1 + zl3) / 2d0
        xlc2 = (xl2 + xl4) / 2d0
        ylc2 = (yl2 + yl4) / 2d0
        zlc2 = (zl2 + zl4) / 2d0
        if (tot_vof(itd, jtd, k_donner) < this%emf) return

        if (tot_vof(itd, jtd, k_donner) < this%emf1) then
            a = a_adv(1, itd, jtd, k_donner)
            b = b_adv(1, itd, jtd, k_donner)
            c = c_adv(1, itd, jtd, k_donner)
            side = side_adv(1, itd, jtd, k_donner)

            ftph = (Volume_fraction_3d (xl1 , yl1 , zl1 , xl2 , yl2 , zl2 , xl3 , yl3 , zl3 , xlc1, ylc1, zlc1,                      &
                x1  , y1  , z1  , x2  , y2  , z2  , x3  , y3  , z3  , xc1 , yc1 , zc1 , a, b, c, side, vol1) &
                + Volume_fraction_3d (xl1 , yl1 , zl1 , xlc1, ylc1, zlc1, xl3 , yl3 , zl3 , xl4 , yl4 , zl4 ,                      &
                x1  , y1  , z1  , xc1 , yc1 , zc1 , x3  , y3  , z3  , x4  , y4  , z4  , a, b, c, side, vol2) &
                + Volume_fraction_3d (xl1 , yl1 , zl1 , xl2 , yl2 , zl2 , xlc2, ylc2, zlc2, xl4 , yl4 , zl4 ,                      &
                x1  , y1  , z1  , x2  , y2  , z2  , xc2 , yc2 , zc2 , x4  , y4  , z4  , a, b, c, side, vol3) &
                + Volume_fraction_3d (xlc2, ylc2, zlc2, xl2 , yl2 , zl2 , xl3 , yl3 , zl3 , xl4 , yl4 , zl4 ,                      &
                xc2 , yc2 , zc2 , x2  , y2  , z2  , x3  , y3  , z3  , x4  , y4  , z4  , a, b, c, side, vol4))&
                * 0.25d0
        else
            ftph = fxt
        end if

        numb = 1
        ftot = 0d0
        vofmax = 0d0
        mat_max = this%n_materials

        do ii = 1, this%n_materials

            if (mat_vof(ii, itd, jtd, k_donner) < this%emf) cycle

            if (numb < num_mat_cell(itd, jtd, k_donner)) then
                numb = numb + 1
                if (mat_vof(ii, itd, jtd, k_donner) > vofmax) then
                    vofmax = mat_vof (ii, itd, jtd, k_donner)
                    mat_max = ii
                end if

                this%adv_mats%fxtm(ii) = 0.25d0*( &
                    Volume_fraction_3d (xl1, yl1, zl1, xl2, yl2, zl2, xl3, yl3, zl3, &
                    xlc1, ylc1, zlc1, x1, y1, z1, x2, y2, z2, x3, y3, z3, xc1, yc1, zc1, &
                    a, b, c, side, vol1)&
                    + &
                    Volume_fraction_3d (xl1, yl1, zl1, xlc1, ylc1, zlc1, xl3, yl3, zl3, &
                    xl4, yl4, zl4, x1, y1, z1, xc1, yc1, zc1, x3, y3, z3, x4, y4, z4, &
                    a, b, c, side, vol2)&
                    + &
                    Volume_fraction_3d (xl1, yl1, zl1, xl2, yl2, zl2, xlc2, ylc2, zlc2, &
                    xl4, yl4, zl4, x1, y1, z1, x2, y2, z2, xc2, yc2, zc2, x4, y4, z4, &
                    a, b, c, side, vol3)&
                    + &
                    Volume_fraction_3d (xlc2, ylc2, zlc2, xl2, yl2, zl2, xl3, yl3, zl3, &
                    xl4, yl4, zl4, xc2, yc2, zc2, x2, y2, z2, x3, y3, z3, x4, y4, z4, &
                    a, b, c, side, vol4))

                sg = sign(1.d0, this%adv_mats%fxtm(ii))
                if (sg  == sign(1.d0, fxt)) then
                    id = itd
                    jd = jtd
                    kd = k_donner
                else
                    id = ita
                    jd = jta
                    kd = k_accept
                end if



                this%adv_mats%fxtm(ii) = sg * min(abs(this%adv_mats%fxtm(ii)), &
                    0.5 * cell_mass_vof(ii, id,jd,kd) / (density_vof(ii, id,jd,kd) + 1d-20))


                this%adv_mats%fxtm(ii) = sg*min(abs(ftph - ftot), abs(this%adv_mats%fxtm(ii)))
                ftot = ftot + this%adv_mats%fxtm(ii)
            else

                this%adv_mats%fxtm(ii) = ftph - ftot
                sg = sign(1.d0,this%adv_mats%fxtm(ii))
                if (sg == sign(1.d0, fxt)) then
                    id = itd
                    jd = jtd
                    kd = k_donner
                else
                    id = ita
                    jd = jta
                    kd = k_accept
                end if

                fxtmt = sg * min(abs(this%adv_mats%fxtm(ii)), 0.5*cell_mass_vof(ii, id, jd, kd) / (density_vof(ii, id, jd, kd) + 1d-20))

                if (mat_max > 0) then
                    this%adv_mats%fxtm(mat_max) = this%adv_mats%fxtm(mat_max) + this%adv_mats%fxtm(ii) - fxtmt

                end if
                this%adv_mats%fxtm(ii) = fxtmt

            end if
        end do

        if ((tot_vof(ita, jta, k_accept) > this%emf1).and. (tot_vof(itd, jtd, k_donner) > this%emf1).and. &
            (mat_id(ita, jta, k_accept) == mat_id(itd, jtd, k_donner)).and. (num_mat_cell(itd, jtd, k_donner) == 1)) then
            ii = mat_id(itd, jtd, k_donner)
            this%adv_mats%sxt1(ii) = this%adv_mats%fxtm(ii) * (1.d0 - axt)
            this%adv_mats%sxt2(ii) = this%adv_mats%fxtm(ii) * (1.d0 + axt)
        else
            do ii = 1, this%n_materials
                this%adv_mats%sgn = sign(1d0, this%adv_mats%fxtm(ii))
                this%adv_mats%sxt1(ii) = this%adv_mats%fxtm(ii) * (1d0 - this%adv_mats%sgn(ii))
                this%adv_mats%sxt2(ii) = this%adv_mats%fxtm(ii) * (1d0 + this%adv_mats%sgn(ii))
            end do
        end if

    end subroutine Volume_material_3d

    subroutine Find_donnor_acceptor(this, f, i_not_j, up_not_down, i, j, ia, ja, id, jd)
        implicit none

        class(advect_t), intent(in out)    :: this

        real(8)        , intent(in)       :: f
        logical        , intent(in)       :: i_not_j
        integer        , intent(in)       :: i, j, up_not_down
        integer        , intent(out)      :: ia, ja, id, jd

        if (i_not_j) then
            jd = j
            ja = j
            if (f <= 0d0) then
                id = i
                ia = i + up_not_down
            else
                ia = i
                id = i + up_not_down
            end if
        else
            id = i
            ia = i
            if (f <= 0d0) then
                jd = j
                ja = j + up_not_down
            else
                ja = j
                jd = j + up_not_down
            end if
        end if
    end subroutine Find_donnor_acceptor

    subroutine Put_quantities_into_cell(this, i, j, k)
        implicit none

        class(advect_t), intent(in out) :: this

        integer, intent(in) :: i, j, k

        real(8), dimension(:, :, :), pointer :: cell_mass
        real(8), dimension(:, :, :), pointer :: cell_mass_adv
        real(8), dimension(:, :, :, :), pointer :: cell_mass_vof
        real(8), dimension(:, :, :, :), pointer :: mat_cell_mass_adv
        real(8), dimension(:, :, :), pointer :: vof
        real(8), dimension(:, :, :), pointer :: vof_adv
        real(8), dimension(:, :, :, :), pointer :: mat_vof
        real(8), dimension(:, :, :, :), pointer :: mat_vof_adv
        real(8), dimension(:, :, :), pointer :: n_materials_in_cell
        real(8), dimension(:, :, :), pointer :: mat_id
        real(8), dimension(:, :, :), pointer :: sie
        real(8), dimension(:, :, :, :), pointer :: sie_vof
        real(8), dimension(:, :, :, :), pointer :: sie_vof_adv
        real(8), dimension(:, :, :, :), pointer :: init_mat_layers
        real(8), dimension(:, :, :, :), pointer :: init_mat_layers_adv

        integer :: tmp_mat
        integer :: index_factor
        real(8) :: vof_correction
        real(8) :: vof_factor

        real(8), dimension(:, :, :), pointer :: mat_vof1

        call this%total_vof %Point_to_data(vof)
        call this%vof_advect%Point_to_data(vof_adv)
        call this%mat_cells %Point_to_data(mat_id)
        call this%total_sie %Point_to_data(sie)
        call this%total_cell_mass %Point_to_data(cell_mass)
        call this%adv_cell_mass%Point_to_data(cell_mass_adv)
        call this%num_mat_cells   %Point_to_data(n_materials_in_cell)
        call this%materials%cell_mass%Point_to_data(cell_mass_vof)
        call this%materials%sie      %Point_to_data(sie_vof)
        call this%materials%vof      %Point_to_data(mat_vof)
        call this%materials%Point_to_initial_layers(init_mat_layers)

        call this%adv_mats%vof       %Point_to_data(mat_vof_adv)
        call this%adv_mats%cell_mass %Point_to_data(mat_cell_mass_adv)
        call this%adv_mats%sie%Point_to_data(sie_vof_adv)
        call this%adv_mats%Point_to_initial_layers(init_mat_layers_adv)
        cell_mass(i, j, 1) = cell_mass_adv(i, j, 1)
        vof(i, j, 1)       = vof_adv(i, j, 1)




        if ((vof(i, j, 1) < this%emf) .or. (cell_mass(i, j, 1) < this%emfm)) then
            cell_mass(i, j, 1) = 0d0
            vof      (i, j, 1) = 0d0
            n_materials_in_cell(i, j, 1) = 0
            mat_id             (i, j, 1) = 0


            sie(i, j, 1) = 0d0
            do tmp_mat = 1, this%n_materials

                cell_mass_vof(tmp_mat, i, j, 1) = 0d0
                sie_vof(tmp_mat, i, j, 1) = 0d0
                mat_vof(tmp_mat, i, j, 1) = 0d0



                if (this%advect_init_layer_mat > 0) init_mat_layers(tmp_mat, i, j, 1) = 0d0


            end do


        else





            sie      (i, j, 1) = 0d0
            cell_mass(i, j, 1) = 0d0
            vof      (i, j, 1) = 0d0
            mat_id   (i, j, 1) = 0
            index_factor = 1
            n_materials_in_cell(i, j, 1) = 0
            vof_correction = 0d0
            do tmp_mat = 1, this%n_materials


                if ((mat_vof_adv(tmp_mat, i, j, 1) > this%emf) .and. (mat_cell_mass_adv(tmp_mat,i, j, 1) > this%emfm)) then

                    mat_vof(tmp_mat,i, j, 1) = mat_vof_adv(tmp_mat,i, j, 1)
                    cell_mass_vof(tmp_mat,i, j, 1) = mat_cell_mass_adv(tmp_mat,i, j, 1)
                    sie_vof(tmp_mat,i, j, 1) = sie_vof_adv(tmp_mat,i, j, 1)
                    if (this%advect_init_layer_mat > 0) init_mat_layers(tmp_mat,i, j, 1) = init_mat_layers_adv(tmp_mat,i, j, 1)
                    sie(i, j, 1)                 = sie(i, j, 1)       + sie_vof(tmp_mat,i, j, 1) * cell_mass_vof(tmp_mat,i, j, 1)
                    cell_mass(i, j, 1)           = cell_mass(i, j, 1) + cell_mass_vof(tmp_mat,i, j, 1)
                    vof(i, j, 1)                 = vof(i, j, 1)       + mat_vof(tmp_mat,i, j, 1)
                    n_materials_in_cell(i, j, 1) = n_materials_in_cell(i, j, 1) + 1
                    mat_id(i, j, 1)              = mat_id(i, j, 1)    + tmp_mat * index_factor
                    index_factor                 = index_factor * max(this%n_materials, 10)
                else
                    vof_correction         = vof_correction + mat_vof_adv(tmp_mat,i, j, 1)
                    mat_vof(tmp_mat, i, j, 1)       = 0.d0
                    cell_mass_vof(tmp_mat, i, j, 1) = 0.d0
                    sie_vof(tmp_mat, i, j, 1)       = 0.d0

                    if (this%advect_init_layer_mat > 0) init_mat_layers(tmp_mat,i, j, 1) = 0d0

                end if

                if ((abs(mat_vof_adv(tmp_mat,i, j, 1) - mat_vof(tmp_mat,i, j, 1)) > this%emf * 1d1) .and. (cell_mass_vof(tmp_mat,i, j, 1) > this%emfm)) then
                    mat_vof(tmp_mat,i, j, 1) = mat_vof_adv(tmp_mat,i, j, 1)
                end if
            end do
            if (cell_mass(i, j, 1) > this%emfm) sie(i, j, 1) = sie(i, j, 1) / cell_mass(i, j, 1)
            if (vof_correction > 0d0) then
                vof_factor = 1d0 / (1d0 - vof_correction)
                do tmp_mat = 1, this%n_materials
                    mat_vof(tmp_mat,i, j, 1) = mat_vof(tmp_mat,i, j, 1) * vof_factor
                end do
                vof(i, j, 1) = vof(i, j, 1) * vof_factor
            end if
        end if
        return
    end subroutine Put_quantities_into_cell

    subroutine Calculate_cell_quantities_in_advect(this, i, j, k)
        implicit none


        class(advect_t), intent(in out) :: this

        integer, intent(in) :: i, j, k

        real(8), dimension(:, :, :), pointer :: sie_adv
        real(8), dimension(:, :, :), pointer :: cell_mass_adv
        real(8), dimension(:, :, :), pointer :: vof_adv

        real(8), dimension(:, :, :, :), pointer :: sie_vof_adv
        real(8), dimension(:, :, :, :), pointer :: mat_cell_mass_adv
        real(8), dimension(:, :, :, :), pointer :: mat_vof_adv

        real(8), dimension(:, :, :), pointer :: sie
        real(8), dimension(:, :, :), pointer :: vol
        real(8), dimension(:, :, :), pointer :: vof
        real(8), dimension(:, :, :), pointer :: cell_mass

        real(8), dimension(:, :, :, :), pointer :: mat_vof
        real(8), dimension(:, :, :, :), pointer :: sie_vof
        real(8), dimension(:, :, :, :), pointer :: density_vof
        real(8), dimension(:, :, :, :), pointer :: cell_mass_vof

        real(8), dimension(:, :, :), pointer :: mat_vol

        real(8), dimension(:, :, :, :), pointer :: init_mat_layers

        real(8), dimension(:, :, :, :), pointer :: init_mat_layers_adv

        real(8), dimension(:, :, :), pointer :: momentum_x
        real(8), dimension(:, :, :), pointer :: momentum_y
        real(8), dimension(:, :, :), pointer :: momentum_x_adv
        real(8), dimension(:, :, :), pointer :: momentum_y_adv

        real(8), dimension(:, :, :, :), pointer :: area_t_in

        real(8), dimension(:), allocatable :: mat_cell_mass_adv_tmp
        real(8), dimension(:), allocatable :: sie_vof_adv_tmp
        real(8), dimension(:), allocatable :: mat_vof_adv_tmp
        real(8), dimension(:), allocatable :: init_mat_layers_tmp

        real(8) :: area_r_in
        real(8) :: area_l_in
        real(8) :: area_b_in
        real(8) :: area_out

        real(8) :: cell_mass_adv_right
        real(8) :: cell_mass_adv_top
        real(8) :: cell_mass_adv_left
        real(8) :: cell_mass_adv_bottom
        real(8) :: cell_mass_adv_out

        real(8) :: total_sie_adv
        real(8) :: total_cell_mass_adv
        real(8) :: total_area_adv

        integer :: tmp_mat

        call this%total_vof       %Point_to_data(vof)
        call this%total_cell_mass %Point_to_data(cell_mass)
        call this%total_sie       %Point_to_data(sie)

        call this%vof_advect      %Point_to_data(vof_adv)
        call this%adv_cell_mass%Point_to_data(cell_mass_adv)
        call this%adv_sie         %Point_to_data(sie_adv)

        call this%volume                %Point_to_data(vol)
        call this%rezone%material_volume%Point_to_data(mat_vol)

        call this%momentum%Point_to_data(momentum_x_adv, momentum_y_adv)
        call this%momentum%Point_to_data(momentum_x    , momentum_y    )


        allocate(mat_cell_mass_adv_tmp(this%n_materials))
        allocate(sie_vof_adv_tmp      (this%n_materials))
        allocate(mat_vof_adv_tmp      (this%n_materials))
        allocate(init_mat_layers_tmp  (this%n_materials))

        total_cell_mass_adv = 0d0
        total_area_adv      = 0d0
        total_sie_adv       = 0d0

        cell_mass_adv_right  = 0d0
        cell_mass_adv_top    = 0d0
        cell_mass_adv_left   = 0d0
        cell_mass_adv_bottom = 0d0
        cell_mass_adv_out    = 0d0

        call this%materials%density%Point_to_data(density_vof)
        call this%materials%sie    %Point_to_data(sie_vof)
        call this%adv_mats%area_top_in%Point_to_data(area_t_in)
        call this%materials%Point_to_initial_layers(init_mat_layers)


        call this%materials%cell_mass%Point_to_data(cell_mass_vof)
        call this%materials%vof      %Point_to_data(mat_vof)
        call this%adv_mats %cell_mass%Point_to_data(mat_cell_mass_adv)
        call this%adv_mats %sie      %Point_to_data(sie_vof_adv)
        call this%adv_mats %vof      %Point_to_data(mat_vof_adv)
        call this%adv_mats%Point_to_initial_layers(init_mat_layers_adv)


        do tmp_mat = 1, this%n_materials

            area_b_in  = this%adv_mats%area_bottom_in(tmp_mat)
            area_r_in  = this%adv_mats%area_right_in(tmp_mat)
            area_l_in  = this%adv_mats%area_left_in(tmp_mat)
            area_out   = this%adv_mats%total_out(tmp_mat)
            mat_cell_mass_adv_tmp(tmp_mat) = area_r_in * density_vof(tmp_mat, i+1, j, 1) + area_t_in(tmp_mat, i, 1 ,1) * density_vof(tmp_mat, i, j+1, 1) + &
                area_l_in * density_vof(tmp_mat, i-1, j, 1) + area_b_in          * density_vof(tmp_mat, i, j-1, 1) + &
                area_out  * density_vof(tmp_mat, i  , j, 1)

            sie_vof_adv_tmp(tmp_mat) = area_r_in          * density_vof(tmp_mat, i+1, j  , 1) * sie_vof(tmp_mat, i+1, j  , 1) + &
                area_t_in(tmp_mat, i, 1 ,1) * density_vof(tmp_mat, i  , j+1, 1) * sie_vof(tmp_mat, i  , j+1, 1) + &
                area_l_in          * density_vof(tmp_mat, i-1, j  , 1) * sie_vof(tmp_mat, i-1, j  , 1) + &
                area_b_in          * density_vof(tmp_mat, i  , j-1, 1) * sie_vof(tmp_mat, i  , j-1, 1) + &
                area_out           * density_vof(tmp_mat, i  , j  , 1) * sie_vof(tmp_mat, i  , j  , 1)

            mat_vof_adv_tmp(tmp_mat)  = area_r_in + area_l_in + area_t_in(tmp_mat, i, 1 ,1) + area_b_in + area_out


            if (this%advect_init_layer_mat > 0) then
                init_mat_layers_tmp(tmp_mat) = area_r_in          * density_vof(tmp_mat,i+1, j  , 1) * init_mat_layers(tmp_mat,i+1, j  , 1) + &
                    area_t_in(tmp_mat,i, 1 ,1) * density_vof(tmp_mat,i  , j+1, 1) * init_mat_layers(tmp_mat,i  , j+1, 1) + &
                    area_l_in          * density_vof(tmp_mat,i-1, j  , 1) * init_mat_layers(tmp_mat,i-1, j  , 1) + &
                    area_b_in          * density_vof(tmp_mat,i  , j-1, 1) * init_mat_layers(tmp_mat,i  , j-1, 1) + &
                    area_out           * density_vof(tmp_mat,i  , j  , 1) * init_mat_layers(tmp_mat,i  , j  , 1)
            else
                init_mat_layers_tmp(tmp_mat) = 0d0
            end if

            total_cell_mass_adv  = total_cell_mass_adv + mat_cell_mass_adv_tmp  (tmp_mat)
            total_area_adv       = total_area_adv      + mat_vof_adv_tmp (tmp_mat)
            total_sie_adv        = total_sie_adv       + sie_vof_adv_tmp(tmp_mat)

            cell_mass_adv_top    = cell_mass_adv_top    + area_t_in(tmp_mat,i, 1 ,1) * density_vof(tmp_mat,i  , j+1, 1)
            cell_mass_adv_right  = cell_mass_adv_right  + area_r_in          * density_vof(tmp_mat,i+1, j  , 1)
            cell_mass_adv_left   = cell_mass_adv_left   + area_l_in          * density_vof(tmp_mat,i-1, j  , 1)
            cell_mass_adv_bottom = cell_mass_adv_bottom + area_b_in          * density_vof(tmp_mat,i  , j-1, 1)
            cell_mass_adv_out    = cell_mass_adv_out    + area_out           * density_vof(tmp_mat,i  , j  , 1)
        end do


        vof_adv(i, j, 1)       = (vof(i, j, 1) * mat_vol(i, j, 1) + total_area_adv) / vol(i, j, 1)
        cell_mass_adv(i, j, 1) = cell_mass(i, j, 1) + total_cell_mass_adv

        do tmp_mat = 1, this%n_materials


            mat_vof_adv(tmp_mat,i, j, 1)       = (mat_vof(tmp_mat,i, j, 1) * mat_vol(i, j, 1) + mat_vof_adv_tmp(tmp_mat)) / vol(i, j, 1)
            mat_cell_mass_adv(tmp_mat,i ,j, 1) = cell_mass_vof(tmp_mat,i, j, 1) + mat_cell_mass_adv_tmp(tmp_mat)

            if (mat_vof_adv(tmp_mat,i, j, 1) > this%emf) then
                sie_vof_adv(tmp_mat,i, j, 1) = (sie_vof(tmp_mat,i, j, 1) * cell_mass_vof(tmp_mat,i, j, 1) + sie_vof_adv_tmp(tmp_mat)) / &
                    (mat_cell_mass_adv(tmp_mat,i ,j, 1) + 1d-30)


                if (this%advect_init_layer_mat > 0) then
                    init_mat_layers_adv(tmp_mat,i, j, 1) = (init_mat_layers(tmp_mat,i, j, 1) * cell_mass_vof(tmp_mat,i, j, 1) + init_mat_layers_tmp(tmp_mat)) / &
                        (mat_cell_mass_adv(tmp_mat,i ,j, 1) + 1d-30)
                end if



            end if
        end do

        if ((vof_adv(i, j, 1) > this%emf) .and. (cell_mass_adv(i, j, 1) > this%emfm)) then
            sie_adv(i, j, 1) = (cell_mass(i, j, 1) * sie(i, j, 1) + total_sie_adv) / cell_mass_adv(i, j, 1)
        else
            sie_adv(i, j, 1) = 0d0

        end if
        momentum_x_adv(i, j, 1) = cell_mass_adv_out   * momentum_x(i  , j, 1) + &
            cell_mass_adv_right * momentum_x(i+1, j, 1) + cell_mass_adv_top    * momentum_x(i, j+1, 1) + &
            cell_mass_adv_left  * momentum_x(i-1, j, 1) + cell_mass_adv_bottom * momentum_x(i, j-1, 1)
        momentum_y_adv(i, j, 1) = cell_mass_adv_out   * momentum_y(i  , j, 1) + &
            cell_mass_adv_right * momentum_y(i+1, j, 1) + cell_mass_adv_top    * momentum_y(i, j+1, 1) + &
            cell_mass_adv_left  * momentum_y(i-1, j, 1) + cell_mass_adv_bottom * momentum_y(i, j-1, 1)
        deallocate(mat_cell_mass_adv_tmp)
        deallocate(sie_vof_adv_tmp      )
        deallocate(mat_vof_adv_tmp      )
        deallocate(init_mat_layers_tmp  )

        return
    end subroutine Calculate_cell_quantities_in_advect

    subroutine Advect_velocities(this)
        use geometry_module, only : Triangle_area
        use general_utils_module, only : Int2d
        implicit none

        class(advect_t), intent(in out) :: this

        real(8), dimension(:, :, :), pointer :: x
        real(8), dimension(:, :, :), pointer :: y
        real(8), dimension(:, :, :), pointer :: material_x
        real(8), dimension(:, :, :), pointer :: material_y
        real(8), dimension(:, :, :), pointer :: density
        real(8), dimension(:, :, :), pointer :: velocity_x
        real(8), dimension(:, :, :), pointer :: velocity_y
        real(8), dimension(:, :, :), pointer :: momentum_x
        real(8), dimension(:, :, :), pointer :: momentum_y
        real(8), dimension(:, :, :), pointer :: vertex_mass_adv

        integer :: i, j, k
        integer :: l
        integer :: iv1, iv2, jv1, jv2
        integer :: i1, j1, i2, j2
        integer :: im, jm, ip, jp
        integer :: iv1m, iv1p, jv1m, jv1p, iv2m, iv2p, jv2m, jv2p

        real(8) :: area1, area2
        real(8) :: factor, factor1, factor2, factor3
        real(8) :: velocity_sq

        integer, dimension(4) :: i_2_add = (/1,  0, -1,  0/)
        integer, dimension(4) :: j_2_add = (/0,  1,  0, -1/)
        integer, dimension(4) :: i_3_add = (/0, -1,  0,  1/)
        integer, dimension(4) :: j_3_add = (/1,  0, -1,  0/)

        real(8), dimension(4) :: xc1, yc1, wc1

        call this%mesh           %Point_to_data (x, y)
        call this%rezone         %Point_to_coordinates_2d(material_x, material_y)
        call this%total_density  %Point_to_data (density)
        call this%velocity       %Point_to_data(velocity_x, velocity_y)
        call this%momentum       %Point_to_data(momentum_x, momentum_y)
        call this%adv_vertex_mass%Point_to_data (vertex_mass_adv)



        call this%vel_adv_weight%Exchange_virtual_space_blocking()

        do j = 1, this%nyp
            do i = 1, this%nxp
                l = 0
                area1 = -1d0
                area2 = -1d0
                do while((l < 4) .and. ((area1 < 0d0) .or. (area2 < 0d0)))
                    l = l + 1
                    i1 = i + i_2_add(l)
                    j1 = j + j_2_add(l)
                    i2 = i + i_3_add(l)
                    j2 = j + j_3_add(l)
                    area1 = Triangle_area(material_x(i, j, 1), material_y(i, j, 1), material_x(i1, j1, 1), material_y(i1, j1, 1), &
                        x(i, j, 1), y(i, j, 1))
                    area2 = Triangle_area(material_x(i, j, 1), material_y(i, j, 1), x(i, j, 1), y(i, j, 1),  &
                        material_x(i2, j2, 1), material_y(i2, j2, 1))
                end do

                iv1 = i + i_2_add(l)
                if (iv1 > this%nxp) iv1 = this%nx
                if (iv1 < 1) iv1 = 2

                iv2 = i + i_3_add(l)
                if (iv2 > this%nxp) iv2 = this%nx
                if (iv2 < 1) iv2 = 2
                jv1 = j + j_2_add(l)
                jv2 = j + j_3_add(l)


                if (jv1 > this%nyp) jv1 = this%ny
                if (jv1 < 1       ) jv1 = 2
                if (jv2 > this%nyp) jv2 = this%ny
                if (jv2 < 1       ) jv2 = 2

                xc1(1) = material_x(i  , j  , 1)
                xc1(2) = material_x(iv1, jv1, 1)
                xc1(3) = material_x(iv2, jv2, 1)
                yc1(1) = material_y(i  , j  , 1)
                yc1(2) = material_y(iv1, jv1, 1)
                yc1(3) = material_y(iv2, jv2, 1)
                im     = max(1, i - 1)
                ip     = min(this%nx, i)
                iv1m   = max(1, iv1 - 1)
                iv1p   = min(this%nx, iv1)
                iv2m   = max(1, iv2 - 1)
                iv2p   = min(this%nx, iv2)


                jm      = max(1, j - 1)
                jp      = min(this%ny, j)
                jv1m    = max(1, jv1 - 1)
                jv1p    = min(this%ny, jv1)
                jv2m    = max(1, jv2 - 1)
                jv2p    = min(this%ny, jv2)




                wc1(1) = density(ip  , jp  , 1) + density(im  , jp  , 1) + density(im  , jm  , 1) + density(ip  , jm  , 1)
                wc1(2) = density(iv1p, jv1p, 1) + density(iv1m, jv1p, 1) + density(iv1m, jv1m, 1) + density(iv1p, jv1m, 1)
                wc1(3) = density(iv2p, jv2p, 1) + density(iv2m, jv2p, 1) + density(iv2m, jv2m, 1) + density(iv2p, jv2m, 1)
                call Int2d(xc1(1), yc1(1), xc1(2), yc1(2), xc1(3), yc1(3), x(i, j, 1), y(i, j, 1), wc1(1), wc1(2), wc1(3),  &
                    factor1, factor2, factor3)
                momentum_x(i, j, 1) = factor1 * velocity_x(i, j, 1) + factor2 * velocity_x(iv1, jv1, 1) &
                    + factor3 * velocity_x(iv2, jv2, 1)
                momentum_y(i, j, 1) = factor1 * velocity_y(i, j, 1) + factor2 * velocity_y(iv1, jv1, 1) &
                    + factor3 * velocity_y(iv2, jv2, 1)
            end do
        end do
        do j = 1, this%nyp
            do i = 1, this%nxp
                if (vertex_mass_adv(i, j, 1) < this%emfm) then
                    velocity_x(i, j, 1) = 0d0
                    velocity_y(i, j, 1) = 0d0
                else
                    velocity_x(i, j, 1) = momentum_x(i, j, 1)
                    velocity_y(i, j, 1) = momentum_y(i, j, 1)
                end if
            end do
        end do


        call this%momentum%Exchange_virtual_space_blocking()
        do j = 1, this%nyp
            do i = 1, this%nxp
                velocity_sq = velocity_y(i, j, 1) * velocity_y(i, j, 1) + velocity_x(i, j, 1) * velocity_x(i, j, 1)
                if ((vertex_mass_adv(i, j, 1) < this%rezone%mass_threshold) .and. (velocity_sq > this%rezone%velocity_limit)) then
                    factor = sqrt(this%rezone%velocity_limit / velocity_sq)
                    velocity_x(i, j, 1) = velocity_x(i, j, 1) * factor
                    velocity_y(i, j, 1) = velocity_y(i, j, 1) * factor
                end if
            end do
        end do
        call this%velocity%Apply_boundary(this%mesh%coordinates%data)
        return
    end subroutine Advect_velocities

    subroutine Fix_overflow_error(this, i, j, k)
        implicit none

        class(advect_t), intent(in out) :: this

        integer, intent(in) :: i, j, k

        real(8), dimension(:, :, :, :), pointer :: mat_vof
        real(8), dimension(:, :, :), pointer :: mat_vol
        real(8), dimension(:, :, :, :), pointer :: area_t_in
        real(8), dimension(:, :, :, :), pointer :: area_t_out_max
        !        real(8), dimension(:, :, :), pointer :: area_t_out_max

        logical :: need_fix

        real(8) :: vof_max

        integer :: tmp_mat
        integer :: mat_max
        integer :: mat_fix

        real(8) :: area_r_in
        real(8) :: area_l_in
        real(8) :: area_b_in
        real(8) :: area_in

        real(8) :: area_r_out
        real(8) :: area_l_out
        real(8) :: area_b_out
        real(8) :: area_out

        real(8) :: fac


        if (.not. this%fix_overflow) then
            return
        end if
        call this%materials%vof%Point_to_data(mat_vof)

        need_fix = .false.
        vof_max = 0d0
        do tmp_mat = 1, this%n_materials
            if (mat_vof(tmp_mat,i, j, 1) < -this%emf) then
                need_fix = .true.
                mat_fix  = tmp_mat
            end if
            if (mat_vof(tmp_mat,i, j, 1) > vof_max) then
                vof_max = mat_vof(tmp_mat,i, j, 1)
                mat_max = tmp_mat
            end if
        end do
        if (.not. need_fix) return



        if (vof_max == 0d0) then
            return
        end if



        call this%adv_mats%area_top_in%Point_to_data(area_t_in)
        call this%adv_mats%area_top_out%Point_to_data(area_t_out_max)

        area_b_in  = this%adv_mats%area_bottom_in(mat_fix)
        area_r_in  = this%adv_mats%area_right_in(mat_fix)
        area_l_in  = this%adv_mats%area_left_in(mat_fix)

        area_b_out = this%adv_mats%area_bottom_out(mat_fix)
        area_r_out = this%adv_mats%area_right_out(mat_fix)
        area_l_out = this%adv_mats%area_left_out(mat_fix)

        area_out   = this%adv_mats%total_out(mat_fix)

        area_in = area_r_in + area_t_in(mat_fix, i, 1 ,1) + area_l_in + area_b_in

        call this%materials%vof%Point_to_data(mat_vof)
        call this%rezone%material_volume%Point_to_data(mat_vol)

        fac = (mat_vof(mat_fix,i, j, 1) * mat_vol(i, j, 1) + area_in) / abs(area_out)

        if (i /= 1 .and. area_l_out < 0d0) then
            call this%Fix_overflow_for_passed_cell(i, j, k, .true.  , &
                this%adv_mats%area_left_out(mat_fix)  , this%adv_mats%area_left_out(mat_max)  , fac, mat_fix, mat_max)
        end if
        if (j /= 0 .and. area_b_out < 0d0) then
            call this%Fix_overflow_for_passed_cell(i, j, k, .false. , &
                this%adv_mats%area_bottom_out(mat_fix), this%adv_mats%area_bottom_out(mat_max), fac, mat_fix, mat_max)
        end if


        !        call this%adv_mats%Point_to_area_top_out(area_t_out_max)

        this%adv_mats%area_bottom_out(mat_max) = this%adv_mats%area_bottom_out(mat_max) + area_b_out * (1d0 - fac)
        this%adv_mats%area_right_out(mat_max)  = this%adv_mats%area_right_out(mat_max)  + area_r_out * (1d0 - fac)
        this%adv_mats%area_left_out(mat_max)   = this%adv_mats%area_left_out(mat_max)   + area_l_out * (1d0 - fac)
        this%adv_mats%total_out(mat_max)       = this%adv_mats%total_out(mat_max)       + area_out   * (1d0 - fac)

        area_t_out_max(mat_max,i, 1 ,1) = area_t_out_max(mat_max,i, 1 ,1) + area_t_out_max(mat_fix,i, 1 ,1) * (1d0 - fac)

        this%adv_mats%area_bottom_out(mat_fix) = this%adv_mats%area_bottom_out(mat_fix) * fac
        this%adv_mats%area_right_out(mat_fix)  = this%adv_mats%area_right_out (mat_fix) * fac
        this%adv_mats%area_left_out (mat_fix)  = this%adv_mats%area_left_out(mat_fix)   * fac
        this%adv_mats%total_out (mat_fix)      = this%adv_mats%total_out(mat_fix)       * fac

        area_t_out_max(mat_fix,i, 1 ,1)     = area_t_out_max(mat_fix,i, 1 ,1) * fac

        call this%Calculate_cell_quantities_in_advect(i, j, k)


        return
    end subroutine Fix_overflow_error

    subroutine Fix_overflow_for_passed_cell(this, i, j, k, acceptor_i_direction, area_out_fix, area_out_max, fac, mat_fix, mat_max)
        implicit none

        class(advect_t), intent(in out) :: this

        logical, intent(in) :: acceptor_i_direction

        integer, intent(in) :: i, j, k
        integer, intent(in) :: mat_fix
        integer, intent(in) :: mat_max

        real(8), intent(in) :: fac
        real(8), intent(in) :: area_out_fix
        real(8), intent(in) :: area_out_max

        real(8), dimension(:, :, :), pointer :: cell_mass_adv
        real(8), dimension(:, :, :), pointer :: sie_adv

        real(8), dimension(:, :, :, :), pointer :: cell_mass_vof_adv_fix
        real(8), dimension(:, :, :, :), pointer :: cell_mass_vof_adv_max
        real(8), dimension(:, :, :, :), pointer :: mat_vof_adv_fix
        real(8), dimension(:, :, :, :), pointer :: mat_vof_adv_max
        real(8), dimension(:, :, :, :), pointer :: sie_vof_adv_fix
        real(8), dimension(:, :, :, :), pointer :: sie_vof_adv_max

        real(8), dimension(:, :, :, :), pointer :: density_vof_fix
        real(8), dimension(:, :, :, :), pointer :: density_vof_max
        real(8), dimension(:, :, :, :), pointer :: sie_vof_fix
        real(8), dimension(:, :, :, :), pointer :: sie_vof_max

        real(8), dimension(:, :, :), pointer :: vol

        real(8), dimension(:, :, :, :), pointer :: init_mat_layers
        real(8), dimension(:, :, :, :), pointer :: init_mat_layers_adv

        real(8), dimension(:, :, :), pointer :: momentum_x
        real(8), dimension(:, :, :), pointer :: momentum_y
        real(8), dimension(:, :, :), pointer :: momentum_x_adv
        real(8), dimension(:, :, :), pointer :: momentum_y_adv

        integer :: ia, ja

        real(8) :: cell_mass_adv_temp
        real(8) :: mat_vof_adv_fix_temp
        real(8) :: mat_vof_adv_max_temp
        real(8) :: cell_mass_mat_vof_adv_fix_temp
        real(8) :: cell_mass_mat_vof_adv_max_temp
        real(8) :: area_out_fix_corrected
        real(8) :: area_out_fix_diff
        real(8) :: area_out_max_diff



        call this%momentum%Point_to_data(momentum_x, momentum_y)

        call this%volume  %Point_to_data(vol)

        call this%materials%density %Point_to_data(density_vof_fix)
        call this%materials%density %Point_to_data(density_vof_max)
        call this%materials%sie     %Point_to_data(sie_vof_fix)
        call this%materials%sie     %Point_to_data(sie_vof_max)

        call this%adv_mats%cell_mass%Point_to_data(cell_mass_vof_adv_fix)
        call this%adv_mats%cell_mass%Point_to_data(cell_mass_vof_adv_max)
        call this%adv_mats%vof      %Point_to_data(mat_vof_adv_fix)
        call this%adv_mats%vof      %Point_to_data(mat_vof_adv_max)
        call this%adv_mats%sie      %Point_to_data(sie_vof_adv_fix)
        call this%adv_mats%sie      %Point_to_data(sie_vof_adv_max)

        call this%adv_cell_mass           %Point_to_data(cell_mass_adv)
        call this%adv_sie                    %Point_to_data(sie_adv)

        call this%momentum%Point_to_data(momentum_x_adv, momentum_y_adv)
        call this%materials%Point_to_initial_layers(init_mat_layers)
        call this%adv_mats%Point_to_initial_layers(init_mat_layers_adv)

        area_out_fix_corrected = area_out_fix * fac
        area_out_fix_diff      = area_out_fix_corrected - area_out_fix
        area_out_max_diff      = -area_out_fix_diff

        if (acceptor_i_direction) then
            ia = i - 1
            ja = j
        else
            ia = i
            ja = j - 1
        end if

        mat_vof_adv_fix_temp = mat_vof_adv_fix(mat_fix, ia, ja, 1)
        mat_vof_adv_max_temp = mat_vof_adv_max(mat_max, ia, ja, 1)

        mat_vof_adv_fix(mat_fix, ia, ja, 1) = mat_vof_adv_fix(mat_fix, ia, ja, 1) + area_out_fix_diff / vol(ia, ja, 1)
        mat_vof_adv_max(mat_max, ia, ja, 1) = mat_vof_adv_max(mat_max, ia, ja, 1) + area_out_max_diff / vol(ia, ja, 1)

        cell_mass_mat_vof_adv_fix_temp = cell_mass_vof_adv_fix(mat_fix, ia, ja, 1)
        cell_mass_mat_vof_adv_max_temp = cell_mass_vof_adv_max(mat_max, ia, ja, 1)


        cell_mass_vof_adv_fix(mat_fix, ia, ja, 1) = cell_mass_vof_adv_fix(mat_fix, ia, ja, 1) + area_out_fix_diff * density_vof_fix(mat_fix, i, j ,1)
        cell_mass_vof_adv_max(mat_max, ia, ja, 1) = cell_mass_vof_adv_max(mat_max, ia, ja, 1) + area_out_max_diff * density_vof_max(mat_max, i, j ,1)

        cell_mass_adv_temp = cell_mass_adv(ia, ja, 1)
        cell_mass_adv(ia, ja, 1) = cell_mass_adv(ia, ja, 1) + area_out_fix_diff * density_vof_fix(mat_fix, i, j ,1) + &
            area_out_max_diff * density_vof_max(mat_max, i, j ,1)

        sie_vof_adv_fix(mat_fix, i, j, 1) = (sie_vof_adv_fix(mat_fix, i, j, 1) * cell_mass_mat_vof_adv_fix_temp + &
            area_out_fix_diff * density_vof_fix(mat_fix, i, j ,1) * sie_vof_fix(mat_fix, i, j, 1)) &
            / (cell_mass_vof_adv_fix(mat_fix, ia, ja, 1) + 1d-30)
        sie_vof_adv_max(mat_max, i, j, 1) = (sie_vof_adv_max(mat_max, i, j, 1) * cell_mass_mat_vof_adv_max_temp + &
            area_out_max_diff * density_vof_max(mat_max, i, j ,1) * sie_vof_max(mat_max, i, j, 1)) &
            / (cell_mass_vof_adv_max(mat_max, ia, ja, 1) + 1d-30)


        if (this%advect_init_layer_mat > 0) then

            init_mat_layers_adv(mat_fix,ia, ja, 1) = (init_mat_layers_adv(mat_fix,ia, ja, 1) * cell_mass_mat_vof_adv_fix_temp + &
                area_out_fix_diff * density_vof_fix(mat_fix, i, j ,1) * init_mat_layers(mat_fix,i, j, 1)) / &
                (cell_mass_vof_adv_fix(mat_fix, ia, ja, 1) + 1d-30)

            init_mat_layers_adv(mat_max, ia, ja, 1) = (init_mat_layers_adv(mat_max, ia, ja, 1) * cell_mass_mat_vof_adv_max_temp + &
                area_out_max_diff * density_vof_max(mat_max, i, j ,1) * init_mat_layers(mat_max, i, j, 1)) / &
                (cell_mass_vof_adv_max(mat_max, ia, ja, 1) + 1d-30)
        end if







        sie_adv(ia, ja, 1) = (sie_adv(ia, ja, 1) * cell_mass_adv_temp + &
            area_out_fix_diff * density_vof_fix(mat_fix,ia, ja ,1) * sie_vof_fix(mat_fix,ia, ja, 1) + &
            area_out_max_diff * density_vof_max(mat_max, ia, ja ,1) * sie_vof_max(mat_max, ia, ja, 1)) / cell_mass_adv(ia, ja, 1)



        momentum_x_adv(ia, ja, 1) = momentum_x_adv(ia, ja, 1) + &
            momentum_x(i, j, 1) * (area_out_fix_diff * density_vof_fix(mat_fix,i, j ,1) + &
            area_out_max_diff * density_vof_max(mat_max, i, j ,1))
        momentum_y_adv(ia, ja, 1) = momentum_y_adv(ia, ja, 1) + &
            momentum_x(i, j, 1) * (area_out_fix_diff * density_vof_fix(mat_fix,i, j ,1) + &
            area_out_max_diff * density_vof_max(mat_max, i, j ,1))
        return
    end subroutine Fix_overflow_for_passed_cell


    subroutine Calculate_advect_3d(this)
        use geometry_module, only : Tetrahederon_volume, Hexahedron_volume, Vertex_interp_3d
        implicit none
        class(advect_t), intent(in out) :: this

        real(8), dimension(:, :, :), pointer :: x
        real(8), dimension(:, :, :), pointer :: y
        real(8), dimension(:, :, :), pointer :: z
        real(8), dimension(:, :, :), pointer :: material_x
        real(8), dimension(:, :, :), pointer :: material_y
        real(8), dimension(:, :, :), pointer :: material_z
        real(8), dimension(:, :, :), pointer :: velocity_x
        real(8), dimension(:, :, :), pointer :: velocity_y
        real(8), dimension(:, :, :), pointer :: velocity_z

        real(8), allocatable, dimension(:, :, :) :: velocity_x_adv
        real(8), allocatable, dimension(:, :, :) :: velocity_y_adv
        real(8), allocatable, dimension(:, :, :) :: velocity_z_adv

        real(8), dimension(:, :, :), pointer :: vol
        real(8), dimension(:, :, :), pointer :: vof
        real(8), dimension(:, :, :), pointer :: sie
        real(8), dimension(:, :, :), pointer :: density
        real(8), dimension(:, :, :), pointer :: cell_mass
        real(8), dimension(:, :, :), pointer :: vertex_mass
        real(8), dimension(:, :, :), pointer :: n_materials_in_cell
        real(8), dimension(:, :, :), pointer :: mat_id

        real(8), dimension(:, :, :, :), pointer :: mat_vof
        real(8), dimension(:, :, :, :), pointer :: sie_vof
        real(8), dimension(:, :, :, :), pointer :: density_vof
        real(8), dimension(:, :, :, :), pointer :: cell_mass_vof

        real(8), dimension(:, :, :), pointer :: vof_adv
        real(8), dimension(:, :, :, :), pointer :: mat_vof_adv
        real(8), dimension(:, :, :, :), pointer :: sie_vof_adv
        real(8), dimension(:, :, :, :), pointer :: mat_cell_mass_adv
        real(8), dimension(:, :, :), pointer :: cell_mass_adv

        real(8), dimension(:, :, :, :), pointer :: init_mat_layers

        real(8), dimension(:, :, :, :), pointer :: init_mat_layers_adv

        real(8), dimension(:, :, :), pointer :: vel_adv_w

        integer :: i, j, k, tmp_mat
        integer :: i_start, i_end
        integer :: j_start, j_end
        integer :: k_start, k_end
        integer :: i1, j1, k1, i2, j2, k2
        integer :: i3, j3, k3, i4, j4, k4
        integer :: i_face, j_face, k_face
        integer :: face
        integer :: tet_counter
        integer :: is_iregular
        integer :: id, jd, kd
        integer :: ia, ja, ka
        integer :: tet_counter_min
        integer :: i_tet1, j_tet1, k_tet1
        integer :: i_tet2, j_tet2, k_tet2
        integer :: i_tet3, j_tet3, k_tet3

        real(8) :: emf1
        real(8) :: x1, x2, x3, x4
        real(8) :: y1, y2, y3, y4
        real(8) :: z1, z2, z3, z4
        real(8) :: x_lag1, x_lag2, x_lag3, x_lag4
        real(8) :: y_lag1, y_lag2, y_lag3, y_lag4
        real(8) :: z_lag1, z_lag2, z_lag3, z_lag4
        real(8) :: x_mid1, y_mid1, z_mid1
        real(8) :: x_mid2, y_mid2, z_mid2
        real(8) :: x_mid_lag1, x_mid_lag2
        real(8) :: y_mid_lag1, y_mid_lag2
        real(8) :: z_mid_lag1, z_mid_lag2
        real(8) :: sign_tet_vol1, sign_tet_vol2
        real(8) :: dx_avg, dy_avg, dz_avg
        real(8) :: hexa_vol1, hexa_vol2
        real(8) :: hexa_vol3, hexa_vol4
        real(8) :: adv_vol
        real(8) :: weight
        real(8) :: hexa_vol_quart
        real(8) :: material_tet_vol
        real(8) :: material_tet_vol_sign
        real(8) :: donnor_mass
        real(8) :: donnor_sie
        real(8) :: donnor_init_mat_layer
        real(8) :: vol_temp
        real(8),dimension(:,:,:),allocatable :: vof_correction
        real(8) :: vof_factor
        real(8) :: tet_vol_diff
        real(8) :: tet_vol_diff_min
        real(8) :: symmetry_factor
        real(8) :: weight0, weight1
        real(8) :: weight2, weight3
        real(8) :: vel_x_min, vel_y_min, vel_z_min
        real(8) :: vel_x_max, vel_y_max, vel_z_max
        real(8) :: vel_factor
        real(8) :: velocity_sq
        real(8) :: dvof
        real(8) :: delete1,delete2,delete3
        logical :: associated_mat_vof
        real(8), dimension(4) :: tet_vol

        integer, dimension(6) :: face_i = [-1,  0,  0, 1, 0, 0]
        integer, dimension(6) :: face_j = [ 0, -1,  0, 0, 1, 0]
        integer, dimension(6) :: face_k = [ 0,  0, -1, 0, 0, 1]

        integer, dimension(4, 6) :: face_vert_i = reshape([ 0, 0, 0, 0, &
            1, 1, 0, 0, &
            0, 1, 1, 0, &
            1, 1, 1, 1, &
            0, 0, 1, 1, &
            0, 1, 1, 0  ], [4, 6])
        integer, dimension(4, 6) :: face_vert_j = reshape([ 0, 1, 1, 0, &
            0, 0, 0, 0, &
            1, 1, 0, 0, &
            0, 1, 1, 0, &
            1, 1, 1, 1, &
            0, 0, 1, 1  ], [4, 6])
        integer, dimension(4, 6) :: face_vert_k = reshape([ 1, 1, 0, 0, &
            0, 1, 1, 0, &
            0, 0, 0, 0, &
            0, 0, 1, 1, &
            0, 1, 1, 0, &
            1, 1, 1, 1  ], [4, 6])

        integer, dimension(8) :: i_1_add = [ 1,  0, -1,  0,  1,  0, -1,  0 ]
        integer, dimension(8) :: j_1_add = [ 0,  1,  0, -1,  0,  1, 0 , -1 ]
        integer, dimension(8) :: k_1_add = [ 0,  0,  0,  0,  0,  0, 0 ,  0 ]
        integer, dimension(8) :: i_2_add = [ 0, -1,  0,  1,  0,  1, 0 , -1 ]
        integer, dimension(8) :: j_2_add = [ 1,  0, -1,  0, -1,  0, 1 ,  0 ]
        integer, dimension(8) :: k_2_add = [ 0,  0,  0,  0,  0,  0, 0 ,  0 ]
        integer, dimension(8) :: i_3_add = [ 0,  0,  0,  0,  0,  0, 0 ,  0 ]
        integer, dimension(8) :: j_3_add = [ 0,  0,  0,  0,  0,  0, 0 ,  0 ]
        integer, dimension(8) :: k_3_add = [ 1,  1,  1,  1, -1, -1, -1, -1 ]

        logical :: wall_x_top, wall_x_bot, wall_y_top, wall_y_bot, wall_z_top, wall_z_bot
        integer :: virt_i_1, virt_i_nxp,virt_j_1, virt_j_nyp,virt_k_1, virt_k_nzp
        integer :: virt_nxp, virt_nyp, virt_nzp

        virt_nxp = this%parallel_params%virt_nxp
        virt_nyp = this%parallel_params%virt_nyp
        virt_nzp = this%parallel_params%virt_nzp
        wall_x_top = this%parallel_params%is_wall_x_top
        wall_x_bot = this%parallel_params%is_wall_x_bot
        wall_y_top = this%parallel_params%is_wall_y_top
        wall_y_bot = this%parallel_params%is_wall_y_bot
        wall_z_top = this%parallel_params%is_wall_z_top
        wall_z_bot = this%parallel_params%is_wall_z_bot


        call this%rezone  %Point_to_coordinates_3d(material_x, material_y, material_z)
        call this%mesh    %Point_to_data          (         x,          y,          z)
        call this%velocity%Point_to_data          (velocity_x, velocity_y, velocity_z)

        call this%total_vof%Point_to_data(vof)
        call this%volume   %Point_to_data(vol)
        call this%mat_cells%Point_to_data(mat_id)
        call this%total_sie        %Point_to_data(sie)
        call this%total_density    %Point_to_data(density)
        call this%total_cell_mass  %Point_to_data(cell_mass)
        call this%adv_cell_mass  %Point_to_data(cell_mass_adv)
        call this%vertex_mass%Point_to_data(vertex_mass)
        call this%num_mat_cells%Point_to_data(n_materials_in_cell)
        call this%vof_advect%Point_to_data(vof_adv)
        call this%adv_mats%vof%Point_to_data(mat_vof_adv)
        call this%adv_mats%Point_to_initial_layers(init_mat_layers_adv)
        call this%adv_mats%cell_mass%Point_to_data(mat_cell_mass_adv)
        call this%adv_mats%sie%Point_to_data(sie_vof_adv)

        call this%materials%vof%Point_to_data(mat_vof)
        call this%materials%density%Point_to_data(density_vof)
        call this%materials%sie%Point_to_data(sie_vof)
        call this%materials%cell_mass%Point_to_data(cell_mass_vof)
        call this%materials%Point_to_initial_layers(init_mat_layers)

        allocate(vof_correction(this%nx, this%ny, this%nz))

        i_start = 0
        j_start = 0
        k_start = 0
        i_end = this%nxp
        j_end = this%nyp
        k_end = this%nzp

        if (this%parallel_params%is_parallel .eqv. .true.) then
            virt_i_1 = this%parallel_params%i_virt(1)
            virt_i_nxp = this%parallel_params%i_virt(this%nxp)

            virt_j_1 = this%parallel_params%j_virt(1)
            virt_j_nyp = this%parallel_params%j_virt(this%nyp)

            virt_k_1 = this%parallel_params%k_virt(1)
            virt_k_nzp = this%parallel_params%k_virt(this%nzp)

            if (virt_i_1 /= 1 .or. wall_x_bot .eqv. .false.) i_start = 1
            if (virt_j_1 /= 1 .or. wall_y_bot .eqv. .false.) j_start = 1
            if (virt_k_1 /= 1 .or. wall_z_bot .eqv. .false.) k_start = 1
            if (virt_i_nxp /= virt_nxp .or. wall_x_top .eqv. .false.) i_end = i_end - 1
            if (virt_j_nyp /= virt_nyp .or. wall_y_top .eqv. .false.) j_end = j_end - 1
            if (virt_k_nzp /= virt_nzp .or. wall_z_top .eqv. .false.) k_end = k_end - 1
        end if









        call this%line_calc_3d(0, i_start, i_end, j_start, j_end, k_start, k_end)
        call this%line_calc_3d(1, i_start, i_end, j_start, j_end, k_start, k_end)





        !        do tmp_mat=1, this%n_materials
        call this%adv_mats%a%Exchange_virtual_space_blocking()
        call this%adv_mats%b%Exchange_virtual_space_blocking()
        call this%adv_mats%c%Exchange_virtual_space_blocking()
        call this%adv_mats%side%Exchange_virtual_space_blocking()
        !        end do

        call this%num_mat_cells%Exchange_virtual_space_blocking()



        vof_adv = 0d0
        mat_vof_adv = 0d0


        do k = 0, this%nzp
            do j = 0, this%nyp
                do i = 0, this%nxp
                    do tmp_mat = 1, this%n_materials

                        init_mat_layers_adv(tmp_mat, i, j, k) = init_mat_layers(tmp_mat, i, j, k) * cell_mass_vof(tmp_mat, i, j, k)
                        mat_cell_mass_adv  (tmp_mat, i, j, k) = cell_mass_vof(tmp_mat, i, j, k)
                        sie_vof_adv        (tmp_mat, i, j, k) = sie_vof(tmp_mat, i, j, k) * cell_mass_vof(tmp_mat, i, j, k)
                    end do
                end do
            end do
        end do

        do k = 1, this%nz
            do j = 1, this%ny
                do i = 1, this%nx

                    if ((vof(i  , j, k  ) < this%emf) .and. (vof(i  , j-1, k) < this%emf) .and.  &
                        (vof(i  , j, k-1) < this%emf) .and. (vof(i  , j+1, k) < this%emf) .and.  &
                        (vof(i  , j, k+1) < this%emf) .and. (vof(i-1, j  , k) < this%emf) .and.  &
                        (vof(i+1, j, k  ) < this%emf)) then

                        cycle
                    end if
                    associated_mat_vof = .false.

                    do face = 1, 6

                        if (this%shorter_advect) then
                            if ((i > 1 .and. face == 1) .or. (j > 1 .and. face == 2) .or. (k > 1 .and. face == 3)) then
                                cycle
                            end if
                        end if
                        i_face = i + face_i(face)
                        j_face = j + face_j(face)
                        k_face = k + face_k(face)

                        i1 = i + face_vert_i(1, face)
                        j1 = j + face_vert_j(1, face)
                        k1 = k + face_vert_k(1, face)

                        i2 = i + face_vert_i(2, face)
                        j2 = j + face_vert_j(2, face)
                        k2 = k + face_vert_k(2, face)

                        i3 = i + face_vert_i(3, face)
                        j3 = j + face_vert_j(3, face)
                        k3 = k + face_vert_k(3, face)

                        i4 = i + face_vert_i(4, face)
                        j4 = j + face_vert_j(4, face)
                        k4 = k + face_vert_k(4, face)

                        x_lag1 = material_x(i1, j1, k1)
                        y_lag1 = material_y(i1, j1, k1)
                        z_lag1 = material_z(i1, j1, k1)

                        x_lag2 = material_x(i2, j2, k2)
                        y_lag2 = material_y(i2, j2, k2)
                        z_lag2 = material_z(i2, j2, k2)

                        x_lag3 = material_x(i3, j3, k3)
                        y_lag3 = material_y(i3, j3, k3)
                        z_lag3 = material_z(i3, j3, k3)

                        x_lag4 = material_x(i4, j4, k4)
                        y_lag4 = material_y(i4, j4, k4)
                        z_lag4 = material_z(i4, j4, k4)

                        x1 = x(i1, j1, k1)
                        y1 = y(i1, j1, k1)
                        z1 = z(i1, j1, k1)

                        x2 = x(i2, j2, k2)
                        y2 = y(i2, j2, k2)
                        z2 = z(i2, j2, k2)

                        x3 = x(i3, j3, k3)
                        y3 = y(i3, j3, k3)
                        z3 = z(i3, j3, k3)

                        x4 = x(i4, j4, k4)
                        y4 = y(i4, j4, k4)
                        z4 = z(i4, j4, k4)


                        is_iregular = 0
                        tet_vol(1) = Tetrahederon_volume(x1, y1, z1, x2, y2, z2, x3, y3, z3, x_lag2, y_lag2, z_lag2)
                        tet_vol(2) = Tetrahederon_volume(x2, y2, z2, x3, y3, z3, x1, y1, z1, x_lag3, y_lag3, z_lag3)
                        tet_vol(3) = Tetrahederon_volume(x3, y3, z3, x1, y1, z1, x2, y2, z2, x_lag1, y_lag1, z_lag1)

                        tet_counter = 1
                        sign_tet_vol1 = 0d0

                        do while(tet_counter < 4 .and. is_iregular == 0)
                            if (abs(sign_tet_vol1) < 1d-10 .and. abs(tet_vol(tet_counter)) > 1d-10) then
                                sign_tet_vol1  = sign(1d0, tet_vol(tet_counter))
                            else if (abs(sign_tet_vol1) > 1d-10 .and. abs(tet_vol(tet_counter)) > 1d-10) then
                                sign_tet_vol2 = sign(1d0, tet_vol(tet_counter))
                                if (sign_tet_vol2 / sign_tet_vol1 < 0d0) is_iregular = 1
                            end if
                            tet_counter = tet_counter + 1
                        end do

                        if (is_iregular == 0) then
                            tet_vol(1) = Tetrahederon_volume(x1, y1, z1, x4, y4, z4, x3, y3, z3, x_lag4, y_lag4, z_lag4)
                            tet_vol(2) = Tetrahederon_volume(x4, y4, z4, x3, y3, z3, x1, y1, z1, x_lag3, y_lag3, z_lag3)
                            tet_vol(3) = Tetrahederon_volume(x3, y3, z3, x1, y1, z1, x4, y4, z4, x_lag1, y_lag1, z_lag1)
                            sign_tet_vol1 = 0d0
                            tet_counter = 1
                            do while(tet_counter < 4 .and. is_iregular == 0)
                                if      (abs(sign_tet_vol1) < 1d-10 .and. abs(tet_vol(tet_counter)) > 1d-10) then
                                    sign_tet_vol1  = sign(1d0, tet_vol(tet_counter))
                                else if (abs(sign_tet_vol1) > 1d-10 .and. abs(tet_vol(tet_counter)) > 1d-10) then
                                    sign_tet_vol2 = sign(1d0, tet_vol(tet_counter))
                                    if (sign_tet_vol2 / sign_tet_vol1 < 0d0) is_iregular = 1
                                end if
                                tet_counter = tet_counter + 1
                            end do
                        end if
                        if(is_iregular == 0) then
                            tet_vol(1) = Tetrahederon_volume(x1, y1, z1, x2, y2, z2, x4, y4, z4, x_lag2, y_lag2, z_lag2)
                            tet_vol(2) = Tetrahederon_volume(x2, y2, z2, x4, y4, z4, x1, y1, z1, x_lag4, y_lag4, z_lag4)
                            tet_vol(3) = Tetrahederon_volume(x4, y4, z4, x1, y1, z1, x2, y2, z2, x_lag1, y_lag1, z_lag1)
                            sign_tet_vol1 = 0d0
                            tet_counter = 1
                            do while(tet_counter < 4 .and. is_iregular == 0)
                                if      (abs(sign_tet_vol1) < 1d-10 .and. abs(tet_vol(tet_counter)) > 1d-10) then
                                    sign_tet_vol1  = sign(1d0, tet_vol(tet_counter))
                                else if (abs(sign_tet_vol1) > 1d-10 .and. abs(tet_vol(tet_counter)) > 1d-10) then
                                    sign_tet_vol2 = sign(1d0, tet_vol(tet_counter))
                                    if(sign_tet_vol2 / sign_tet_vol1 < 0d0) is_iregular = 1
                                end if
                                tet_counter = tet_counter + 1
                            end do
                        end if
                        if(is_iregular == 0) then
                            tet_vol(1) = Tetrahederon_volume(x4, y4, z4, x2, y2, z2, x3, y3, z3, x_lag2, y_lag2, z_lag2)
                            tet_vol(2) = Tetrahederon_volume(x2, y2, z2, x3, y3, z3, x4, y4, z4, x_lag3, y_lag3, z_lag3)
                            tet_vol(3) = Tetrahederon_volume(x3, y3, z3, x4, y4, z4, x2, y2, z2, x_lag4, y_lag4, z_lag4)
                            sign_tet_vol1 = 0d0
                            tet_counter = 1
                            do while(tet_counter < 4 .and. is_iregular == 0)
                                if      (abs(sign_tet_vol1) < 1d-10 .and. abs(tet_vol(tet_counter)) > 1d-10) then
                                    sign_tet_vol1  = sign(1d0, tet_vol(tet_counter))
                                else if (abs(sign_tet_vol1) > 1d-10 .and. abs(tet_vol(tet_counter)) > 1d-10) then
                                    sign_tet_vol2 = sign(1d0, tet_vol(tet_counter))
                                    if(sign_tet_vol2 / sign_tet_vol1 < 0d0) is_iregular = 1
                                end if
                                tet_counter = tet_counter + 1
                            end do
                        end if


                        if(is_iregular == 1) then
                            dx_avg = ((x_lag1 - x1) + (x_lag2 - x2) + (x_lag3 - x3) + (x_lag4 - x4)) / 4d0
                            dy_avg = ((y_lag1 - y1) + (y_lag2 - y2) + (y_lag3 - y3) + (y_lag4 - y4)) / 4d0
                            dz_avg = ((z_lag1 - z1) + (z_lag2 - z2) + (z_lag3 - z3) + (z_lag4 - z4)) / 4d0

                            x1 = x_lag1 - dx_avg
                            y1 = y_lag1 - dy_avg
                            z1 = z_lag1 - dz_avg

                            x2 = x_lag2 - dx_avg
                            y2 = y_lag2 - dy_avg
                            z2 = z_lag2 - dz_avg

                            x3 = x_lag3 - dx_avg
                            y3 = y_lag3 - dy_avg
                            z3 = z_lag3 - dz_avg

                            x4 = x_lag4 - dx_avg
                            y4 = y_lag4 - dy_avg
                            z4 = z_lag4 - dz_avg
                        end if



                        x_mid1  = (x1 + x3) / 2d0
                        y_mid1  = (y1 + y3) / 2d0
                        z_mid1  = (z1 + z3) / 2d0
                        x_mid2  = (x2 + x4) / 2d0
                        y_mid2  = (y2 + y4) / 2d0
                        z_mid2  = (z2 + z4) / 2d0

                        x_mid_lag1 = (x_lag1 + x_lag3) / 2d0
                        y_mid_lag1 = (y_lag1 + y_lag3) / 2d0
                        z_mid_lag1 = (z_lag1 + z_lag3) / 2d0
                        x_mid_lag2 = (x_lag2 + x_lag4) / 2d0
                        y_mid_lag2 = (y_lag2 + y_lag4) / 2d0
                        z_mid_lag2 = (z_lag2 + z_lag4) / 2d0

                        hexa_vol1 = Hexahedron_volume(x_lag1    , y_lag1    , z_lag1    , x_lag2    , y_lag2    , z_lag2    , &
                            x_lag3    , y_lag3    , z_lag3    , x_mid_lag1, y_mid_lag1, z_mid_lag1, &
                            x1        , y1        , z1        , x2        , y2        , z2        , &
                            x3        , y3        , z3        , x_mid1    , y_mid1    , z_mid1      )

                        hexa_vol2 = Hexahedron_volume(x_lag1    , y_lag1    , z_lag1    , x_mid_lag1, y_mid_lag1, z_mid_lag1, &
                            x_lag3    , y_lag3    , z_lag3    , x_lag4    , y_lag4    , z_lag4    , &
                            x1        , y1        , z1        , x_mid1    , y_mid1    , z_mid1    , &
                            x3        , y3        , z3        , x4        , y4        , z4          )

                        hexa_vol3 = Hexahedron_volume(x_lag1    , y_lag1    , z_lag1    , x_lag2    , y_lag2    , z_lag2    , &
                            x_mid_lag2, y_mid_lag2, z_mid_lag2, x_lag4    , y_lag4    , z_lag4    , &
                            x1        , y1        , z1        , x2        , y2        , z2        , &
                            x_mid2    , y_mid2    , z_mid2    , x4        , y4        , z4          )

                        hexa_vol4 = Hexahedron_volume(x_mid_lag2, y_mid_lag2, z_mid_lag2, x_lag2    , y_lag2    , z_lag2    , &
                            x_lag3    , y_lag3    , z_lag3    , x_lag4    , y_lag4    , z_lag4    , &
                            x_mid2    , y_mid2    , z_mid2    , x2        , y2        , z2        , &
                            x3        , y3        , z3        , x4        , y4        , z4          )

                        hexa_vol_quart = 0.25d0 * (hexa_vol1 + hexa_vol2 + hexa_vol3 + hexa_vol4)
                        weight = this%a0 * sign(1d0, hexa_vol_quart) + this%b0 * 2d0 * hexa_vol_quart / vol(i, j, k)

                        if (abs(hexa_vol_quart) > 1d-15) then
                            if (hexa_vol_quart < 0d0) then
                                id = i
                                jd = j
                                kd = k
                                ia = i_face
                                ja = j_face
                                ka = k_face
                            else
                                id = i_face
                                jd = j_face
                                kd = k_face
                                ia = i
                                ja = j
                                ka = k
                            end if
                            adv_vol = hexa_vol_quart * 2d0



                            call this%Volume_material_3d(i, j, k, id, jd, kd, ia, ja, ka               , &
                                hexa_vol_quart, weight                        , &
                                x_lag1, y_lag1, z_lag1, x_lag2, y_lag2, z_lag2, &
                                x_lag3, y_lag3, z_lag3, x_lag4, y_lag4, z_lag4, &
                                x1    , y1    , z1    , x2    , y2    , z2    , &
                                x3    , y3    , z3    , x4    , y4    , z4    , &
                                hexa_vol1, hexa_vol2, hexa_vol3, hexa_vol4, 1)

                            vof_adv(i, j, k) = vof_adv(i, j, k) + adv_vol
                            if (this%shorter_advect) vof_adv(i_face, j_face, k_face) = vof_adv(i_face, j_face, k_face) - adv_vol


                            do tmp_mat = 1, this%n_materials
                                dvof = this%adv_mats%fxtm(tmp_mat) * 2d0
                                mat_vof_adv(tmp_mat, i, j, k) = mat_vof_adv(tmp_mat, i, j, k) + dvof

                                if (this%shorter_advect) then
                                    mat_vof_adv(tmp_mat, i_face, j_face, k_face) = mat_vof_adv(tmp_mat, i_face, j_face, k_face) - dvof
                                end if
                            end do

                            do tmp_mat = 1, this%n_materials
                                dvof = this%adv_mats%fxtm(tmp_mat) * 2d0



                                donnor_mass = density_vof(tmp_mat, id, jd, kd) * dvof
                                donnor_sie  = sie_vof    (tmp_mat, id, jd, kd) * donnor_mass


                                if (abs(dvof) > vol(i,j,k)*this%emf /100d0 .and. tmp_mat /= mat_id(i,j,k)) then
                                    associated_mat_vof = .true.
                                    sie_vof_adv(tmp_mat, i, j, k) = sie_vof_adv(tmp_mat, i, j, k) + donnor_sie
                                    mat_cell_mass_adv(tmp_mat, i, j, k) = mat_cell_mass_adv(tmp_mat, i, j, k) + donnor_mass
                                else
                                    if (mat_id(i,j,k) == tmp_mat) then
                                        sie_vof_adv(tmp_mat, i, j, k) = sie_vof_adv(tmp_mat, i, j, k) + donnor_sie
                                        mat_cell_mass_adv(tmp_mat, i, j, k) = mat_cell_mass_adv(tmp_mat, i, j, k) + donnor_mass
                                    else
                                        if (n_materials_in_cell(i,j,k) > 1 .or. associated_mat_vof) then
                                            sie_vof_adv(tmp_mat, i, j, k) = sie_vof_adv(tmp_mat, i, j, k) + donnor_sie

                                            mat_cell_mass_adv(tmp_mat, i, j, k) = mat_cell_mass_adv(tmp_mat, i, j, k) + donnor_mass
                                        end if
                                    end if
                                end if


                                if (mat_id(i, j, k) == tmp_mat) then
                                    cell_mass_adv(i, j, k) = cell_mass_adv(i, j, k) + donnor_mass
                                end if

                                if (this%shorter_advect) then
                                    mat_cell_mass_adv(tmp_mat, i_face, j_face, k_face) = mat_cell_mass_adv(tmp_mat, i_face, j_face, k_face) - donnor_mass
                                end if





                                if (this%shorter_advect) then
                                    sie_vof_adv(tmp_mat, i_face, j_face, k_face) = sie_vof_adv(tmp_mat, i_face, j_face, k_face) - donnor_sie

                                end if

                                donnor_init_mat_layer = donnor_mass * init_mat_layers(tmp_mat, id, jd, kd)

                                init_mat_layers_adv(tmp_mat, i, j, k) = init_mat_layers_adv(tmp_mat, i, j, k) + donnor_init_mat_layer

                                if (this%shorter_advect) then
                                    init_mat_layers_adv(tmp_mat, i_face,j_face,k_face) =        init_mat_layers_adv(tmp_mat, i_face,j_face,k_face) - &
                                        donnor_init_mat_layer
                                end if

                            end do
                        end if
                    end do
                end do
            end do
        end do



        vof_correction = 0d0
        do k = 1, this%nz
            do j = 1, this%ny
                do i = 1, this%nx
                    do tmp_mat=1,this%n_materials

                        mat_vof(tmp_mat,i, j, k) = (mat_vof(tmp_mat,i, j, k) * (vol(i, j, k) - vof_adv(i, j, k)) + mat_vof_adv(tmp_mat,i, j, k)) / vol(i, j, k)
                        if (mat_vof(tmp_mat,i, j, k) > this%emf .and. mat_cell_mass_adv(tmp_mat, i, j, k) > this%emfm) then

                            sie_vof(tmp_mat,i, j, k) = max(0d0, sie_vof_adv(tmp_mat,i, j, k) / mat_cell_mass_adv(tmp_mat,i, j, k))
                            init_mat_layers(tmp_mat,i, j, k) = init_mat_layers_adv(tmp_mat,i, j, k) / mat_cell_mass_adv(tmp_mat,i, j, k)
                        else
                            vof_correction(i,j,k) = vof_correction(i,j,k) + mat_vof(tmp_mat,i, j, k)
                            mat_vof(tmp_mat,i, j, k) = 0d0
                            sie_vof(tmp_mat,i, j, k) = 0d0
                            init_mat_layers(tmp_mat,i, j, k) = 0d0
                        end if
                    end do
                end do
            end do
        end do

        do k = 1, this%nz
            do j = 1, this%ny
                do i = 1, this%nx
                    do tmp_mat = 1, this%n_materials

                        if (vof_correction(i,j,k) > 0d0) then
                            vof_factor = 1d0 / (1d0 - vof_correction(i,j,k))
                            mat_vof(tmp_mat,i, j, k) = mat_vof(tmp_mat,i, j, k) * vof_factor
                        end if
                    end do
                end do
            end do
        end do

        vof = 0d0
        mat_id = 0
        cell_mass = 0d0
        sie = 0d0
        n_materials_in_cell = 0d0
        do k = 1, this%nz
            do j = 1, this%ny
                do i = 1, this%nx
                    do tmp_mat=1,this%n_materials

                        if (mat_vof(tmp_mat,i, j, k) < this%emf .or. mat_cell_mass_adv(tmp_mat,i, j, k) < this%emfm) then
                            mat_cell_mass_adv(tmp_mat,i, j, k) = 0d0
                            sie_vof(tmp_mat,i, j, k) = 0d0
                            mat_vof(tmp_mat,i, j, k) = 0d0
                        else
                            cell_mass_vof(tmp_mat,i, j, k) = mat_cell_mass_adv(tmp_mat,i, j, k)
                            cell_mass    (i, j, k) = cell_mass(i, j, k) + mat_cell_mass_adv(tmp_mat,i, j, k)
                            vof          (i, j, k) = vof(i, j, k) + mat_vof(tmp_mat,i, j, k)
                            sie(i, j, k)           = sie(i, j, k) + sie_vof(tmp_mat,i, j, k) * cell_mass_vof(tmp_mat,i, j, k)
                            mat_id(i, j, k)        = mat_id(i, j, k) + tmp_mat * 10 ** n_materials_in_cell(i, j, k)
                            n_materials_in_cell(i, j, k) = n_materials_in_cell(i, j, k) + 1
                        end if
                    end do
                end do
            end do
        end do

        do k = 1, this%nz
            do j = 1, this%ny
                do i = 1, this%nx
                    if(cell_mass(i, j, k) > this%emfm) sie(i, j, k) = sie(i, j, k) / cell_mass(i, j, k)
                end do
            end do
        end do



        call this%vertex_mass%Calculate_vertex_mass_3d(this%mesh%coordinates, this%total_density, &
            this%total_cell_mass)
        call this%vertex_mass%Exchange_virtual_space_blocking()

        call this%materials%vof%Exchange_virtual_space_blocking()

        do k = 1, this%nzp
            do j = 1, this%nyp
                do i = 1, this%nxp
                    symmetry_factor = 1d0


                    this%vel_adv_weight%values(i, j, k) = (density(i  , j  , k  ) + density(i-1, j  , k  ) + &
                        density(i-1, j-1, k  ) + density(i  , j-1, k  ) + &
                        density(i  , j  , k-1) + density(i-1, j  , k-1) + &
                        density(i-1, j-1, k-1) + density(i  , j-1, k-1))

                end do
            end do
        end do

        velocity_x_adv = velocity_x
        velocity_y_adv = velocity_y
        velocity_z_adv = velocity_z


        call this%vel_adv_weight%Exchange_virtual_space_blocking()
        do k = 1, this%nzp
            do j = 1, this%nyp
                do i = 1, this%nxp
                    if (material_x(i, j, k) == x(i, j, k) .and. &
                        material_y(i, j, k) == y(i, j, k) .and. &
                        material_z(i, j, k) == z(i, j, k)) cycle
                    tet_counter=0
                    tet_vol = -1d0
                    tet_counter_min = 1
                    tet_vol_diff_min = 1d13

                    do while ((tet_vol(1) < 0d0) .or. (tet_vol(2) < 0d0) .or. (tet_vol(3) < 0d0))
                        tet_counter = tet_counter + 1
                        if(tet_counter == 9) then
                            tet_counter = tet_counter_min
                            tet_vol_diff_min = -1d0
                        end if
                        i1 = i + i_1_add(tet_counter)
                        j1 = j + j_1_add(tet_counter)
                        k1 = k + k_1_add(tet_counter)

                        i2 = i + i_2_add(tet_counter)
                        j2 = j + j_2_add(tet_counter)
                        k2 = k + k_2_add(tet_counter)

                        i3 = i + i_3_add(tet_counter)
                        j3 = j + j_3_add(tet_counter)
                        k3 = k + k_3_add(tet_counter)

                        material_tet_vol = Tetrahederon_volume(material_x(i ,  j, k ), material_y(i , j , k ), material_z(i , j , k ), &
                            material_x(i1, j1, k1), material_y(i1, j1, k1), material_z(i1, j1, k1), &
                            material_x(i2, j2, k2), material_y(i2, j2, k2), material_z(i2, j2, k2), &
                            material_x(i3, j3, k3), material_y(i3, j3, k3), material_z(i3, j3, k3))

                        material_tet_vol_sign = sign(1d0, material_tet_vol)

                        tet_vol(1)       = Tetrahederon_volume(material_x(i , j , k ), material_y(i , j , k ), material_z(i , j , k ), &
                            x(i , j , k ),          y(i , j , k ),          z(i , j , k ), &
                            material_x(i2, j2, k2), material_y(i2, j2, k2), material_z(i2, j2, k2), &
                            material_x(i3, j3, k3), material_y(i3, j3, k3), material_z(i3, j3, k3)) &
                            * material_tet_vol_sign

                        tet_vol(2)       = Tetrahederon_volume(material_x(i , j , k ), material_y(i , j , k ), material_z(i , j , k ), &
                            material_x(i1, j1, k1), material_y(i1, j1, k1), material_z(i1, j1, k1), &
                            x(i , j , k ),          y(i , j , k ),          z(i , j , k ), &
                            material_x(i3, j3, k3), material_y(i3, j3, k3), material_z(i3, j3, k3)) &
                            * material_tet_vol_sign

                        tet_vol(3)       = Tetrahederon_volume(material_x(i , j , k ), material_y(i , j , k ), material_z(i , j , k ), &
                            material_x(i1, j1, k1), material_y(i1, j1, k1), material_z(i1, j1, k1), &
                            material_x(i2, j2, k2), material_y(i2, j2, k2), material_z(i2, j2, k2), &
                            x(i , j , k ),          y(i , j , k ),          z(i , j , k )) &
                            * material_tet_vol_sign

                        tet_vol(4)       = Tetrahederon_volume(         x(i , j , k ),          y(i , j , k ),          z(i , j , k ), &
                            material_x(i1, j1, k1), material_y(i1, j1, k1), material_z(i1, j1, k1), &
                            material_x(i2, j2, k2), material_y(i2, j2, k2), material_z(i2, j2, k2), &
                            material_x(i3, j3, k3), material_y(i3, j3, k3), material_z(i3, j3, k3)) &
                            * material_tet_vol_sign

                        tet_vol_diff = abs(tet_vol(1) + tet_vol(2) + tet_vol(3) + tet_vol(4) - abs(material_tet_vol))
                        if (tet_vol_diff < tet_vol_diff_min) then
                            tet_counter_min = tet_counter
                            tet_vol_diff_min = tet_vol_diff
                        end if
                        if (tet_vol_diff_min < -0.5d0) then
                            tet_vol = 1d0
                        end if
                    end do


                    i_tet1 = i + i_1_add(tet_counter)
                    j_tet1 = j + j_1_add(tet_counter)
                    k_tet1 = k + k_1_add(tet_counter)

                    i_tet2 = i + i_2_add(tet_counter)
                    j_tet2 = j + j_2_add(tet_counter)
                    k_tet2 = k + k_2_add(tet_counter)

                    i_tet3 = i + i_3_add(tet_counter)
                    j_tet3 = j + j_3_add(tet_counter)
                    k_tet3 = k + k_3_add(tet_counter)

                    if (wall_x_top .eqv. .true.) then
                        if (i_tet1 > this%nxp) i_tet1 = this%nx
                        if (i_tet2 > this%nxp) i_tet2 = this%nx
                        if (i_tet3 > this%nxp) i_tet3 = this%nx
                    end if
                    if (wall_x_bot .eqv. .true.) then
                        if (i_tet1 < 1) i_tet1 = 2
                        if (i_tet2 < 1) i_tet2 = 2
                        if (i_tet3 < 1) i_tet3 = 2
                    end if
                    if (wall_y_top .eqv. .true.) then
                        if (j_tet1 > this%nyp) j_tet1 = this%ny
                        if (j_tet2 > this%nyp) j_tet2 = this%ny
                        if (j_tet3 > this%nyp) j_tet3 = this%ny
                    end if
                    if (wall_y_bot .eqv. .true.) then
                        if (j_tet1 < 1) j_tet1 = 2
                        if (j_tet2 < 1) j_tet2 = 2
                        if (j_tet3 < 1) j_tet3 = 2
                    end if
                    if (wall_z_top .eqv. .true.) then
                        if (k_tet1 > this%nzp) k_tet1 = this%nz
                        if (k_tet2 > this%nzp) k_tet2 = this%nz
                        if (k_tet3 > this%nzp) k_tet3 = this%nz
                    end if
                    if (wall_z_bot .eqv. .true.) then
                        if (k_tet1 < 1) k_tet1 = 2
                        if (k_tet2 < 1) k_tet2 = 2
                        if (k_tet3 < 1) k_tet3 = 2
                    end if
                    weight0 = this%vel_adv_weight%values(i     , j     , k     )
                    weight1 = this%vel_adv_weight%values(i_tet1, j_tet1, k_tet1)
                    weight2 = this%vel_adv_weight%values(i_tet2, j_tet2, k_tet2)
                    weight3 = this%vel_adv_weight%values(i_tet3, j_tet3, k_tet3)

                    if(weight0 < this%emf) then
                        velocity_x_adv(i, j, k) = 0d0
                        velocity_y_adv(i, j, k) = 0d0
                        velocity_z_adv(i, j, k) = 0d0
                    else
                        call Vertex_interp_3d(x(i, j, k), y(i, j, k), z(i, j, k), i, j, k                            , &
                            i_tet1, j_tet1, k_tet1, i_tet2, j_tet2, k_tet2, i_tet3, j_tet3, k_tet3 , &
                            velocity_x_adv(i, j, k), velocity_x, material_x, material_y, material_z, &
                            weight0, weight1, weight2, weight3)
                        call Vertex_interp_3d(x(i, j, k), y(i, j, k), z(i, j, k), i, j, k                            , &
                            i_tet1, j_tet1, k_tet1, i_tet2, j_tet2, k_tet2, i_tet3, j_tet3, k_tet3 , &
                            velocity_y_adv(i, j, k), velocity_y, material_x, material_y, material_z, &
                            weight0, weight1, weight2, weight3)
                        call Vertex_interp_3d(x(i, j, k), y(i, j, k), z(i, j, k), i, j, k                            , &
                            i_tet1, j_tet1, k_tet1, i_tet2, j_tet2, k_tet2, i_tet3, j_tet3, k_tet3 , &
                            velocity_z_adv(i, j, k), velocity_z, material_x, material_y, material_z, &
                            weight0, weight1, weight2, weight3)



                        vel_x_min = min(velocity_x(i_tet1, j_tet1, k_tet1), velocity_x(i_tet2, j_tet2, k_tet2), &
                            velocity_x(i_tet3, j_tet3, k_tet3), velocity_x(i     , j     , k     ))
                        vel_x_max = max(velocity_x(i_tet1, j_tet1, k_tet1), velocity_x(i_tet2, j_tet2, k_tet2), &
                            velocity_x(i_tet3, j_tet3, k_tet3), velocity_x(i     , j     , k     ))
                        vel_y_min = min(velocity_y(i_tet1, j_tet1, k_tet1), velocity_y(i_tet2, j_tet2, k_tet2), &
                            velocity_y(i_tet3, j_tet3, k_tet3), velocity_y(i     , j     , k     ))
                        vel_y_max = max(velocity_y(i_tet1, j_tet1, k_tet1), velocity_y(i_tet2, j_tet2, k_tet2), &
                            velocity_y(i_tet3, j_tet3, k_tet3), velocity_y(i     , j     , k     ))
                        vel_z_min = min(velocity_z(i_tet1, j_tet1, k_tet1), velocity_z(i_tet2, j_tet2, k_tet2), &
                            velocity_z(i_tet3, j_tet3, k_tet3), velocity_z(i     , j     , k     ))
                        vel_z_max = max(velocity_z(i_tet1, j_tet1, k_tet1), velocity_z(i_tet2, j_tet2, k_tet2), &
                            velocity_z(i_tet3, j_tet3, k_tet3), velocity_z(i     , j     , k     ))
                        velocity_x_adv(i, j, k) = max(velocity_x_adv(i, j, k), vel_x_min)
                        velocity_x_adv(i, j, k) = min(velocity_x_adv(i, j, k), vel_x_max)
                        velocity_y_adv(i, j, k) = max(velocity_y_adv(i, j, k), vel_y_min)
                        velocity_y_adv(i, j, k) = min(velocity_y_adv(i, j, k), vel_y_max)
                        velocity_z_adv(i, j, k) = max(velocity_z_adv(i, j, k), vel_z_min)
                        velocity_z_adv(i, j, k) = min(velocity_z_adv(i, j, k), vel_z_max)

                    end if
                end do
            end do
        end do

        velocity_x = velocity_x_adv
        velocity_y = velocity_y_adv
        velocity_z = velocity_z_adv

        do k = 1, this%nzp
            do j = 1, this%nyp
                do i = 1, this%nxp
                    velocity_sq = velocity_z(i, j, k) * velocity_z(i, j, k) + &
                        velocity_y(i, j, k) * velocity_y(i, j, k) + &
                        velocity_x(i, j, k) * velocity_x(i, j, k)
                    if ((vertex_mass(i, j, k) < this%rezone%mass_threshold) .and. (velocity_sq > this%rezone%velocity_limit)) then
                        vel_factor = sqrt(this%rezone%velocity_limit / velocity_sq)
                        velocity_x(i, j, k) = velocity_x(i, j, k) * vel_factor
                        velocity_y(i, j, k) = velocity_y(i, j, k) * vel_factor
                        velocity_z(i, j, k) = velocity_z(i, j, k) * vel_factor
                    end if
                end do
            end do
        end do

        call this%velocity%Apply_boundary(this%mesh%coordinates%data)

        deallocate(vof_correction)
        return
    end subroutine Calculate_advect_3d

    subroutine Line_calc_3d(this, nm, i_start, i_end, j_start, j_end, k_start, k_end)
        use geometry_module , only : Hexahedron_volume, Vector_grad,Vector_grad_planes, Vector_grad_vec, Volume_fraction_3d
        implicit none

        class(advect_t), intent(inout)    :: this
        integer        , intent(in)       :: nm, i_start, i_end, j_start, j_end, k_start, k_end

        real(8), dimension(:, :, :), pointer :: x_lag
        real(8), dimension(:, :, :), pointer :: y_lag
        real(8), dimension(:, :, :), pointer :: z_lag
        real(8), dimension(:, :, :), pointer :: vof
        real(8), dimension(:, :, :), pointer :: mat_id
        real(8), dimension(:, :, :), pointer :: n_materials_in_cell
        real(8), dimension(:, :, :, :), pointer :: a1
        real(8), dimension(:, :, :, :), pointer :: b1
        real(8), dimension(:, :, :, :), pointer :: c2
        real(8), dimension(:, :, :, :), pointer :: side1

        real(8), dimension(:, :, :, :), pointer :: mat_vof
        real(8), dimension(:, :, :, :), pointer :: density_vof
        real(8), dimension(:, :, :, :), pointer :: cell_mass_vof
        real(8) :: x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4, x5, y5, z5, x6, y6, z6, x7, y7, z7, x8, y8, z8, vijk, vipjk&
            , vipjpk, &
            vijpk, vijkp, vipjkp, vipjpkp, vijpkp, volijk, dfdx, dfdy, dfdz, venc, cmax, cmin, fac, cn, cold, v1, v2, delc,&
            c1, dvdc, vvmin, vvmax, vend
        integer :: im, jm, km, ip, jp, kp, nnii, nnic, tmp_mat, jj, kk, iv, jv, kv, i, j, k, fake_num_mats, numit, ii
        logical :: c_flag
        logical :: wall_x_top, wall_x_bot, wall_y_top, wall_y_bot, wall_z_top, wall_z_bot
        integer :: virt_nx, virt_ny, virt_nz

        real(8) :: u1u2u4, u1u4u5, u1u2u5, u2u3u1, u2u1u6, u2u3u6, u3u4u2, &
            u3u2u7, u3u4u7, u4u1u3, u4u3u8, u4u1u8, u5u8u6, u5u6u1, &
            u5u8u1, u6u5u7, u6u7u2, u6u5u2, u7u6u8, u7u8u3, u7u6u3, &
            u8u7u5, u8u5u4, u8u7u4

        real(8) :: y4y1z2z1z4z1y2y1, y5y1z4z1z5z1y4y1, y2y1z5z1z2z1y5y1, &
            y1y2z3z2z1z2y3y2, y6y2z1z2z6z2y1y2, y3y2z6z2z3z2y6y2, &
            y2y3z4z3z2z3y4y3, y7y3z2z3z7z3y2y3, y4y3z7z3z4z3y7y3, &
            y3y4z1z4z3z4y1y4, y8y4z3z4z8z4y3y4, y1y4z8z4z1z4y8y4, &
            y6y5z8z5z6z5y8y5, y1y5z6z5z1z5y6y5, y8y5z1z5z8z5y1y5, &
            y7y6z5z6z7z6y5y6, y2y6z7z6z2z6y7y6, y5y6z2z6z5z6y2y6, &
            y8y7z6z7z8z7y6y7, y3y7z8z7z3z7y8y7, y6y7z3z7z6z7y3y7, &
            y5y8z7z8z5z8y7y8, y4y8z5z8z4z8y5y8, y7y8z4z8z7z8y4y8

        real(8) :: z4z1x2x1x4x1z2z1, z5z1x4x1x5x1z4z1, z2z1x5x1x2x1z5z1, &
            z1z2x3x2x1x2z3z2, z6z2x1x2x6x2z1z2, z3z2x6x2x3x2z6z2, &
            z2z3x4x3x2x3z4z3, z7z3x2x3x7x3z2z3, z4z3x7x3x4x3z7z3, &
            z3z4x1x4x3x4z1z4, z8z4x3x4x8x4z3z4, z1z4x8x4x1x4z8z4, &
            z6z5x8x5x6x5z8z5, z1z5x6x5x1x5z6z5, z8z5x1x5x8x5z1z5, &
            z7z6x5x6x7x6z5z6, z2z6x7x6x2x6z7z6, z5z6x2x6x5x6z2z6, &
            z8z7x6x7x8x7z6z7, z3z7x8x7x3x7z8z7, z6z7x3x7x6x7z3z7, &
            z5z8x7x8x5x8z7z8, z4z8x5x8x4x8z5z8, z7z8x4x8x7x8z4z8

        real(8) :: x4x1y2y1y4y1x2x1, x5x1y4y1y5y1x4x1, x2x1y5y1y2y1x5x1, &
            x1x2y3y2y1y2x3x2, x6x2y1y2y6y2x1x2, x3x2y6y2y3y2x6x2, &
            x2x3y4y3y2y3x4x3, x7x3y2y3y7y3x2x3, x4x3y7y3y4y3x7x3, &
            x3x4y1y4y3y4x1x4, x8x4y3y4y8y4x3x4, x1x4y8y4y1y4x8x4, &
            x6x5y8y5y6y5x8x5, x1x5y6y5y1y5x6x5, x8x5y1y5y8y5x1x5, &
            x7x6y5y6y7y6x5x6, x2x6y7y6y2y6x7x6, x5x6y2y6y5y6x2x6, &
            x8x7y6y7y8y7x6x7, x3x7y8y7y3y7x8x7, x6x7y3y7y6y7x3x7, &
            x5x8y7y8y5y8x7x8, x4x8y5y8y4y8x5x8, x7x8y4y8y7y8x4x8


        virt_nx = this%parallel_params%virt_nx
        virt_ny = this%parallel_params%virt_ny
        virt_nz = this%parallel_params%virt_nz
        wall_x_top = this%parallel_params%is_wall_x_top
        wall_x_bot = this%parallel_params%is_wall_x_bot
        wall_y_top = this%parallel_params%is_wall_y_top
        wall_y_bot = this%parallel_params%is_wall_y_bot
        wall_z_top = this%parallel_params%is_wall_z_top
        wall_z_bot = this%parallel_params%is_wall_z_bot

        numit = 0
        call this%rezone   %Point_to_coordinates_3d (x_lag, y_lag, z_lag)
        call this%num_mat_cells%Point_to_data (n_materials_in_cell)


        if (nm == 0) then
            fake_num_mats = 1
        else
            fake_num_mats = this%n_materials
        end if

        !        do tmp_mat = 1, fake_num_mats
        if (nm /= 0) then
            call this%materials%vof%Point_to_data(mat_vof)
            call this%adv_mats%a%Point_to_data (a1)
            call this%adv_mats%b%Point_to_data (b1)
            call this%adv_mats%c%Point_to_data (c2)
            call this%adv_mats%side%Point_to_data(side1)
        else
            call this%total_vof%Point_to_data(vof)
            call this%a%Point_to_data (a1)
            call this%b%Point_to_data (b1)
            call this%c%Point_to_data (c2)
            call this%side%Point_to_data(side1)
        end if
        do k = k_start, k_end
            do j = j_start, j_end
                do i = i_start, i_end

                    if (nm /= 0) then
                        if( n_materials_in_cell(i, j, k) <= 1) cycle

                    else
                        if (vof(i, j, k) > this%emf1) cycle
                        if (vof(i,j,k) < this%emf) cycle
                    end if


                    do tmp_mat = 1, fake_num_mats
                        if (nm /= 0) then
                            if (mat_vof(tmp_mat, i,j,k) < this%emf) cycle
                        end if

                        c_flag = .true.
                        !                        if (nm /= 0 ) then
                        !                            if (n_materials_in_cell(i, j, k) <= 1) cycle
                        !                        else
                        !                            if (mat_vof(i, j, k) > this%emf1) then
                        !                                cycle
                        !                            end if
                        !                        end if

                        !                        if (mat_vof(i, j, k) < this%emf) cycle
                        ip = i + 1
                        im = i - 1
                        jp = j + 1
                        jm = j - 1
                        kp = k + 1
                        km = k - 1

                        x1 = x_lag(i, j, k)
                        y1 = y_lag(i, j, k)
                        z1 = z_lag(i, j, k)
                        x2 = x_lag(ip, j, k)
                        y2 = y_lag(ip, j, k)
                        z2 = z_lag(ip, j, k)
                        x3 = x_lag(ip, jp, k)
                        y3 = y_lag(ip, jp, k)
                        z3 = z_lag(ip, jp, k)
                        x4 = x_lag(i, jp, k)
                        y4 = y_lag(i, jp, k)
                        z4 = z_lag(i, jp, k)
                        x5 = x_lag(i, j, kp)
                        y5 = y_lag(i, j, kp)
                        z5 = z_lag(i, j, kp)
                        x6 = x_lag(ip, j, kp)
                        y6 = y_lag(ip, j, kp)
                        z6 = z_lag(ip, j, kp)
                        x7 = x_lag(ip, jp, kp)
                        y7 = y_lag(ip, jp, kp)
                        z7 = z_lag(ip, jp, kp)
                        x8 = x_lag(i, jp, kp)
                        y8 = y_lag(i, jp, kp)
                        z8 = z_lag(i, jp, kp)

                        volijk = Hexahedron_volume(x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4, &
                            x5, y5, z5, x6, y6, z6, x7, y7, z7, x8, y8, z8)

                        if (ip > this%nx .and. wall_x_top .eqv. .true.) ip = this%nx
                        if (im < 1       .and. wall_x_bot .eqv. .true.) im = 1
                        if (jp > this%ny .and. wall_y_top .eqv. .true.) jp = this%ny
                        if (jm < 1       .and. wall_y_bot .eqv. .true.) jm = 1
                        if (kp > this%nz .and. wall_z_top .eqv. .true.) kp = this%nz
                        if (km < 1       .and. wall_z_bot .eqv. .true.) km = 1

                        if (nm == 0) then
                            vijk    =  (vof(i ,  j, k)   + vof(im, j, k)   + vof(i, j, km)+ &
                                vof(im,  j, km)  + vof(i , jm, k)  + vof(im, jm, k)+ &
                                vof(i , jm, km)  + vof(im, jm, km)) / 8d0

                            vipjk   =  (vof(ip, j, k)   + vof(i, j, k)   + vof(ip, j, km)+ &
                                vof(i, j, km)   + vof(ip, jm, k) + vof(i, jm, k)+ &
                                vof(ip, jm, km) + vof(i, jm, km)) / 8d0

                            vipjpk  =  (vof(ip, jp, k)  + vof(i, jp, k)  + &
                                vof(ip, jp, km) + vof(i, jp, km) + vof(ip, j, k)+ &
                                vof(i, j, k)    + vof(ip, j, km) + vof(i, j, km))/ 8d0

                            vijpk   =  (vof(i, jp, k)   + vof(im, jp, k) + vof(i, jp, km) + &
                                vof(im, jp, km) + vof(i, j, k)   + vof(im, j, k)  + &
                                vof(i, j, km)   + vof(im, j, km))/ 8d0

                            vijkp   =  (vof(i, j, kp) + vof(im, j, kp) + vof(i, j, k)+ &
                                vof(im, j, k) + vof(i, jm, kp) + vof(im, jm, kp)+ &
                                vof(i, jm, k) + vof(im, jm, k))/ 8d0

                            vipjkp  =  (vof(ip, j, kp) + vof(i, j, kp)   + vof(ip, j, k)  + &
                                vof(i, j, k)   + vof(ip, jm, kp) + vof(i, jm, kp) + &
                                vof(ip, jm, k) + vof(i, jm, k))/ 8d0

                            vipjpkp =  (vof(ip, jp, kp) + vof(i, jp, kp) + &
                                vof(ip, jp, k)  + vof(i, jp, k)  + vof(ip, j, kp)+ &
                                vof(i, j, kp)   + vof(ip, j, k)  + vof(i, j, k))/ 8d0

                            vijpkp  =  (vof(i, jp, kp) + vof(im, jp, kp) + &
                                vof(i, jp, k)  + vof(im, jp, k)  + vof(i, j, kp)+ &
                                vof(im, j, kp) + vof(i, j, k)    + vof(im, j, k)) / 8d0
                                vend = vof(i, j, k) * volijk

                        else
                            vijk    =  (mat_vof(tmp_mat, i ,  j, k)   + mat_vof(tmp_mat, im, j, k)   + mat_vof(tmp_mat, i, j, km)+ &
                                mat_vof(tmp_mat, im,  j, km)  + mat_vof(tmp_mat, i , jm, k)  + mat_vof(tmp_mat, im, jm, k)+ &
                                mat_vof(tmp_mat, i , jm, km)  + mat_vof(tmp_mat, im, jm, km)) / 8d0

                            vipjk   =  (mat_vof(tmp_mat, ip, j, k)   + mat_vof(tmp_mat, i, j, k)   + mat_vof(tmp_mat, ip, j, km)+ &
                                mat_vof(tmp_mat, i, j, km)   + mat_vof(tmp_mat, ip, jm, k) + mat_vof(tmp_mat, i, jm, k)+ &
                                mat_vof(tmp_mat, ip, jm, km) + mat_vof(tmp_mat, i, jm, km)) / 8d0

                            vipjpk  =  (mat_vof(tmp_mat, ip, jp, k)  + mat_vof(tmp_mat, i, jp, k)  + &
                                mat_vof(tmp_mat, ip, jp, km) + mat_vof(tmp_mat, i, jp, km) + mat_vof(tmp_mat, ip, j, k)+ &
                                mat_vof(tmp_mat, i, j, k)    + mat_vof(tmp_mat, ip, j, km) + mat_vof(tmp_mat, i, j, km))/ 8d0

                            vijpk   =  (mat_vof(tmp_mat, i, jp, k)   + mat_vof(tmp_mat, im, jp, k) + mat_vof(tmp_mat, i, jp, km) + &
                                mat_vof(tmp_mat, im, jp, km) + mat_vof(tmp_mat, i, j, k)   + mat_vof(tmp_mat, im, j, k)  + &
                                mat_vof(tmp_mat, i, j, km)   + mat_vof(tmp_mat, im, j, km))/ 8d0

                            vijkp   =  (mat_vof(tmp_mat, i, j, kp) + mat_vof(tmp_mat, im, j, kp) + mat_vof(tmp_mat, i, j, k)+ &
                                mat_vof(tmp_mat, im, j, k) + mat_vof(tmp_mat, i, jm, kp) + mat_vof(tmp_mat, im, jm, kp)+ &
                                mat_vof(tmp_mat, i, jm, k) + mat_vof(tmp_mat, im, jm, k))/ 8d0

                            vipjkp  =  (mat_vof(tmp_mat, ip, j, kp) + mat_vof(tmp_mat, i, j, kp)   + mat_vof(tmp_mat, ip, j, k)  + &
                                mat_vof(tmp_mat, i, j, k)   + mat_vof(tmp_mat, ip, jm, kp) + mat_vof(tmp_mat, i, jm, kp) + &
                                mat_vof(tmp_mat, ip, jm, k) + mat_vof(tmp_mat, i, jm, k))/ 8d0

                            vipjpkp =  (mat_vof(tmp_mat, ip, jp, kp) + mat_vof(tmp_mat, i, jp, kp) + &
                                mat_vof(tmp_mat, ip, jp, k)  + mat_vof(tmp_mat, i, jp, k)  + mat_vof(tmp_mat, ip, j, kp)+ &
                                mat_vof(tmp_mat, i, j, kp)   + mat_vof(tmp_mat, ip, j, k)  + mat_vof(tmp_mat, i, j, k))/ 8d0

                            vijpkp  =  (mat_vof(tmp_mat, i, jp, kp) + mat_vof(tmp_mat, im, jp, kp) + &
                                mat_vof(tmp_mat, i, jp, k)  + mat_vof(tmp_mat, im, jp, k)  + mat_vof(tmp_mat, i, j, kp)+ &
                                mat_vof(tmp_mat, im, j, kp) + mat_vof(tmp_mat, i, j, k)    + mat_vof(tmp_mat, im, j, k)) / 8d0
                                vend = mat_vof(tmp_mat, i, j, k) * volijk
                        end if



                        call Vector_grad_planes(y1, y2, y3, y4, y5, y6, y7, y8, &
                            z1, z2, z3, z4, z5, z6, z7, z8, y4y1z2z1z4z1y2y1, y5y1z4z1z5z1y4y1, y2y1z5z1z2z1y5y1, &
                            y1y2z3z2z1z2y3y2, y6y2z1z2z6z2y1y2, y3y2z6z2z3z2y6y2, &
                            y2y3z4z3z2z3y4y3, y7y3z2z3z7z3y2y3, y4y3z7z3z4z3y7y3, &
                            y3y4z1z4z3z4y1y4, y8y4z3z4z8z4y3y4, y1y4z8z4z1z4y8y4, &
                            y6y5z8z5z6z5y8y5, y1y5z6z5z1z5y6y5, y8y5z1z5z8z5y1y5, &
                            y7y6z5z6z7z6y5y6, y2y6z7z6z2z6y7y6, y5y6z2z6z5z6y2y6, &
                            y8y7z6z7z8z7y6y7, y3y7z8z7z3z7y8y7, y6y7z3z7z6z7y3y7, &
                            y5y8z7z8z5z8y7y8, y4y8z5z8z4z8y5y8, y7y8z4z8z7z8y4y8)
                        call Vector_grad_planes(z1, z2, z3, z4, z5, z6, z7, z8, &
                            x1, x2, x3, x4, x5, x6, x7, x8, z4z1x2x1x4x1z2z1, z5z1x4x1x5x1z4z1, z2z1x5x1x2x1z5z1, &
                            z1z2x3x2x1x2z3z2, z6z2x1x2x6x2z1z2, z3z2x6x2x3x2z6z2, &
                            z2z3x4x3x2x3z4z3, z7z3x2x3x7x3z2z3, z4z3x7x3x4x3z7z3, &
                            z3z4x1x4x3x4z1z4, z8z4x3x4x8x4z3z4, z1z4x8x4x1x4z8z4, &
                            z6z5x8x5x6x5z8z5, z1z5x6x5x1x5z6z5, z8z5x1x5x8x5z1z5, &
                            z7z6x5x6x7x6z5z6, z2z6x7x6x2x6z7z6, z5z6x2x6x5x6z2z6, &
                            z8z7x6x7x8x7z6z7, z3z7x8x7x3x7z8z7, z6z7x3x7x6x7z3z7, &
                            z5z8x7x8x5x8z7z8, z4z8x5x8x4x8z5z8, z7z8x4x8x7x8z4z8)

                        call Vector_grad_planes(x1, x2, x3, x4, x5, x6, x7, x8, &
                            y1, y2, y3, y4, y5, y6, y7, y8, x4x1y2y1y4y1x2x1, x5x1y4y1y5y1x4x1, x2x1y5y1y2y1x5x1, &
                            x1x2y3y2y1y2x3x2, x6x2y1y2y6y2x1x2, x3x2y6y2y3y2x6x2, &
                            x2x3y4y3y2y3x4x3, x7x3y2y3y7y3x2x3, x4x3y7y3y4y3x7x3, &
                            x3x4y1y4y3y4x1x4, x8x4y3y4y8y4x3x4, x1x4y8y4y1y4x8x4, &
                            x6x5y8y5y6y5x8x5, x1x5y6y5y1y5x6x5, x8x5y1y5y8y5x1x5, &
                            x7x6y5y6y7y6x5x6, x2x6y7y6y2y6x7x6, x5x6y2y6y5y6x2x6, &
                            x8x7y6y7y8y7x6x7, x3x7y8y7y3y7x8x7, x6x7y3y7y6y7x3x7, &
                            x5x8y7y8y5y8x7x8, x4x8y5y8y4y8x5x8, x7x8y4y8y7y8x4x8)

                        call Vector_grad_vec(vijk, vipjk, vipjpk, vijpk, vijkp, vipjkp, vipjpkp, vijpkp,&
                            u1u2u4, u1u4u5, u1u2u5, &
                            u2u3u1, u2u1u6, u2u3u6, &
                            u3u4u2, u3u2u7, u3u4u7, &
                            u4u1u3, u4u3u8, u4u1u8, &
                            u5u8u6, u5u6u1, u5u8u1, &
                            u6u5u7, u6u7u2, u6u5u2, &
                            u7u6u8, u7u8u3, u7u6u3, &
                            u8u7u5, u8u5u4, u8u7u4)


                        dfdx = 1 / (12*volijk)*( &
                            (u1u2u4 * y4y1z2z1z4z1y2y1 + u1u4u5 * y5y1z4z1z5z1y4y1 + u1u2u5 * y2y1z5z1z2z1y5y1) +&
                            (u2u3u1 * y1y2z3z2z1z2y3y2 + u2u1u6 * y6y2z1z2z6z2y1y2 + u2u3u6 * y3y2z6z2z3z2y6y2) +&
                            (u3u4u2 * y2y3z4z3z2z3y4y3 + u3u2u7 * y7y3z2z3z7z3y2y3 + u3u4u7 * y4y3z7z3z4z3y7y3) +&
                            (u4u1u3 * y3y4z1z4z3z4y1y4 + u4u3u8 * y8y4z3z4z8z4y3y4 + u4u1u8 * y1y4z8z4z1z4y8y4) +&
                            (u5u8u6 * y6y5z8z5z6z5y8y5 + u5u6u1 * y1y5z6z5z1z5y6y5 + u5u8u1 * y8y5z1z5z8z5y1y5) +&
                            (u6u5u7 * y7y6z5z6z7z6y5y6 + u6u7u2 * y2y6z7z6z2z6y7y6 + u6u5u2 * y5y6z2z6z5z6y2y6) +&
                            (u7u6u8 * y8y7z6z7z8z7y6y7 + u7u8u3 * y3y7z8z7z3z7y8y7 + u7u6u3 * y6y7z3z7z6z7y3y7) +&
                            (u8u7u5 * y5y8z7z8z5z8y7y8 + u8u5u4 * y4y8z5z8z4z8y5y8 + u8u7u4 * y7y8z4z8z7z8y4y8))

                        dfdy = 1 / (12*volijk)*( &
                            (u1u2u4 * z4z1x2x1x4x1z2z1 + u1u4u5 * z5z1x4x1x5x1z4z1 + u1u2u5 * z2z1x5x1x2x1z5z1) +&
                            (u2u3u1 * z1z2x3x2x1x2z3z2 + u2u1u6 * z6z2x1x2x6x2z1z2 + u2u3u6 * z3z2x6x2x3x2z6z2) +&
                            (u3u4u2 * z2z3x4x3x2x3z4z3 + u3u2u7 * z7z3x2x3x7x3z2z3 + u3u4u7 * z4z3x7x3x4x3z7z3) +&
                            (u4u1u3 * z3z4x1x4x3x4z1z4 + u4u3u8 * z8z4x3x4x8x4z3z4 + u4u1u8 * z1z4x8x4x1x4z8z4) +&
                            (u5u8u6 * z6z5x8x5x6x5z8z5 + u5u6u1 * z1z5x6x5x1x5z6z5 + u5u8u1 * z8z5x1x5x8x5z1z5) +&
                            (u6u5u7 * z7z6x5x6x7x6z5z6 + u6u7u2 * z2z6x7x6x2x6z7z6 + u6u5u2 * z5z6x2x6x5x6z2z6) +&
                            (u7u6u8 * z8z7x6x7x8x7z6z7 + u7u8u3 * z3z7x8x7x3x7z8z7 + u7u6u3 * z6z7x3x7x6x7z3z7) +&
                            (u8u7u5 * z5z8x7x8x5x8z7z8 + u8u5u4 * z4z8x5x8x4x8z5z8 + u8u7u4 * z7z8x4x8x7x8z4z8))

                        dfdz = 1 / (12.*volijk)*( &
                            (u1u2u4 * x4x1y2y1y4y1x2x1 + u1u4u5 * x5x1y4y1y5y1x4x1 + u1u2u5 * x2x1y5y1y2y1x5x1) +&
                            (u2u3u1 * x1x2y3y2y1y2x3x2 + u2u1u6 * x6x2y1y2y6y2x1x2 + u2u3u6 * x3x2y6y2y3y2x6x2) +&
                            (u3u4u2 * x2x3y4y3y2y3x4x3 + u3u2u7 * x7x3y2y3y7y3x2x3 + u3u4u7 * x4x3y7y3y4y3x7x3) +&
                            (u4u1u3 * x3x4y1y4y3y4x1x4 + u4u3u8 * x8x4y3y4y8y4x3x4 + u4u1u8 * x1x4y8y4y1y4x8x4) +&
                            (u5u8u6 * x6x5y8y5y6y5x8x5 + u5u6u1 * x1x5y6y5y1y5x6x5 + u5u8u1 * x8x5y1y5y8y5x1x5) +&
                            (u6u5u7 * x7x6y5y6y7y6x5x6 + u6u7u2 * x2x6y7y6y2y6x7x6 + u6u5u2 * x5x6y2y6y5y6x2x6) +&
                            (u7u6u8 * x8x7y6y7y8y7x6x7 + u7u8u3 * x3x7y8y7y3y7x8x7 + u7u6u3 * x6x7y3y7y6y7x3x7) +&
                            (u8u7u5 * x5x8y7y8y5y8x7x8 + u8u5u4 * x4x8y5y8y4y8x5x8 + u8u7u4 * x7x8y4y8y7y8x4x8))

                        if (abs(dfdz) < 1d-8) dfdz = sign(1d-8, dfdz)
                        side1(tmp_mat, i, j, k) = sign(1d0, dfdz)
                        a1(tmp_mat, i, j, k) = -dfdx / dfdz
                        b1(tmp_mat, i, j, k) = -dfdy / dfdz



                        cmax = -1d20
                        cmin = 1d20
                        do kk = 0, 1
                            do jj = 0, 1
                                do ii = 0, 1
                                    cmax = max(cmax, z_lag(i + ii,j + jj, k + kk) - a1(tmp_mat, i, j, k)*x_lag(i + ii, j + jj, k + kk)- &
                                        b1(tmp_mat, i, j, k)*y_lag(i + ii, j + jj, k + kk))
                                    cmin = min(cmin, z_lag(i + ii, j + jj, k + kk) - a1(tmp_mat, i, j, k)*x_lag(i + ii, j + jj, k + kk)- &
                                        b1(tmp_mat, i, j, k)*y_lag(i + ii, j + jj, k + kk))
                                end do
                            end do
                        end do
                        if(cmax == cmin) then
                            c2(tmp_mat, i, j, k) = cmax
                        end if


                        nnii = 0
                        cn = 0.5d0*(cmin + cmax)
                        do
                            cold = cn
                            c2(tmp_mat, i, j, k) = cn
                            nnii = nnii + 1
                            if (nnii > 15) then
                                exit
                            end if
                            v1 = Volume_fraction_3d (x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4, &
                                x5, y5, z5, x6, y6, z6, x7, y7, z7, x8, y8, z8, &
                                a1(tmp_mat, i, j, k), b1(tmp_mat, i, j, k), c2(tmp_mat, i, j, k), side1(tmp_mat, i, j, k), volijk)
                            delc = 1d-7*(cmax - cmin)
                            if (c2(tmp_mat, i, j, k) + delc >= cmax) then
                                if (abs(c2(tmp_mat, i, j, k) - cmax) > 1d-8) then
                                    delc = 0.9d0*(cmax - c2(tmp_mat, i, j, k))
                                else
                                    delc = 0.1d0*(cmin - c2(tmp_mat, i, j, k))
                                end if
                            end if
                            if (c2(tmp_mat, i, j, k) + delc <= cmin) then
                                if (abs(c2(tmp_mat, i, j, k) - cmin) > 1d-8) then
                                    delc = 0.9d0*(cmin - c2(tmp_mat, i, j, k))
                                else
                                    delc = 0.1d0*(cmax - c2(tmp_mat, i, j, k))
                                end if
                            end if
                            c1 = c2(tmp_mat, i, j, k) + delc
                            v2 = Volume_fraction_3d (x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4, &
                                x5, y5, z5, x6, y6, z6, x7, y7, z7, x8, y8, z8, a1(tmp_mat, i, j, k), b1(tmp_mat, i, j, k), c1, side1(tmp_mat, i, j, k), volijk)
                            if (abs(delc) < 1d-20) then
                                exit
                            end if
                            dvdc = (v2 - v1) / delc
                            if (dvdc == 0) then
                                exit
                            end if
                            cn = c2(tmp_mat, i, j, k) - (v1 - vend) / dvdc
                            if (cn > cmax) then
                                cn = 2*cmax - cn
                                if (cn < cmin) cn = 0.5d0*(cmax + cmin)
                            else if (cn < cmin) then
                                cn = 2*cmin - cn
                                if (cn >= cmax) cn = 0.5d0*(cmax + cmin)
                            end if
                            if (abs(v1 - vend) > 1d-1 * this%emf * volijk) cycle
                            c2(tmp_mat, i, j, k) = cn
                            numit = max(nnii, numit)
                            c_flag = .false.
                            exit
                        end do
                        if (c_flag) then
                            vvmin = Volume_fraction_3d(x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4, &
                                x5, y5, z5, x6, y6, z6, x7, y7, z7, x8, y8, z8, &
                                a1(tmp_mat, i, j, k), b1(tmp_mat, i, j, k), cmin, side1(tmp_mat, i, j, k), volijk)

                            vvmax = Volume_fraction_3d(x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4, &
                                x5, y5, z5, x6, y6, z6, x7, y7, z7, x8, y8, z8, &
                                a1(tmp_mat, i, j, k), b1(tmp_mat, i, j, k), cmax, side1(tmp_mat, i, j, k), volijk)
                            nnii = 0
                            do
                                cn = 0.5d0*(cmin + cmax)
                                nnii = nnii + 1
                                if (nnii > 50) then
                                    c2(tmp_mat, i, j, k) = cn
                                    exit
                                end if
                                v1 = Volume_fraction_3d(x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4, &
                                    x5, y5, z5, x6, y6, z6, x7, y7, z7, x8, y8, z8, &
                                    a1(tmp_mat, i, j, k), b1(tmp_mat, i, j, k), cn, side1(tmp_mat, i, j, k), volijk)
                                if ((v1 < vend .and. vvmax > vvmin) .or.(v1 > vend .and. vvmax < vvmin)) then
                                    cmin = cn
                                else
                                    cmax = cn
                                end if
                                if (abs(v1 - vend) <= 1d-1*this%emf*volijk) then
                                    exit
                                end if
                            end do
                        end if
                        c2(tmp_mat, i, j, k) = cn
                    end do
                end do
            end do

        end do
        return
    end subroutine Line_calc_3d

    subroutine Point_to_material (this, material_ptr, mat_index)
        class (advect_t), intent(in out) :: this
        type(material_t), pointer, intent(out) :: material_ptr
        integer, intent(in) :: mat_index

!        material_ptr => this%materials(mat_index)

    end subroutine Point_to_material

    subroutine Point_to_adv_material (this, adv_material_ptr, mat_index)
        class (advect_t), intent(in out) :: this
        type(material_advect_t), pointer, intent(out) :: adv_material_ptr
        integer, intent(in) :: mat_index

!        adv_material_ptr => this%adv_mats(mat_index)

    end subroutine Point_to_adv_material

    subroutine Point_to_adv_vof_data (this, adv_vof_ptr)
        class (advect_t), intent(in out) :: this
        real(8), dimension(:,:,:), pointer, intent(out) :: adv_vof_ptr

!        call this%vof_advect%Point_to_data(adv_vof_ptr)

    end subroutine Point_to_adv_vof_data

    subroutine Point_to_adv_cell_mass_data (this, adv_cell_mass_ptr)
        class (advect_t), intent(in out) :: this
        real(8), dimension(:,:,:), pointer, intent(out) :: adv_cell_mass_ptr

!        call this%adv_cell_mass%Point_to_data(adv_cell_mass_ptr)

    end subroutine Point_to_adv_cell_mass_data

    subroutine Point_to_adv_sie_data (this, adv_sie_ptr)
        class (advect_t), intent(in out) :: this
        real(8), dimension(:,:,:), pointer, intent(out) :: adv_sie_ptr

!        call this%adv_sie%Point_to_data(adv_sie_ptr)

    end subroutine Point_to_adv_sie_data

    subroutine Point_to_momentum_data (this, momentum_i, momentum_j)
        class (advect_t), intent(in out) :: this
        real(8), dimension(:,:,:), pointer, intent(out) :: momentum_i, momentum_j

!        call this%momentum%Point_to_data(momentum_i, momentum_j)

    end subroutine Point_to_momentum_data

    subroutine Point_to_adv_momentum_data (this, adv_momentum_i, adv_momentum_j)
        class (advect_t), intent(in out) :: this
        real(8), dimension(:,:,:), pointer, intent(out) :: adv_momentum_i, adv_momentum_j

        call this%adv_momentum%Point_to_data(adv_momentum_i, adv_momentum_j)

    end subroutine Point_to_adv_momentum_data

    subroutine Point_to_adv_vertex_mass_data (this, vertex_mass_adv_ptr)
        class (advect_t), intent(in out) :: this
        real(8), dimension(:,:,:), pointer, intent(out) :: vertex_mass_adv_ptr

        call this%adv_vertex_mass%Point_to_data(vertex_mass_adv_ptr)

    end subroutine Point_to_adv_vertex_mass_data

    subroutine Set_advect_init_layer_mat (this, new_advect_init_layer_mat)
        class (advect_t), intent(in out) :: this
        integer, intent(in) :: new_advect_init_layer_mat

        this%advect_init_layer_mat = new_advect_init_layer_mat

    end subroutine Set_advect_init_layer_mat

    subroutine Get_all_areas (this, area_t_out, area_t_in, t_out, t_in, t_l_out, t_l_in,t_r_out, t_r_in, t_b_out, t_b_in, &
        t_t_out, t_t_in)
        class (advect_t), intent(in out) :: this
        real(8), dimension(:,:,:), pointer, intent(out)  :: area_t_out, area_t_in
        real(8), intent(out)                             :: t_out, t_in, t_l_out, t_l_in, t_r_out, t_r_in, t_b_out, t_b_in, t_t_out&
            , t_t_in

        t_out   = this%total_out
        t_in    = this%total_in
        t_l_out = this%total_left_out
        t_l_in  = this%total_left_in
        t_r_out = this%total_right_out
        t_r_in  = this%total_right_in
        t_b_out = this%total_bottom_out
        t_b_in  = this%total_bottom_in
        t_t_out = this%total_top_out
        t_t_in  = this%total_top_in

    end subroutine Get_all_areas

    subroutine Set_communication(this, comm, comm_params_cell, comm_params_vertex, comm_material)
        class (advect_t)            :: this
        type(communication_t), pointer            :: comm
        type(communication_parameters_t), pointer :: comm_params_cell, comm_params_vertex, comm_material
        integer :: i, num_mat


        call this%adv_cell_mass%Set_communication(comm, comm_params_cell)
        call this%adv_vertex_mass%Set_communication(comm, comm_params_vertex)
        call this%adv_sie%Set_communication(comm, comm_params_cell)
        call this%vel_adv_weight%Set_communication(comm, comm_params_vertex)

        call this%momentum%Set_communication(comm, comm_params_vertex)

            call this%adv_mats%Set_communication_material_advect(comm, comm_material)

    end subroutine Set_communication

    subroutine Clean_advect (this)
        implicit none
        class (advect_t), intent(in out) :: this

    end subroutine Clean_advect

    subroutine Write_advect(this, unit, iostat, iomsg)
        class (advect_t), intent(in) :: this
        integer,      intent(in)     :: unit
        integer,      intent(out)    :: iostat
        character(*), intent(in out)  :: iomsg

        write(*,*) '@@@ in Write_advect @@@'

        write(*,*) '@@@ end Write_data @@@'

    end subroutine Write_advect

    subroutine Read_advect(this, unit, iostat, iomsg)
        class (advect_t), intent(in out) :: this
        integer,      intent(in)     :: unit
        integer,      intent(out)    :: iostat
        character(*), intent(in out)  :: iomsg
        real(8), dimension(:,:,:), pointer :: tmp1

        write(*,*) "@@@ in Read_advect @@@"
        tmp1 => this%area_top_in%values

        write(*,*) "@@@ end Read_advect @@@"

    end subroutine Read_advect

end module advect_module
