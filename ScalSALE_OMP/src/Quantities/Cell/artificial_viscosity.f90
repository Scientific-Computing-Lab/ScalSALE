
module artificial_viscosity_module
    use cell_quantity_module, only : cell_quantity_t
    use cell_boundary_condition_module, only : cell_boundary_condition_t, cell_bc_wrapper_t
    use boundary_parameters_module      , only : boundary_parameters_t

    implicit none
    private
    public :: artificial_viscosity_t

    type, extends(cell_quantity_t) :: artificial_viscosity_t
        private

        real(8)                                         :: q_visc                          
        real(8)                                         :: xl_visc                         
        real(8), dimension(:,:,:), pointer, public      :: q_visc_mat                      
        real(8), dimension(:,:,:), pointer, public      :: xl_visc_mat                     
        integer, public                                 :: to_radial_index_sphere          
        integer, public                                 :: from_radial_index_sphere        
        integer                                         :: start_no_xl_visc
        integer                                         :: end_no_xl_visc
        logical                                         :: no_xl_visc_flag
        logical                                         :: update_visc_factor

    contains

        procedure, public :: Clean_artificial_viscosity
        procedure, public :: Get_quad_visc_factor
        procedure, public :: Get_linear_visc_factor

        procedure, public :: Update_visc_factors


        procedure, public :: Initialize

    end type


    interface artificial_viscosity_t

        module procedure Constructor
        module procedure Constructor1
    end interface artificial_viscosity_t

contains

    type(artificial_viscosity_t) function Constructor(initial_val, d1, d2, d3, bc, bc_params)
        implicit none
        real(8)                                       , intent(in) :: initial_val 
        integer                                       , intent(in) :: d1          
        integer                                       , intent(in) :: d2          
        integer                                       , intent(in) :: d3          
        type(cell_bc_wrapper_t), dimension(:), pointer, intent(in) :: bc           
        type(boundary_parameters_t), pointer, intent(in) :: bc_params

        Constructor%start_no_xl_visc = 0
        Constructor%end_no_xl_visc = 0

        call Constructor%Init_cell_quantity_init_val (initial_val, d1, d2, d3, bc, bc_params)
    end function

    type(artificial_viscosity_t) function Constructor1(initial_val, q_v, xl_v, d1, d2, d3, bc, bc_params)
        implicit none
        real(8)                                       , intent(in) :: initial_val 
        integer                                       , intent(in) :: d1          
        integer                                       , intent(in) :: d2          
        integer                                       , intent(in) :: d3          
        real(8)                                       , intent(in) :: q_v         
        real(8)                                       , intent(in) :: xl_v        
        type(cell_bc_wrapper_t), dimension(:), pointer, intent(in) :: bc          
        type(boundary_parameters_t), pointer, intent(inout) :: bc_params

        Constructor1%start_no_xl_visc = 0
        Constructor1%end_no_xl_visc = 0
        Constructor1%q_visc = q_v
        Constructor1%xl_visc = xl_v

        Constructor1%from_radial_index_sphere = 0
        Constructor1%to_radial_index_sphere   = -1

        call Constructor1%Init_cell_quantity_init_val (initial_val, d1, d2, d3, bc, bc_params)
    end function

    subroutine Clean_artificial_viscosity (this)
        class (artificial_viscosity_t), intent(in out) :: this 

        call this%Clean_cell_quantity ()
    end subroutine Clean_artificial_viscosity

    function Get_quad_visc_factor(this) result (q_v)
        class (artificial_viscosity_t), intent(in out) :: this
        real(8)                                        :: q_v

        q_v = this%q_visc

    end function Get_quad_visc_factor

    function Get_linear_visc_factor(this) result (xl_v)
        class (artificial_viscosity_t), intent(in out) :: this
        real(8)                                        :: xl_v

        xl_v = this%xl_visc

    end function Get_linear_visc_factor


    subroutine Update_visc_factors(this, total_rho, vof)
        use density_module, only: density_t
        use vof_module    , only: vof_t
        implicit none
        class (artificial_viscosity_t), intent(in out) :: this
        class (density_t), optional   , intent(in out) :: total_rho
        class (vof_t)    ,optional    , intent(in out) :: vof

        integer :: i,j,k, nx, ny, nz
        nx = this%d1
        ny = this%d2
        if (this%update_visc_factor .eqv. .false.) return
        if (.not. present(total_rho)) then
            do k = this%start_no_xl_visc, this%end_no_xl_visc
                do j = 1, ny
                    do i = 1, nx
                        this%xl_visc_mat(i,j,k) = 0d0
                    end do
                end do
            end do
            this%update_visc_factor = .false.
        else
        end if

    end subroutine Update_visc_factors

    subroutine Initialize(this, q_visc, xl_visc, i_sphere, i_sphere_up, start_index&
    ,no_xl_flag, start_no_xl_visc_layer, end_no_xl_visc_layer, offset_no_xl_visc)
        use parallel_parameters_module      , only : parallel_parameters_t
        implicit none
        class (artificial_viscosity_t), intent(in out) :: this
        real(8)                                       , intent(in) :: q_visc               
        real(8)                                       , intent(in) :: xl_visc              
        integer                                       , intent(in) :: i_sphere             
        integer                                       , intent(in) :: i_sphere_up          
        integer                                       , intent(in) :: start_no_xl_visc_layer     
        integer                                       , intent(in) :: end_no_xl_visc_layer       
        integer                                       , intent(in) :: offset_no_xl_visc    
        integer, dimension(:,:), allocatable          , intent(in) :: start_index
        logical                                       , intent(in) :: no_xl_flag
        integer :: k_start, k_end, inside
        integer :: k, i_sphere2



        this%xl_visc = xl_visc
        this%q_visc = q_visc

        inside = this%parallel_params%Inside_physical_world(i_sphere, axis=3)
        select case (inside)
            case (1)
                this%to_radial_index_sphere = this%d3 + 2
            case (0)
                this%to_radial_index_sphere = this%parallel_params%get_physical_index(i_sphere, 3)
            case (-1)
                this%to_radial_index_sphere = -1
        end select

        i_sphere2 = this%parallel_params%virt_nzp + 1 - i_sphere_up
        inside = this%parallel_params%Inside_physical_world(i_sphere2, axis=3)
        select case (inside)
            case (1)
                this%from_radial_index_sphere = huge(1)
            case (0)
                this%from_radial_index_sphere = this%parallel_params%get_physical_index(i_sphere2, 3)
            case (-1)
                this%from_radial_index_sphere = 1
        end select

        if (no_xl_flag .eqv. .false.) then
            this%no_xl_visc_flag = .false.
            this%start_no_xl_visc = 0
            this%end_no_xl_visc = 0
        else
            this%start_no_xl_visc = start_index(start_no_xl_visc_layer, 1) - offset_no_xl_visc
            this%end_no_xl_visc   = start_index(end_no_xl_visc_layer, 1)   + offset_no_xl_visc

            if (this%start_no_xl_visc == this%end_no_xl_visc) then
                this%no_xl_visc_flag = .false.
            else
                inside = this%parallel_params%Inside_physical_world(this%start_no_xl_visc, axis=3)
                select case(inside)
                    case(1)
                        this%start_no_xl_visc = 0
                        this%no_xl_visc_flag = .false.
                    case(0)
                        this%no_xl_visc_flag = .true.
                    case(-1)
                        this%no_xl_visc_flag  = .true.
                        this%start_no_xl_visc = 1
                end select

                inside = this%parallel_params%Inside_physical_world(this%end_no_xl_visc, axis=3)
                select case(inside)
                    case(1)
                        this%no_xl_visc_flag = .true.
                        this%end_no_xl_visc  = this%d3 + 1
                    case(0)
                        this%no_xl_visc_flag = .true.
                    case(-1)
                        this%end_no_xl_visc = 0
                        this%no_xl_visc_flag = .false.
                end select
            end if
        end if

        this%update_visc_factor = .false.
        if (this%no_xl_visc_flag .eqv. .true.) then
            this%update_visc_factor = .true.
        end if

        allocate(this%xl_visc_mat(1:this%d1, 1:this%d2, 1:this%d3), source=this%xl_visc)

    end subroutine Initialize

end module artificial_viscosity_module
