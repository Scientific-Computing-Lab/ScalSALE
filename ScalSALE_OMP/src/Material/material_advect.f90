
module material_advect_module

    use energy_module                 , only : energy_t
    use cell_mass_module              , only : cell_mass_t
    use vof_module                    , only : vof_t
    use cell_boundary_condition_module, only : cell_bc_wrapper_t
    use data_module                   , only : data_t
    use material_base_module          , only : material_base_t
    use communication_module, only : communication_t
    use communication_parameters_module, only : communication_parameters_t
    use material_quantity_module, only : material_quantity_t

    implicit none
    private

    type, public, extends(material_base_t) :: material_advect_t
        private

        type(material_quantity_t), pointer, public        :: area_top_in
        type(material_quantity_t), pointer, public         :: area_top_out
        real(8), dimension(:), pointer, public      :: area_left_in
        real(8), dimension(:), pointer, public      :: area_left_out
        real(8), dimension(:), pointer, public      :: area_right_in
        real(8), dimension(:), pointer, public      :: area_right_out
        real(8), dimension(:), pointer, public      :: area_bottom_in
        real(8), dimension(:), pointer, public      :: area_bottom_out
        real(8), dimension(:), pointer, public      :: total_out
        real(8), dimension(:), pointer, public      :: total_in
        real(8), dimension(:), pointer, public      :: fxtm
        real(8), dimension(:), pointer, public      :: sxt1
        real(8), dimension(:), pointer, public      :: sxt2
        integer, dimension(:), pointer, public      :: sgn
        type(material_quantity_t),pointer, public      :: a
        type(material_quantity_t),pointer, public      :: b
        type(material_quantity_t),pointer, public      :: c
        type(material_quantity_t),pointer, public      :: side



    contains

        procedure, public :: Set_communication_material_advect

        procedure, public :: Clean_data

        procedure, public :: Point_to_areas

        procedure, public :: Point_to_area_top_out

        procedure, public :: Point_to_area_top_in

    end type material_advect_t

    interface material_advect_t
        module procedure Constructor
    end interface material_advect_t


contains

    type(material_advect_t) function Constructor(nxp, nyp, nzp, nmats, mat_ids, bc_cell, bc_params)
        use boundary_parameters_module, only : boundary_parameters_t
        implicit none
        integer                               , intent(in)           :: nmats
        integer,dimension(:), allocatable         , intent(in)           :: mat_ids

        integer                               , intent(in)           :: nxp
        integer                               , intent(in)           :: nyp
        integer                               , intent(in)           :: nzp

        type(cell_bc_wrapper_t), dimension(:), pointer,  intent(in) :: bc_cell
        type(boundary_parameters_t), pointer, intent(in) :: bc_params

        call Constructor%Init_material_base(nxp, nyp, nzp, nmats, mat_ids, bc_cell, bc_params)

        allocate(Constructor%area_left_in(nmats))
        allocate(Constructor%area_left_out(nmats))
        allocate(Constructor%area_right_in(nmats))
        allocate(Constructor%area_right_out(nmats))
        allocate(Constructor%area_bottom_in(nmats))
        allocate(Constructor%area_bottom_out(nmats))
        allocate(Constructor%total_out(nmats))
        allocate(Constructor%total_in(nmats))
        allocate(Constructor%fxtm(nmats))
        allocate(Constructor%sxt1(nmats))
        allocate(Constructor%sxt2(nmats))
        allocate(Constructor%sgn(nmats))

        allocate(Constructor%a)
        allocate(Constructor%b)
        allocate(Constructor%c)
        allocate(Constructor%side)
        allocate(Constructor%area_top_out)
        allocate(Constructor%area_top_in)
        Constructor%area_top_in = material_quantity_t (0d0, nxp, 1, 1, nmats)
        Constructor%area_top_out = material_quantity_t (0d0, nxp, 1, 1, nmats)
        Constructor%a = material_quantity_t(0d0, nxp, nyp, nzp, nmats)
        Constructor%b = material_quantity_t(0d0, nxp, nyp, nzp, nmats)
        Constructor%c = material_quantity_t(0d0, nxp, nyp, nzp, nmats)
        Constructor%side = material_quantity_t(0d0, nxp, nyp, nzp, nmats)
        Constructor%area_left_in = 0
        Constructor%area_left_out = 0
        Constructor%area_right_in = 0
        Constructor%area_right_out = 0
        Constructor%area_bottom_in = 0
        Constructor%area_bottom_out = 0

    end function

    subroutine Point_to_areas (this, ptr_top_in, ptr_top_out)
        implicit none
        class (material_advect_t)                , intent(in out) :: this
        real(8)       , dimension(:,:,:), pointer, intent(out)    :: ptr_top_in
        real(8)       , dimension(:,:,:), pointer, intent(out)    :: ptr_top_out

!        call this%area_top_in %Point_to_data (ptr_top_in)
!        call this%area_top_out%Point_to_data (ptr_top_out)

    end subroutine Point_to_areas


    subroutine Point_to_area_top_in (this, ptr_top_in)
        implicit none
        class (material_advect_t)                , intent(in out) :: this
        real(8)       , dimension(:,:,:), pointer, intent(out)    :: ptr_top_in

!        call this%area_top_in %Point_to_data (ptr_top_in)
    end subroutine Point_to_area_top_in

    subroutine Point_to_area_top_out (this, ptr_top_out)
        implicit none
        class (material_advect_t)                , intent(in out) :: this
        real(8)       , dimension(:,:,:), pointer, intent(out)    :: ptr_top_out
!
!        call this%area_top_out%Point_to_data (ptr_top_out)
    end subroutine Point_to_area_top_out


    subroutine Set_communication_material_advect(this, comm, comm_params)
        class (material_advect_t)            :: this
        type(communication_t), pointer            :: comm
        type(communication_parameters_t), pointer :: comm_params

        call this%Set_communication_material_base(comm, comm_params)

        call this%a%Set_communication(comm, comm_params)
        call this%b%Set_communication(comm, comm_params)
        call this%c%Set_communication(comm, comm_params)
        call this%side%Set_communication(comm, comm_params)

    end subroutine Set_communication_material_advect



    subroutine Clean_data(this)
        class (material_advect_t), intent(in out) :: this

!        call this%Clean_material_base()
!
!        call this%area_top_in%Clean_data ()
!        call this%area_top_out%Clean_data ()

    end subroutine Clean_data

end module material_advect_module

