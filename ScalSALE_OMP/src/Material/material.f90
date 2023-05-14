
module material_module
    use data_module                   , only : data_t
    use equation_of_state_module      , only : equation_of_state_t, eos_wrapper_t
    use energy_module                 , only : energy_t
    use pressure_module               , only : pressure_t
    use cell_mass_module              , only : cell_mass_t
    use materials_in_cells_module     , only : materials_in_cells_t
    use vof_module                    , only : vof_t
    use density_module                , only : density_t
    use temperature_module            , only : temperature_t
    use sound_velocity_module         , only : sound_velocity_t
    use ideal_gas_module                , only : ideal_gas_t

    use cell_boundary_condition_module, only : cell_bc_wrapper_t
    use volume_module                 , only : volume_t
    use material_quantity_module      , only : material_quantity_t
    use material_base_module          , only : material_base_t
    use communication_module, only : communication_t
    use communication_parameters_module, only : communication_parameters_t

    implicit none
    private

    type, public, extends(material_base_t) :: material_t
        private
        real(8), dimension(:),public, pointer           :: atomic_mass
        real(8), dimension(:),public, pointer                        :: num_protons
        real(8), dimension(:), public,pointer                        :: num_protons_2
        real(8), dimension(:), public,pointer                        :: gamma_gas
        type (material_quantity_t),pointer, public       :: pressure
        type (material_quantity_t),pointer, public :: sound_vel
        type (material_quantity_t),pointer, public    :: temperature
        type (material_quantity_t),pointer, public    :: temperature_old
        type (material_quantity_t),pointer, public        :: density
        type (material_quantity_t),pointer, public           :: dp_de
        type (material_quantity_t),pointer, public           :: dp_drho
        type (material_quantity_t),pointer, public           :: dt_de
        type (material_quantity_t),pointer, public           :: dt_drho
        type (eos_wrapper_t), dimension(:), pointer           :: equation_of_state

    contains


        procedure, public :: Apply_eos

        procedure, public :: Clean_data

        procedure, public :: Set_communication_material
      
      
        procedure, public :: Write_material
        procedure, public :: Write_material_abstract => Write_material

        procedure, public :: Read_material
        procedure, public :: Read_material_abstract => Read_material
    end type material_t

    interface material_t
        module procedure Constructor
    end interface material_t

contains

    type(material_t) function Constructor(nxp, nyp, nzp, nmats, mat_ids, gamma_gas, atomic_mass,&
        num_protons, num_protons_2, rho_0, temperature_init, sie_0, mat_cells, bc_cell&
        , bc_params)
        use boundary_parameters_module, only : boundary_parameters_t

        implicit none
        !        type (eos_wrapper_t), dimension(:), pointer                  , intent(in out)       :: eos
        integer, dimension(:), allocatable        , intent(in)           :: mat_ids
        integer                               , intent(in)           :: nxp
        integer                               , intent(in)           :: nyp
        integer                               , intent(in)           :: nzp
        integer                               , intent(in)           :: nmats
        real(8), dimension(:), allocatable                               , intent(in)           :: atomic_mass
        real(8) , dimension(:), allocatable                              , intent(in)           :: gamma_gas
        real(8), dimension(:), allocatable                               , intent(in)           :: rho_0
        real(8), dimension(:), allocatable                               , intent(in)           :: sie_0
        real(8)                               , intent(in)           :: temperature_init
        real(8), dimension(:), allocatable                               , intent(in)           :: num_protons
        real(8), dimension(:), allocatable                               , intent(in)           :: num_protons_2
        type(materials_in_cells_t)   ,pointer         , intent(in out)       :: mat_cells
        type(cell_bc_wrapper_t), dimension(:), pointer,  intent(inout) :: bc_cell
        type(boundary_parameters_t), pointer, intent(in) :: bc_params

        !        type(materials_in_cells_t)   ,pointer                :: mat_cells
        !        type(cell_bc_wrapper_t), dimension(:), pointer:: bc_cell
        !        type(boundary_parameters_t), pointer :: bc_params
        !        real(8), dimension(:), pointer                                          :: atomic_mass
        !        real(8) , dimension(:), pointer                                         :: gamma_gas
        !        real(8), dimension(:), pointer                                         :: rho_0
        !        real(8), dimension(:), pointer                                          :: sie_0
        !        real(8)                                          :: temperature_init
        !        real(8), dimension(:), pointer                                          :: num_protons
        !        real(8), dimension(:), pointer                                          :: num_protons_2

        real(8), dimension (:,:,:,:), pointer                          :: density_vof
        real(8), dimension (:,:,:,:), pointer                          :: sie_vof
        real(8), dimension (:,:,:,:), pointer                          :: mat_vof
        real(8), dimension (:,:,:,:), pointer                          :: temp, temp_old
        real(8), dimension (:,:,:), pointer                          :: mat_cell
        real(8), dimension (:,:,:,:), pointer                          :: delete
        real(8), dimension (:,:,:,:), pointer                          :: delete2
        type(eos_wrapper_t), allocatable                                :: eos_c_wrap
        type(ideal_gas_t), target                                       :: ig_eos_c


        integer                                                      :: i, j, k, m

        allocate(Constructor%dp_de)
        allocate(Constructor%dp_drho)
        allocate(Constructor%dt_de)
        allocate(Constructor%dt_drho)
        allocate(Constructor%pressure)
        allocate(Constructor%temperature)
        allocate(Constructor%temperature_old)
        allocate(Constructor%sound_vel)
        allocate(Constructor%density)
        allocate(Constructor%nrg_calc(nmats))
        allocate(Constructor%num_protons_2(nmats))
        allocate(Constructor%num_protons(nmats))
        allocate(Constructor%gamma_gas(nmats))
        allocate(Constructor%atomic_mass(nmats))
        allocate(eos_wrapper_t :: Constructor%equation_of_state (nmats))


        call Constructor%Init_material_base(nxp, nyp, nzp, nmats, mat_ids, bc_cell, bc_params)

        Constructor%dp_de   = material_quantity_t(0d0, nxp, nyp, nzp, nmats)
        Constructor%dp_drho = material_quantity_t(0d0, nxp, nyp, nzp, nmats)
        Constructor%dt_de   = material_quantity_t(0d0, nxp, nyp, nzp, nmats)
        Constructor%dt_drho = material_quantity_t(0d0, nxp, nyp, nzp, nmats)
        Constructor%density = material_quantity_t(0d0, nxp, nyp, nzp, nmats, bc_cell, bc_params)

        Constructor%pressure = material_quantity_t (0d0, nxp, nyp, nzp, nmats, bc_cell, bc_params)
        Constructor%temperature = material_quantity_t (0d0, nxp, nyp, nzp, nmats, bc_cell, bc_params)
        Constructor%temperature_old = material_quantity_t (0d0, nxp, nyp, nzp, nmats, bc_cell, bc_params)
        Constructor%sound_vel = material_quantity_t (0d0, nxp, nyp, nzp, nmats, bc_cell, bc_params)

        allocate(eos_c_wrap)
        eos_c_wrap%eos => ig_eos_c

        do i = 1, nmats
            Constructor%equation_of_state(i) = eos_c_wrap
            Constructor%atomic_mass(i) = atomic_mass(i)
            Constructor%num_protons(i) = num_protons(i)
            Constructor%num_protons_2(i) = num_protons_2(i)
            Constructor%gamma_gas(i) = gamma_gas(i)
        end do



        call Constructor%density%Point_to_data(density_vof)
        call Constructor%vof%Point_to_data(mat_vof)
        call Constructor%temperature%Point_to_data(temp)
        call Constructor%temperature_old%Point_to_data(temp_old)
        call Constructor%sie%Point_to_data(sie_vof)
        call mat_cells%Point_to_data(mat_cell)

        do k = 1, nzp
            do j = 1, nyp
                do i = 1, nxp
                    do m = 1, nmats
                        if (mat_cell(i, j, k) == mat_ids(m)) then
                            density_vof(m, i, j, k) = rho_0(m)
                            temp(m, i, j, k)        = temperature_init
                            temp_old(m, i, j, k)    = temperature_init
                            mat_vof(m, i, j, k)     = 1d0
                            sie_vof(m,i,j,k) = sie_0(m)
                            if (sie_0(m) == 0) then
                                Constructor%nrg_calc(m) = 1
                            else
                                Constructor%nrg_calc(m) = 0
                            end if
                        end if
                    end do
                end do
            end do
        end do
    end function

    subroutine Apply_eos(this, nx, ny, nz, emf, is_old_temperature)
        class(material_t)                  , intent(in out)    :: this

        !        integer                            , intent(in    ) :: mat_num

        !        integer                            , intent(in    ) :: nrg_or_tmp
        integer                            , intent(in    ) :: nx
        integer                            , intent(in    ) :: ny
        integer                            , intent(in    ) :: nz
        real(8)                            , intent(in    ) :: emf
        logical                            , intent(in    ) :: is_old_temperature

        real(8), dimension (:,:,:,:), pointer :: rho
        real(8), dimension (:,:,:,:), pointer :: p
        real(8), dimension (:,:,:,:), pointer :: t
        real(8), dimension (:,:,:,:), pointer :: e
        real(8), dimension (:,:,:,:), pointer :: dp_de_p
        real(8), dimension (:,:,:,:), pointer :: dp_drho_p
        real(8), dimension (:,:,:,:), pointer :: dt_de_p
        real(8), dimension (:,:,:,:), pointer :: dt_drho_p
        real(8), dimension (:,:,:,:), pointer :: mat_vof
        real(8), dimension (:,:,:,:), pointer :: sound_vel

        integer :: m
        call this%pressure%Point_to_data(p)
        call this%density%Point_to_data(rho)
        call this%sie%Point_to_data(e)
        call this%dp_de%Point_to_data(dp_de_p)
        call this%dp_drho%Point_to_data(dp_drho_p)
        call this%dt_de%Point_to_data(dt_de_p)
        call this%dt_drho%Point_to_data(dt_drho_p)
        call this%vof%Point_to_data(mat_vof)
        call this%sound_vel%Point_to_data(sound_vel)
        if (is_old_temperature .eqv. .true.) then
            call this%temperature_old%Point_to_data(t)
        else
            call this%temperature%Point_to_data(t)
        end if
        do m = 1, this%nmats
            call this%equation_of_state(m)%eos%Calculate(p, sound_vel, rho, e, &
                dp_de_p, dp_drho_p, dt_de_p, dt_drho_p, this%gamma_gas(m), this%atomic_mass(m), &
                t, m, this%nrg_calc(m),&
                nx, ny, nz, mat_vof, emf)
        end do
        this%nrg_calc = 0

    end subroutine Apply_eos


    subroutine Set_communication_material(this, comm, comm_params)
        class (material_t)            :: this
        type(communication_t), pointer            :: comm
        type(communication_parameters_t), pointer :: comm_params

        call this%Set_communication_material_base(comm, comm_params)

        call this%dp_de%Set_communication(comm, comm_params)
        call this%dp_drho%Set_communication(comm, comm_params)
        call this%dt_de%Set_communication(comm, comm_params)
        call this%dt_drho%Set_communication(comm, comm_params)
        call this%pressure%Set_communication(comm, comm_params)
        call this%temperature%Set_communication(comm, comm_params)
        call this%temperature_old%Set_communication(comm, comm_params)
        call this%sound_vel%Set_communication(comm, comm_params)
        call this%density%Set_communication(comm, comm_params)
    end subroutine Set_communication_material


    subroutine Clean_data(this)
        class (material_t)            :: this

    !        call this%Clean_material_base()
    !        call this%dp_de          %Clean_data
    !        call this%dp_drho        %Clean_data
    !        call this%dt_de          %Clean_data
    !        call this%dt_drho        %Clean_data
    !        call this%pressure       %Clean_pressure
    !        call this%density        %Clean_density
    !        call this%temperature    %Clean_temperature
    !        call this%temperature_old%Clean_temperature
    end subroutine Clean_data

    subroutine Write_material(this, unit, iostat, iomsg)
        class (material_t), intent(in) :: this
        integer,      intent(in)     :: unit
        integer,      intent(out)    :: iostat
        character(*), intent(in out)  :: iomsg
    !
    !#ifdef DEBUG
    !        write(*,*) '@@@ in Write_material @@@'
    !#endif
    !
    !        call this%Write_material_base(unit, iostat=iostat, iomsg=iomsg)
    !
    !        call this%pressure%Write_quantity_abstract(unit, iostat=iostat, iomsg=iomsg)
    !        call this%sound_vel%Write_quantity_abstract(unit, iostat=iostat, iomsg=iomsg)
    !        call this%temperature%Write_quantity_abstract(unit, iostat=iostat, iomsg=iomsg)
    !        call this%temperature_old%Write_quantity_abstract(unit, iostat=iostat, iomsg=iomsg)
    !        call this%density%Write_quantity_abstract(unit, iostat=iostat, iomsg=iomsg)
    !
    !        write(unit, iostat=iostat, iomsg=iomsg) &
    !            this%atomic_mass, &
    !            this%num_protons, &
    !            this%num_protons_2, &
    !            this%gamma_gas, &
    !            this%dp_de, &
    !            this%dp_drho, &
    !            this%dt_de, &
    !            this%dt_drho
    !
    !#ifdef DEBUG
    !        write(*,*) &
    !            'atomic_mass', &
    !            this%atomic_mass, &
    !            'num_protons', &
    !            this%num_protons, &
    !            'num_protons_2', &
    !            this%num_protons_2, &
    !            'gamma_gas', &
    !            this%gamma_gas, &
    !            '###'
    !
    !        write(*,*) '@@@ end Write_material @@@'
    !#endif

    end subroutine Write_material

    subroutine Read_material(this, unit, iostat, iomsg)
        class (material_t), intent(in out) :: this
        integer,      intent(in)     :: unit
        integer,      intent(out)    :: iostat
        character(*), intent(in out)  :: iomsg
    !
    !#ifdef DEBUG
    !        write(*,*) '@@@ in Read_material @@@'
    !#endif
    !
    !        call this%Read_material_base(unit, iostat=iostat, iomsg=iomsg)
    !
    !        call this%pressure%Read_quantity_abstract(unit, iostat=iostat, iomsg=iomsg)
    !        call this%sound_vel%Read_quantity_abstract(unit, iostat=iostat, iomsg=iomsg)
    !        call this%temperature%Read_quantity_abstract(unit, iostat=iostat, iomsg=iomsg)
    !        call this%temperature_old%Read_quantity_abstract(unit, iostat=iostat, iomsg=iomsg)
    !        call this%density%Read_quantity_abstract(unit, iostat=iostat, iomsg=iomsg)
    !
    !        read(unit, iostat=iostat, iomsg=iomsg) &
    !            this%atomic_mass, &
    !            this%num_protons, &
    !            this%num_protons_2, &
    !            this%gamma_gas, &
    !            this%dp_de, &
    !            this%dp_drho, &
    !            this%dt_de, &
    !            this%dt_drho
    !
    !#ifdef DEBUG
    !        write(*,*) &
    !            'atomic_mass', &
    !            this%atomic_mass, &
    !            'num_protons', &
    !            this%num_protons, &
    !            'num_protons_2', &
    !            this%num_protons_2, &
    !            'gamma_gas', &
    !            this%gamma_gas, &
    !            '###'
    !
    !        write(*,*) '@@@ end Read_material @@@'
    !#endif

    end subroutine Read_material

end module material_module

