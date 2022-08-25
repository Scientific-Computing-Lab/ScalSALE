
module hydro_step_module
    use vertex_boundary_condition_module, only : vertex_boundary_condition_t
    use cell_boundary_condition_module  , only : cell_boundary_condition_t
    use num_materials_in_cells_module   , only : num_materials_in_cells_t
    use artificial_viscosity_module     , only : artificial_viscosity_t
    use materials_in_cells_module       , only : materials_in_cells_t
    use equation_of_state_module        , only : equation_of_state_t
    use vertex_quantity_module          , only : vertex_quantity_t
    use sound_velocity_module           , only : sound_velocity_t
    use cell_quantity_module            , only : cell_quantity_t
    use acceleration_module             , only : acceleration_t
    use free_surface_module             , only : free_surface_t
    use temperature_module              , only : temperature_t
    use vertex_mass_module              , only : vertex_mass_t
    use ideal_gas_module                , only : ideal_gas_t
    use cell_mass_module                , only : cell_mass_t
    use mesh_base_module                , only : mesh_base_t
    use pressure_module                 , only : pressure_t
    use velocity_module                 , only : velocity_t
    use material_module                 , only : material_t
    use datafile_module                 , only : datafile_t
    use density_module                  , only : density_t
    use energy_module                   , only : energy_t
    use volume_module                   , only : volume_t
    use rezone_module                   , only : rezone_t
    use advect_module                   , only : advect_t
    use data_module                     , only : data_t
    use time_module                     , only : time_t
    use vof_module                      , only : vof_t
    use parallel_parameters_module      , only : parallel_parameters_t
    use communication_parameters_module , only : communication_parameters_t
    use communication_module            , only : communication_t
    use boundary_parameters_module      , only : boundary_parameters_t

    implicit none
    private

    type, public :: hydro_step_t
        private
        type (num_materials_in_cells_t), pointer :: num_mat_cells        
        type (artificial_viscosity_t)  , pointer :: a_visc               
        type (materials_in_cells_t)    , pointer :: mat_id               
        type (sound_velocity_t)        , pointer :: total_sound_vel      
        type (acceleration_t)          , pointer :: acceleration         
        type (temperature_t)           , pointer :: total_temperature    
        type (vertex_mass_t)           , pointer :: vertex_mass    
        type (vertex_mass_t)           , pointer :: previous_vertex_mass 
        class (mesh_base_t)            , pointer :: mesh                 
        type (cell_mass_t)             , pointer :: total_cell_mass      
        type (cell_mass_t)             , pointer :: previous_cell_mass   
        type (pressure_t)              , pointer :: total_pressure       
        type (pressure_t)              , pointer :: total_pressure_sum   
        type (velocity_t)              , pointer :: velocity             
        type (density_t)               , pointer :: total_density        
        type (volume_t)                , pointer :: total_volume         
        type (energy_t)                , pointer :: total_sie            
        type (rezone_t)                , pointer :: rezone               
        type (advect_t)                , pointer :: advect               
        type (data_t)                  , pointer :: inversed_vertex_mass
        type (vof_t)                   , pointer :: total_vof            
        type (data_t)                  , pointer :: total_dt_drho        
        type (data_t)                  , pointer :: total_dp_drho        
        type (data_t)                  , pointer :: total_dt_de          
        type (data_t)                  , pointer :: total_dp_de          
        type (data_t)                  , pointer :: dvelocity_x_dx       
        type (data_t)                  , pointer :: dvelocity_y_dx       
        type (data_t)                  , pointer :: dvelocity_z_dx       
        type (data_t)                  , pointer :: dvelocity_x_dy       
        type (data_t)                  , pointer :: dvelocity_y_dy       
        type (data_t)                  , pointer :: dvelocity_z_dY       
        type (data_t)                  , pointer :: dvelocity_x_dz       
        type (data_t)                  , pointer :: dvelocity_y_dz       
        type (data_t)                  , pointer :: dvelocity_z_dz       
        type (parallel_parameters_t)   , pointer :: parallel_params      
        type (material_t),  pointer :: materials

        integer :: wilkins_scheme    
        integer, public                           :: nmats
        integer, public                           :: nxp, nyp, nzp, nz, nx, ny
        real(8), public                           :: emf, emfm 
        real(8) :: init_temperature


        integer, public                           :: ncyc 
        integer, public                           :: cyc_delete
    contains

        procedure, public  :: do_time_step_2d

        procedure, public  :: do_time_step_3d

        procedure, public  :: debug_check_nan

        procedure, private :: Calculate_thermodynamics

        procedure, private :: Calculate_density

        procedure, private :: Calculate_inversed_vertex_mass

        procedure, private :: Calculate_acceleration_2d

        procedure, private :: Calculate_acceleration_3d

        procedure, private :: Calculate_artificial_viscosity_2d

        procedure, private :: Calculate_artificial_viscosity_3d

        procedure, private :: Fix_vof_2d

        procedure, private :: Fix_vof_3d

        procedure, private :: Calculate_velocity_2d

        procedure, private :: Calculate_velocity_3d

        procedure, private :: Calculate_energy_2d

        procedure, private :: Calculate_energy_3d

        procedure, private :: Calculate_mesh_2d

        procedure, private :: Calculate_mesh_3d

        procedure, private :: Boundary_mirror_image

        procedure, private :: Calculate_stresd

        procedure, private :: Apply_wilkins_mass

        procedure, private :: Restore_wilkins_mass

        procedure, public :: TEST_Fix_vof_2d => Fix_vof_2d

        procedure, public :: TEST_Restore_wilkins_mass => Restore_wilkins_mass

        procedure, public :: TEST_Apply_wilkins_mass => Apply_wilkins_mass

        procedure, public :: TEST_Calculate_thermodynamics => Calculate_thermodynamics

        procedure, public :: TEST_Calculate_inversed_vertex_mass => Calculate_inversed_vertex_mass

        procedure, public :: TEST_Calculate_density => Calculate_density

        procedure, public :: TEST_Calculate_acceleration_test => Calculate_acceleration_2d

        procedure, public :: TEST_Calculate_artificial_viscosity => Calculate_artificial_viscosity_2d

        procedure, public :: TEST_Calculate_velocity => Calculate_velocity_2d

        procedure, public :: TEST_Calculate_energy => Calculate_energy_2d

        procedure, public :: TEST_Calculate_mesh => Calculate_mesh_2d

        procedure, public :: TEST_Calculate_stresd => Calculate_stresd

        procedure, public :: TEST_Wilkins_mass => Apply_wilkins_mass

        procedure, public :: Set_communication

        procedure, public :: Point_to_velocity_data_2d

        procedure, public :: Point_to_velocity_data_3d

        generic, public :: Point_to_velocity_data =>  &
            Point_to_velocity_data_2d, &
            Point_to_velocity_data_3d

        procedure, public :: Point_to_acceleration_data_2d

        procedure, public :: Point_to_r_factor_data

        procedure, public :: Point_to_cell_mass_data

        procedure, public :: Point_to_vol_data

        procedure, public :: Point_to_artificial_viscosity_data

        procedure, public :: Point_to_density_data

        procedure, public :: Point_to_density_vof_data

        procedure, public :: Point_to_vof_data

        procedure, public :: Point_to_cell_mass_vof_data

        procedure, public :: Point_to_mat_vof_data

        procedure, public :: Point_to_temperature_data

        procedure, public :: Point_to_pressure_data

        procedure, public :: Point_to_sie_data

        procedure, public :: Point_to_sound_velocity_data

        procedure, public :: Point_to_temperature_vof_data

        procedure, public :: Point_to_temperature_old_vof_data

        procedure, public :: Point_to_pressure_vof_data

        procedure, public :: Point_to_sie_vof_data

        procedure, public :: Point_to_sound_velocity_vof_data

        procedure, public :: Point_to_deriv_vof_data

        procedure, public :: Point_to_mesh_velocity_data

        procedure, public :: Point_to_num_mat_cells_data

        procedure, public :: Point_to_mat_id_data

        procedure, public :: Point_to_mesh_data_2d

        procedure, public :: Point_to_mesh_data_3d

        generic,   public    :: Point_to_mesh_data  =>         &
            Point_to_mesh_data_2d,         &
            Point_to_mesh_data_3d

        procedure, public :: Point_to_vertex_mass_data

        procedure, public :: Point_to_inversed_vertex_mass_data

        procedure, public :: Point_to_advect

        procedure, public :: Point_to_rezone
      
      
        procedure, public :: Write_hydro_step
        generic :: write(unformatted) => Write_hydro_step

        procedure, public :: Read_hydro_step
        generic :: read(unformatted) => Read_hydro_step
      
    end type hydro_step_t

    interface hydro_step_t
        module procedure Constructor
    end interface hydro_step_t
  
contains

    type(hydro_step_t) function Constructor (df, nx, ny, nz, nxp, nyp, nzp, wilkins_scheme, mesh, velocity, acceleration,&
        total_volume, total_vof, total_sie, total_pressure, total_pressure_sum, total_density, &
        total_temperature, total_cell_mass, previous_cell_mass, vertex_mass, &
        previous_vertex_mass, inversed_vertex_mass, total_sound_vel,&
        a_visc, total_dp_de, total_dp_drho, total_dt_de, total_dt_drho, init_temperature, &
        nmats, materials, num_mat_cells, mat_id, emf, emfm, parallel_params, mat_ids)
        implicit none
        type(datafile_t), intent(in) :: df 

        integer, intent(in) :: nx   
        integer, intent(in) :: ny   
        integer, intent(in) :: nz   
        integer, intent(in) :: nxp  
        integer, intent(in) :: nyp  
        integer, intent(in) :: nzp  

        integer, intent(in) :: wilkins_scheme

        integer, intent(in) :: nmats
        real(8), intent(in) :: init_temperature 
        real(8), intent(in) :: emf              
        real(8), intent(in) :: emfm             

        class (mesh_base_t)          , pointer, intent(in out) :: mesh             

        type (artificial_viscosity_t), pointer, intent(in) :: a_visc               
        type (sound_velocity_t)      , pointer, intent(in) :: total_sound_vel      
        type (acceleration_t)        , pointer, intent(in) :: acceleration         
        type (temperature_t)         , pointer, intent(in) :: total_temperature    
        type (vertex_mass_t)         , pointer, intent(in) :: vertex_mass    
        type (vertex_mass_t)         , pointer, intent(in) :: previous_vertex_mass 
        type (cell_mass_t)           , pointer, intent(inout) :: total_cell_mass      
        type (cell_mass_t)           , pointer, intent(in) :: previous_cell_mass   
        type (pressure_t)            , pointer, intent(in) :: total_pressure       
        type (pressure_t)            , pointer, intent(in) :: total_pressure_sum   
        type (velocity_t)            , pointer, intent(inout) :: velocity             
        type (density_t)             , pointer, intent(in) :: total_density        
        type (energy_t)              , pointer, intent(in) :: total_sie            
        type (volume_t)              , pointer, intent(in) :: total_volume         
        type (vof_t)                 , pointer, intent(in) :: total_vof
        type (data_t)                , pointer, intent(in) :: inversed_vertex_mass 
        type (data_t)                , pointer, intent(in) :: total_dt_drho   
        type (data_t)                , pointer, intent(in) :: total_dp_drho   
        type (data_t)                , pointer, intent(in) :: total_dp_de     
        type (data_t)                , pointer, intent(in) :: total_dt_de     
        type (parallel_parameters_t) , pointer, intent(in) :: parallel_params     

        type (num_materials_in_cells_t), pointer, intent(in) :: num_mat_cells 
        type (materials_in_cells_t)    , pointer, intent(in) :: mat_id     

        type (material_t), pointer, intent(inout) :: materials
        integer,dimension(:), allocatable         , intent(in)           :: mat_ids


        Constructor%init_temperature     =  init_temperature
        Constructor%nx                   =  nx
        Constructor%ny                   =  ny
        Constructor%nz                   =  nz
        Constructor%nxp                  =  nxp
        Constructor%nyp                  =  nyp
        Constructor%nzp                  =  nzp
        Constructor%emf                  =  emf
        Constructor%emfm                 =  emfm
        Constructor%wilkins_scheme       =  wilkins_scheme
        Constructor%nmats          =  nmats
        Constructor%num_mat_cells        => num_mat_cells
        Constructor%mat_id               => mat_id
        Constructor%mesh                 => mesh
        Constructor%materials            => materials
        Constructor%acceleration         => acceleration
        Constructor%velocity             => velocity
        Constructor%vertex_mass    => vertex_mass
        Constructor%previous_vertex_mass => previous_vertex_mass
        Constructor%total_temperature    => total_temperature
        Constructor%total_vof            => total_vof
        Constructor%total_sie            => total_sie
        Constructor%total_pressure       => total_pressure
        Constructor%total_pressure_sum   => total_pressure_sum
        Constructor%total_cell_mass      => total_cell_mass
        Constructor%previous_cell_mass   => previous_cell_mass
        Constructor%total_sound_vel      => total_sound_vel
        Constructor%total_volume         => total_volume
        Constructor%a_visc               => a_visc
        Constructor%total_density        => total_density
        Constructor%inversed_vertex_mass => inversed_vertex_mass
        Constructor%total_dp_de          => total_dp_de
        Constructor%total_dp_drho        => total_dp_drho
        Constructor%total_dt_de          => total_dt_de
        Constructor%total_dt_drho        => total_dt_drho
        Constructor%cyc_delete = 0
        Constructor%parallel_params => parallel_params
        allocate(Constructor%rezone)

        Constructor%rezone = rezone_t(df%rezone_type, nxp, nyp, nzp, total_cell_mass%boundary_conditions,&
            velocity%boundary_conditions, mesh, velocity, vertex_mass, &
            total_cell_mass%boundary_params)
        allocate(Constructor%advect)



        Constructor%advect = advect_t(nxp, nyp, nzp,nmats, mat_ids, total_cell_mass%boundary_conditions, velocity%boundary_conditions, &
            total_cell_mass%boundary_params&
            ,df%line_calc, df%shorter_advect, df%fix_overflow, Constructor%rezone, mesh, &
            materials, mat_id, num_mat_cells, &
            total_sie, total_vof, &
            total_density, total_cell_mass, vertex_mass, &
            total_volume, velocity, emf, emfm , wilkins_scheme, parallel_params)
        Constructor%ncyc = 0
    end function


    subroutine do_time_step_2d(this, time)
        class (hydro_step_t)         , intent (in out) :: this 
        type (time_t)                , intent (in out) :: time 
        real(8), dimension(:, :, :), pointer :: reem
        integer :: i,j


        call this%Calculate_thermodynamics()

        call time%Calculate_dt(this%mesh, this%velocity, this%rezone%mesh_velocity, this%vertex_mass, this%total_vof, this%emfm)

        this%cyc_delete = this%cyc_delete+1

        if  (this%wilkins_scheme == 1) call this%Apply_wilkins_mass()

        call this%Calculate_acceleration_2d(this%wilkins_scheme, time%dt_mid)


        call this%Calculate_artificial_viscosity_2d(time, this%wilkins_scheme)

        call this%Calculate_velocity_2d(time%dt_mid, this%wilkins_scheme)

        if  (this%wilkins_scheme == 1) call this%Restore_wilkins_mass()

        call this%Calculate_energy_2d(this%init_temperature, time%dt)

        call this%Fix_vof_2d(time%dt)

        call this%rezone%Calculate_rezone_2d(time%dt)

        call this%Calculate_mesh_2d(time)

        if (this%rezone%rezone_type /= 0) call this%advect%Calculate_advect_2d()

        this%ncyc = this%ncyc  + 1
    end subroutine do_time_step_2d


    subroutine do_time_step_3d(this, time)
        use omp_lib
        class (hydro_step_t)         , intent (in out) :: this 
        type (time_t)                , intent (in out) :: time 

        real(8), dimension(:, :, :), pointer :: p,x              
        integer :: i,j,k
        call this%total_sie%Point_to_data(1, x)


        call this%Calculate_thermodynamics()

        call time%Calculate_dt(this%mesh, this%velocity, this%rezone%mesh_velocity, this%vertex_mass, this%total_vof, this%emfm)
        this%cyc_delete = this%cyc_delete+1

        call this%Calculate_acceleration_3d(time%dt_mid)

        call this%Calculate_artificial_viscosity_3d(time)

        call this%Calculate_velocity_3d(time%dt_mid)
        call this%rezone%Calculate_rezone_3d(time%dt)

        call this%Calculate_energy_3d(this%init_temperature, time%dt)

        call this%Fix_vof_3d(time%dt, 0)

        call this%Calculate_mesh_3d(time)

        if (this%rezone%rezone_type /= 0) call this%advect%Calculate_advect_3d()

    end subroutine do_time_step_3d


    subroutine Calculate_thermodynamics(this)
        use omp_lib
        implicit none
        class (hydro_step_t) , intent(inout)   :: this

        integer :: vert_mass_scheme = 1 

        real(8), dimension(:, :, :), pointer :: x   
        real(8), dimension(:, :, :), pointer :: y   
        real(8), dimension(:, :, :), pointer :: vol 

        real(8), dimension(:, :, :), pointer :: inversed_vertex_mass  
        real(8), dimension(:, :, :), pointer :: pressure_sum          
        real(8), dimension(:, :, :), pointer :: temperature           
        real(8), dimension(:, :, :), pointer :: vertex_mass           
        real(8), dimension(:, :, :), pointer :: cell_mass             
        real(8), dimension(:, :, :), pointer :: sound_vel             
        real(8), dimension(:, :, :), pointer :: pressure              
        real(8), dimension(:, :, :), pointer :: density               
        real(8), dimension(:, :, :), pointer :: dp_drho               
        real(8), dimension(:, :, :), pointer :: dp_de                 
        real(8), dimension(:, :, :), pointer :: sie                   
        real(8), dimension(:, :, :), pointer :: vof                   

        real(8), dimension(:, :, :, :), pointer :: temperature_vof_old
        real(8), dimension(:, :, :, :), pointer :: temperature_vof
        real(8), dimension(:, :, :, :), pointer :: cell_mass_vof
        real(8), dimension(:, :, :, :), pointer :: pressure_vof
        real(8), dimension(:, :, :, :), pointer :: density_vof
        real(8), dimension(:, :, :, :), pointer :: dp_drho_vof
        real(8), dimension(:, :, :, :), pointer :: dp_de_vof
        real(8), dimension(:, :, :, :), pointer :: dt_de_vof
        real(8), dimension(:, :, :, :), pointer :: sie_vof
        real(8), dimension(:, :, :, :), pointer :: mat_vof
        real(8), dimension(:, :, :, :), pointer :: sound_vel_vof
        real(8), dimension(:, :, :), pointer :: reem,reem1,reem2,reem3

        real(8), dimension(:, :, :), allocatable :: dt_de_temp 
        integer,save  :: counter = 0
        real(8),save :: timer1 = 0
        real(8) :: tmp
        integer :: tmp_mat  
        integer :: i, j, k     










        tmp = omp_get_wtime()
        call this%total_sie%Exchange_virtual_space_nonblocking()
        call this%total_cell_mass%Exchange_virtual_space_nonblocking()
       ! call this%vertex_mass%Exchange_virtual_space_nonblocking()
       ! call this%inversed_vertex_mass%Exchange_virtual_space_nonblocking()
        call this%mat_id%Exchange_virtual_space_nonblocking()
        call this%num_mat_cells%Exchange_virtual_space_nonblocking()
        call this%total_vof%Exchange_virtual_space_nonblocking()
        call this%materials%vof%Exchange_virtual_space_nonblocking()
        call this%materials%sie%Exchange_virtual_space_nonblocking()
        call this%materials%cell_mass%Exchange_virtual_space_nonblocking()


        call this%total_pressure_sum  %Point_to_data(pressure_sum)
        call this%total_temperature   %Point_to_data(temperature)
        call this%total_cell_mass     %Point_to_data(cell_mass)
        call this%total_sound_vel     %Point_to_data(sound_vel)
        call this%total_pressure      %Point_to_data(pressure)
        call this%total_density       %Point_to_data(density)
        call this%total_dp_drho       %Point_to_data(dp_drho)
        call this%total_dp_de         %Point_to_data(dp_de)
        call this%total_volume        %Point_to_data(vol)
        call this%total_vof           %Point_to_data(vof)
        call this%total_sie           %Point_to_data(sie)

        call this%materials%temperature%Point_to_data(temperature_vof)
        call this%materials%temperature_old%Point_to_data(temperature_vof_old)
        call this%materials%pressure   %Point_to_data(pressure_vof)
        call this%materials%dp_drho    %Point_to_data(dp_drho_vof)
        call this%materials%dp_de      %Point_to_data(dp_de_vof)
        call this%materials%cell_mass      %Point_to_data(cell_mass_vof)
        call this%materials%sound_vel      %Point_to_data(sound_vel_vof)
        call this%materials%dt_de          %Point_to_data(dt_de_vof)
        call this%materials%vof            %Point_to_data(mat_vof)


        call this%Calculate_density(this%total_volume)


        sound_vel   = 1d-20
        temperature = 0d0
        pressure    = 0d0
        dp_drho     = 0d0
        dp_de       = 0d0

        do k = 1, this%nz
            do j = 1, this%ny
                do i = 1, this%nx
                    do tmp_mat = 1, this%nmats

                        if (vof(i, j, k) >= this%emf) then
                            temperature_vof(tmp_mat, i, j, k) = 0d0
                            pressure_vof   (tmp_mat, i, j, k) = 0d0
                            dp_drho_vof    (tmp_mat, i, j, k) = 0d0
                            dp_de_vof      (tmp_mat, i, j, k) = 0d0
                        end if
                    end do
                end do
            end do
        end do



        allocate(dt_de_temp(this%nx, this%ny, this%nz))
        dt_de_temp = 0d0
        call this%total_sie%Exchange_end()
        call this%total_cell_mass%Exchange_end()
!        call this%vertex_mass%Exchange_end()
!        call this%inversed_vertex_mass%Exchange_end()
        call this%mat_id%Exchange_end()
        call this%num_mat_cells%Exchange_end()
        call this%total_vof%Exchange_end()
        call this%materials%vof%Exchange_end()
        call this%materials%sie%Exchange_end()
        call this%materials%cell_mass%Exchange_end()

        call this%materials%Apply_eos(this%nx, this%ny, this%nz, this%emf, .true.)



        do k = 1, this%nz
            do j = 1, this%ny
                do i = 1, this%nx
                    do tmp_mat = 1, this%nmats
                        if (mat_vof(tmp_mat, i, j, k) <= this%emf) cycle
                        temperature(i, j, k) = temperature(i, j, k) + &
                            temperature_vof_old(tmp_mat, i, j, k) * cell_mass_vof(tmp_mat, i, j, k) / (dt_de_vof(tmp_mat, i, j, k) + 1d-30)
                        sound_vel  (i, j, k) = sound_vel  (i, j, k) + sound_vel_vof(tmp_mat, i, j, k) * mat_vof(tmp_mat, i, j, k)
                        pressure   (i, j, k) = pressure   (i, j, k) + pressure_vof (tmp_mat, i, j, k) * mat_vof(tmp_mat, i, j, k)
                        dp_drho    (i, j, k) = dp_drho    (i, j, k) + dp_drho_vof  (tmp_mat, i, j, k) * mat_vof(tmp_mat, i, j, k)
                        dp_de      (i, j, k) = dp_de      (i, j, k) + cell_mass_vof(tmp_mat, i, j, k) / (dp_de_vof(tmp_mat, i, j, k) + 1d-30)
                        dt_de_temp (i, j, k) = dt_de_temp (i, j, k) + cell_mass_vof(tmp_mat, i, j, k) / (dt_de_vof(tmp_mat, i, j ,k) + 1d-30)
                    end do
                end do
            end do
        end do

        call this%total_pressure%Exchange_virtual_space_nonblocking()
        call this%total_density%Apply_boundary(.false.)

        do k = 1, this%nz
            do j = 1, this%ny
                do i = 1, this%nx
                    if (vof(i, j, k) <= this%emf) cycle
                    temperature(i, j, k) = temperature(i, j, k) / dt_de_temp(i, j, k)
                    dp_drho    (i, j, k) = dp_drho    (i, j, k) / (vof      (i, j, k) + 1d-30)
                    dp_de      (i, j, k) = cell_mass  (i, j, k) / dp_de     (i, j, k)
                end do
            end do
        end do

        deallocate(dt_de_temp)

        call this%total_pressure%Exchange_end()
        do k = 0, this%nzp
            do j = 0, this%nyp
                do i = 0, this%nxp
                    pressure_sum(i,j,k) = pressure(i,j,k)
                end do
            end do
        end do








        call this%total_density%Exchange_end()

        call this%num_mat_cells%Apply_boundary(.false.)
        call this%total_volume %Apply_boundary(.false.)
        call this%total_sie    %Apply_boundary(.false.)
        call this%total_vof    %Apply_boundary(.false.)
        call this%materials%cell_mass%Apply_boundary(.false.)
        call this%materials%density  %Apply_boundary(.false.)
        call this%materials%sie      %Apply_boundary(.false.)
        call this%materials%vof      %Apply_boundary(.false.)


        if (this%mesh%dimension == 2) then
            call this%vertex_mass%Calculate_vertex_mass_2d(this%mesh%coordinates, this%total_density, this%total_cell_mass, &
                0, this%mesh%cyl)
        else if (this%mesh%dimension == 3) then
            call this%vertex_mass%Calculate_vertex_mass_3d(this%mesh%coordinates, this%total_density, this%total_cell_mass)
        end if

        call this%num_mat_cells%Exchange_end()
        call this%total_volume %Exchange_end()
        call this%total_sie    %Exchange_end()
        call this%total_vof    %Exchange_end()
        call this%materials%cell_mass%Exchange_end()
        call this%materials%density  %Exchange_end()
        call this%materials%sie      %Exchange_end()
        call this%materials%vof      %Exchange_end()

        call this%vertex_mass%Exchange_virtual_space_blocking()

        call this%Calculate_inversed_vertex_mass()
        call this%inversed_vertex_mass%Exchange_virtual_space_blocking()
    end subroutine Calculate_thermodynamics

    subroutine Calculate_density(this, volume)
        class (hydro_step_t), intent(in out) :: this
        type (volume_t), pointer, intent(in) :: volume         

        real(8), dimension(:, :, :), pointer :: cell_mass      
        real(8), dimension(:, :, :), pointer :: density        
        real(8), dimension(:, :, :), pointer :: vol            
        real(8), dimension(:, :, :), pointer :: vof            
        real(8), dimension(:, :, :, :), pointer :: cell_mass_vof
        real(8), dimension(:, :, :, :), pointer :: density_vof
        real(8), dimension(:, :, :, :), pointer :: mat_vof

        integer :: i, j, k     
        integer :: tmp_mat  

        call this%total_density  %Point_to_data(density)
        call this%total_cell_mass%Point_to_data(cell_mass)
        call this%total_vof      %Point_to_data(vof)
        call this%materials%cell_mass%Point_to_data(cell_mass_vof)
        call this%materials%density  %Point_to_data(density_vof)
        call this%materials%vof      %Point_to_data(mat_vof)
        call volume%Point_to_data(vol)

        do k = 1, this%nz
            do j = 1, this%ny
                do i = 1, this%nx
                    density(i, j, k) = cell_mass(i, j, k) / vol(i, j, k)
                end do
            end do
        end do

        do k = 1, this%nz
            do j = 1, this%ny
                do i = 1, this%nx
                    do tmp_mat = 1, this%nmats
                        density_vof(tmp_mat, i, j, k) = 0d0
                        if (mat_vof(tmp_mat, i, j, k) > this%emf) then
                            density_vof(tmp_mat, i, j, k) = cell_mass_vof(tmp_mat, i, j, k) / (vol(i, j, k) * mat_vof(tmp_mat, i, j, k))
!write(*,*) "Calculating", tmp_mat,i,j,density_vof(tmp_mat,i,j,k),cell_mass_vof(tmp_mat,i,j,k), mat_vof(tmp_mat,i,j,k)
                        end if
                    end do
                end do
            end do
        end do
!        stop
    end subroutine Calculate_density

    subroutine Calculate_inversed_vertex_mass(this)
        class (hydro_step_t), intent(in out) :: this 

        real(8), dimension(:, :, :), pointer :: inversed_vertex_mass 
        real(8), dimension(:, :, :), pointer :: vertex_mass          

        integer :: i, j, k     

        call this%inversed_vertex_mass%Point_to_data(inversed_vertex_mass)
        call this%vertex_mass   %Point_to_data(vertex_mass)
        do k = 1, this%nzp
            do j = 1, this%nyp
                do i = 1, this%nxp
                    if (vertex_mass(i, j, k) > this%emfm) then
                        inversed_vertex_mass(i, j, k) = 1d0 / vertex_mass(i, j, k)
                    else
                        inversed_vertex_mass(i, j, k) = 0
                    end if
                end do
            end do
        end do
    end subroutine Calculate_inversed_vertex_mass

    subroutine Calculate_acceleration_2d(this, wilkins_scheme, dt_mid)
        class (hydro_step_t), intent(in out) :: this
        integer, intent(in) :: wilkins_scheme 
        real(8), intent(in) :: dt_mid         

        real(8), dimension(:, :, :), pointer :: acceleration_x        
        real(8), dimension(:, :, :), pointer :: acceleration_y        
        real(8), dimension(:, :, :), pointer :: r_factor              
        real(8), dimension(:, :, :), pointer :: x                     
        real(8), dimension(:, :, :), pointer :: y                     

        real(8), dimension(:, :, :), pointer :: inversed_vertex_mass  
        real(8), dimension(:, :, :), pointer :: pressure_sum          

        real(8) :: r_i_j, r_ip_j, r_i_jp, r_im_j, r_i_jm  
        integer :: i, j  








        call this%acceleration        %Point_to_data(acceleration_x, acceleration_y)
        call this%mesh                %Point_to_r_factor(r_factor)
        call this%mesh                %Point_to_data(x, y)
        call this%total_pressure_sum  %Point_to_data(pressure_sum)
        call this%inversed_vertex_mass%Point_to_data(inversed_vertex_mass)


        do j = 1, this%nyp
            do i = 1, this%nxp
                r_i_j  = dble(1 - wilkins_scheme) * r_factor(i    , j    , 1) + dble(wilkins_scheme)
                r_ip_j = dble(1 - wilkins_scheme) * r_factor(i + 1, j    , 1) + dble(wilkins_scheme)
                r_i_jp = dble(1 - wilkins_scheme) * r_factor(i    , j + 1, 1) + dble(wilkins_scheme)
                r_im_j = dble(1 - wilkins_scheme) * r_factor(i - 1, j    , 1) + dble(wilkins_scheme)
                r_i_jm = dble(1 - wilkins_scheme) * r_factor(i    , j - 1, 1) + dble(wilkins_scheme)
                acceleration_x(i, j, 1) = 0.5d0 * dt_mid * r_i_j * inversed_vertex_mass(i, j, 1) * ( &
                    pressure_sum(i - 1, j    , 1) * (y(i    , j + 1, 1) - y(i - 1, j    , 1))    &
                    - pressure_sum(i    , j    , 1) * (y(i    , j + 1, 1) - y(i + 1, j    , 1))    &
                    - pressure_sum(i    , j - 1, 1) * (y(i + 1, j    , 1) - y(i    , j - 1, 1))    &
                    + pressure_sum(i - 1, j - 1, 1) * (y(i - 1, j    , 1) - y(i    , j - 1, 1)))
                acceleration_y(i, j, 1) = 0.25d0 * dt_mid        * inversed_vertex_mass(i, j ,1) * (                  &
                    - pressure_sum(i - 1, j    , 1) * (x(i    , j + 1, 1) - x(i - 1, j    , 1)) * (r_im_j + r_i_jp) &
                    + pressure_sum(i    , j    , 1) * (x(i    , j + 1, 1) - x(i + 1, j    , 1)) * (r_i_jp + r_ip_j) &
                    + pressure_sum(i    , j - 1, 1) * (x(i + 1, j    , 1) - x(i    , j - 1, 1)) * (r_ip_j + r_i_jm) &
                    - pressure_sum(i - 1, j - 1, 1) * (x(i - 1, j    , 1) - x(i    , j - 1, 1)) * (r_i_jm + r_im_j))
            end do
        end do




        call this%acceleration%Apply_boundary(this%mesh%coordinates%data)

        return


    end subroutine Calculate_acceleration_2d

    subroutine Calculate_acceleration_3d(this, dt_mid)
        class (hydro_step_t), intent(in out) :: this
        real(8), intent(in) :: dt_mid         

        real(8), dimension(:, :, :), pointer :: acceleration_x        
        real(8), dimension(:, :, :), pointer :: acceleration_y        
        real(8), dimension(:, :, :), pointer :: acceleration_z        
        real(8), dimension(:, :, :), pointer :: x                     
        real(8), dimension(:, :, :), pointer :: y                     
        real(8), dimension(:, :, :), pointer :: z                     

        real(8), dimension(:, :, :), pointer :: inversed_vertex_mass  
        real(8), dimension(:, :, :), pointer :: pressure_sum          

        integer :: i, j, k            
        integer :: ip,im,jp,jm,kp,km            
        real(8) :: x1, x2, x3, x4, x5, x6  
        real(8) :: y1, y2, y3, y4, y5, y6  
        real(8) :: z1, z2, z3, z4, z5, z6  






        call this%acceleration%Point_to_data(acceleration_x, acceleration_y, acceleration_z)
        call this%mesh        %Point_to_data(x, y, z)

        call this%total_pressure_sum  %Point_to_data(pressure_sum)
        call this%inversed_vertex_mass%Point_to_data(inversed_vertex_mass)





        do k = 1, this%nzp
            do j = 1, this%nyp
                do i = 1, this%nxp
                    ip = i + 1
                    im = i - 1
                    jp = j + 1
                    jm = j - 1
                    kp = k + 1
                    km = k - 1

                    x1 = x(ip, j, k)
                    x2 = x(i, jp, k)
                    x3 = x(i, j, kp)
                    x4 = x(im, j, k)
                    x5 = x(i, jm, k)
                    x6 = x(i, j, km)

                    y1 = y(ip, j, k)
                    y2 = y(i, jp, k)
                    y3 = y(i, j, kp)
                    y4 = y(im, j, k)
                    y5 = y(i, jm, k)
                    y6 = y(i, j, km)

                    z1 = z(ip, j, k)
                    z2 = z(i, jp, k)
                    z3 = z(i, j, kp)
                    z4 = z(im, j, k)
                    z5 = z(i, jm, k)
                    z6 = z(i, j, km)

                    acceleration_x(i, j, k) = -0.25d0 * inversed_vertex_mass(i, j, k) * dt_mid * &
                        (pressure_sum(i  , j  , k ) * ((y2 - y1) * (z3 - z1) - (z2 - z1) * (y3 - y1)) &
                        +pressure_sum(i  , j  , km) * ((y6 - y1) * (z2 - z1) - (z6 - z1) * (y2 - y1)) &
                        +pressure_sum(im , j  , k ) * ((y3 - y4) * (z2 - z4) - (z3 - z4) * (y2 - y4)) &
                        +pressure_sum(im , j  , km) * ((y2 - y4) * (z6 - z4) - (z2 - z4) * (y6 - y4)) &
                        +pressure_sum(i  , jm , k ) * ((y1 - y5) * (z3 - z5) - (z1 - z5) * (y3 - y5)) &
                        +pressure_sum(i  , jm , km) * ((y6 - y5) * (z1 - z5) - (z6 - z5) * (y1 - y5)) &
                        +pressure_sum(im , jm , k ) * ((y5 - y4) * (z3 - z4) - (z5 - z4) * (y3 - y4)) &
                        +pressure_sum(im , jm , km) * ((y6 - y4) * (z5 - z4) - (z6 - z4) * (y5 - y4)))


                    acceleration_y(i, j, k) = -0.25d0 * inversed_vertex_mass(i, j, k) * dt_mid * &
                        (pressure_sum(i  , j , k ) * ((z2 - z1) * (x3 - x1) - (x2 - x1) * (z3 - z1)) &
                        +pressure_sum(i  , j , km) * ((z6 - z1) * (x2 - x1) - (x6 - x1) * (z2 - z1)) &
                        +pressure_sum(im , j , k ) * ((z3 - z4) * (x2 - x4) - (x3 - x4) * (z2 - z4)) &
                        +pressure_sum(im , j , km) * ((z2 - z4) * (x6 - x4) - (x2 - x4) * (z6 - z4)) &
                        +pressure_sum(i  , jm, k ) * ((z1 - z5) * (x3 - x5) - (x1 - x5) * (z3 - z5)) &
                        +pressure_sum(i  , jm, km) * ((z6 - z5) * (x1 - x5) - (x6 - x5) * (z1 - z5)) &
                        +pressure_sum(im , jm, k ) * ((z5 - z4) * (x3 - x4) - (x5 - x4) * (z3 - z4)) &
                        +pressure_sum(im , jm, km) * ((z6 - z4) * (x5 - x4) - (x6 - x4) * (z5 - z4)))

                    acceleration_z(i, j, k) = 0.25d0 * inversed_vertex_mass(i, j, k) * dt_mid * &
                        (-pressure_sum(i , j , k ) * ((x2 - x1) * (y3 - y1) - (y2 - y1) * (x3 - x1)) &
                        -pressure_sum(i , j , km) * ((x6 - x1) * (y2 - y1) - (y6 - y1) * (x2 - x1)) &
                        -pressure_sum(im, j , k ) * ((x3 - x4) * (y2 - y4) - (y3 - y4) * (x2 - x4)) &
                        -pressure_sum(im, j , km) * ((x2 - x4) * (y6 - y4) - (y2 - y4) * (x6 - x4)) &
                        -pressure_sum(i , jm, k ) * ((x1 - x5) * (y3 - y5) - (y1 - y5) * (x3 - x5)) &
                        -pressure_sum(i , jm, km) * ((x6 - x5) * (y1 - y5) - (y6 - y5) * (x1 - x5)) &
                        -pressure_sum(im, jm, k ) * ((x5 - x4) * (y3 - y4) - (y5 - y4) * (x3 - x4)) &
                        -pressure_sum(im, jm, km) * ((x6 - x4) * (y5 - y4) - (y6 - y4) * (x5 - x4)))



                end do
            end do
        end do




        call this%acceleration%Apply_boundary(this%mesh%coordinates%data, is_blocking=.false.)


        call this%velocity%Calculate_derivatives(this%mesh%coordinates, &
            this%total_vof, this%total_volume, this%nx, this%ny, this%nz, this%emf)
    end subroutine Calculate_acceleration_3d

    subroutine Calculate_artificial_viscosity_2d(this, time, wilkins_scheme)
        use geometry_module, only : Quad_area
        implicit none
        class (hydro_step_t), intent(in out) :: this
        class (time_t)      , intent(in out) :: time

        integer, intent(in)  :: wilkins_scheme   

        real(8), dimension(:, :, :), pointer :: inversed_vertex_mass 
        real(8), dimension(:, :, :), pointer :: acceleration_x       
        real(8), dimension(:, :, :), pointer :: acceleration_y       
        real(8), dimension(:, :, :), pointer :: velocity_x           
        real(8), dimension(:, :, :), pointer :: velocity_y           
        real(8), dimension(:, :, :), pointer :: sound_vel            
        real(8), dimension(:, :, :), pointer :: pressure             
        real(8), dimension(:, :, :), pointer :: r_factor             
        real(8), dimension(:, :, :), pointer :: density              
        real(8), dimension(:, :, :), pointer :: a_visc               
        real(8), dimension(:, :, :), pointer :: vof                  
        real(8), dimension(:, :, :), pointer :: x                    
        real(8), dimension(:, :, :), pointer :: y                    

        real(8) :: area                        
        real(8) :: inv_area                    
        real(8) :: du_dx, du_dy, dv_dx, dv_dy  
        real(8) :: tot_acc_i ,tot_acc_j        
        real(8) :: tot_acc                     
        real(8) :: cos_acc, sin_acc            
        real(8) :: x_avg, y_avg                
        real(8) :: x1, x2, x3, x4              
        real(8) :: y1, y2, y3, y4              
        real(8) :: u1, u2, u3, u4              
        real(8) :: v1, v2, v3, v4              
        real(8) :: dt_cour_temp                
        real(8) :: linear_visc_fac_temp        
        real(8) :: quad_visc_fac_temp          
        real(8) :: eff_length                  
        real(8) :: eff_vel_div                 
        real(8) :: dd1, dd2, dd3, dd4          
        real(8) :: r_i_j, r_i_jp, r_ip_j       
        real(8) :: r_i_jm, r_im_j              
        real(8) :: linear_visc_fac             
        real(8) :: quad_visc_fac               

        integer :: i, j                        

        character(len=*), parameter :: fmt_1002 = "('NEGATIVE SVELS=',1PE11.3,'I=',I5,'J=',I5)", &
            fmt_1001 = "('NEGATIVE CMAX=',1PE11.3,'I=',I5,'J=',I5,'SVELS=',1PE11.3,'AL=',1PE11.3)", &
            fmt_1101 = "('I_SPHERE=',I7)"

        call this%acceleration   %Point_to_data(acceleration_x, acceleration_y)
        call this%velocity       %Point_to_data(velocity_x, velocity_y)
        call this%mesh           %Point_to_data(x, y)
        call this%mesh           %Point_to_r_factor(r_factor)

        call this%inversed_vertex_mass%Point_to_data(inversed_vertex_mass)
        call this%total_pressure %Point_to_data(pressure)
        call this%total_density  %Point_to_data(density)
        call this%total_sound_vel%Point_to_data(sound_vel)
        call this%total_vof      %Point_to_data(vof)
        call this%a_visc         %Point_to_data(a_visc)

        linear_visc_fac = this%a_visc%Get_linear_visc_factor()
        quad_visc_fac = this%a_visc%Get_quad_visc_factor()

        time%dt_cour = 1d20

        do j = 1, this%ny
            do i = 1, this%nx
                if (vof(i, j, 1) < this%emf .or. pressure(i, j, 1) == 0d0) then  
                    a_visc(i, j, 1) = 0d0
                else

                    x1 = x(i + 1, j    , 1)              
                    x2 = x(i + 1, j + 1, 1)
                    x3 = x(i    , j + 1, 1)
                    x4 = x(i    , j    , 1)

                    y1 = y(i + 1, j    , 1)
                    y2 = y(i + 1, j + 1, 1)
                    y3 = y(i    , j + 1, 1)
                    y4 = y(i    , j    , 1)

                    u1 = velocity_x(i + 1, j    , 1)
                    u2 = velocity_x(i + 1, j + 1, 1)
                    u3 = velocity_x(i    , j + 1, 1)
                    u4 = velocity_x(i    , j    , 1)

                    v1 = velocity_y(i + 1, j    , 1)
                    v2 = velocity_y(i + 1, j + 1, 1)
                    v3 = velocity_y(i    , j + 1, 1)
                    v4 = velocity_y(i    , j    , 1)

                    area = Quad_area(x4, y4, x1, y1, x2, y2, x3, y3)   


                    inv_area = 1d0 / area
                    du_dx  = 0.5d0 * inv_area * ((u2 - u4) * (y3 - y1) - (u1 - u3) * (y4 - y2))      
                    du_dy  = 0.5d0 * inv_area * ((u4 - u2) * (x3 - x1) - (u1 - u3) * (x2 - x4))      
                    dv_dx  = 0.5d0 * inv_area * ((v2 - v4) * (y3 - y1) - (v1 - v3) * (y4 - y2))
                    dv_dy  = 0.5d0 * inv_area * ((v4 - v2) * (x3 - x1) - (v1 - v3) * (x2 - x4))

                    tot_acc_i = acceleration_x(i  , j  , 1) + acceleration_x(i+1, j  , 1) + &           
                        acceleration_x(i+1, j+1, 1) + acceleration_x(i  , j+1, 1)
                    tot_acc_j = acceleration_y(i  , j  , 1) + acceleration_y(i+1, j  , 1) + &
                        acceleration_y(i+1, j+1, 1) + acceleration_y(i  , j+1, 1)
                    tot_acc   = sqrt(tot_acc_i * tot_acc_i + tot_acc_j * tot_acc_j)                

                    if (tot_acc < 1d-8) then           
                        eff_length  = sqrt(area)
                        eff_vel_div = 0d0
                    else
                        cos_acc = tot_acc_i / tot_acc
                        sin_acc = tot_acc_j / tot_acc

                        eff_vel_div = du_dx * cos_acc * cos_acc + dv_dy * sin_acc * sin_acc + (du_dy + dv_dx) * sin_acc * cos_acc   
                        x_avg = 0.25d0 * (x1 + x2 + x3 + x4)
                        y_avg = 0.25d0 * (y1 + y2 + y3 + y4)
                        dd1   = abs((x1 - x_avg) * sin_acc - (y1 - y_avg) * cos_acc)
                        dd2   = abs((x2 - x_avg) * sin_acc - (y2 - y_avg) * cos_acc)
                        dd3   = abs((x3 - x_avg) * sin_acc - (y3 - y_avg) * cos_acc)
                        dd4   = abs((x4 - x_avg) * sin_acc - (y4 - y_avg) * cos_acc)
                        eff_length = 2d0 * area / (dd1 + dd2 + dd3 + dd4)


                    end if



                    quad_visc_fac_temp   = quad_visc_fac
                    linear_visc_fac_temp = linear_visc_fac





                    a_visc(i, j, 1) = density(i, j, 1) * eff_length * &    
                        (quad_visc_fac_temp * eff_length * eff_vel_div - linear_visc_fac_temp * sqrt(sound_vel(i, j, 1)))


                    a_visc(i, j, 1) = min(0d0, eff_vel_div) * a_visc(i, j, 1) * vof(i, j, 1)  






                    dt_cour_temp = eff_length * eff_length / sound_vel(i, j, 1)   
                    if (dt_cour_temp < time%dt_cour) then
                        time%dt_cour  = dt_cour_temp
                        time%i_cour   = i
                        time%j_cour   = j
                    end if
                end if



            end do
        end do




        time%dt_cour = sqrt(time%dt_cour) / time%dt_cour_fac
        call this%a_visc%Exchange_virtual_space_blocking()









    end subroutine Calculate_artificial_viscosity_2d


    subroutine Calculate_artificial_viscosity_3d(this, time)
        use geometry_module,  only : Line_length_in_cell
        use constants_module, only : THIRD
        use general_utils_module, only : Deter
        implicit none
        class (hydro_step_t), intent(in out) :: this
        class (time_t)      , intent(in out) :: time


        real(8), dimension(:, :, :), pointer :: inversed_vertex_mass 
        real(8), dimension(:, :, :), pointer :: acceleration_x       
        real(8), dimension(:, :, :), pointer :: acceleration_y       
        real(8), dimension(:, :, :), pointer :: acceleration_z       
        real(8), dimension(:, :, :), pointer :: velocity_x           
        real(8), dimension(:, :, :), pointer :: velocity_y           
        real(8), dimension(:, :, :), pointer :: velocity_z           
        real(8), dimension(:, :, :), pointer :: sound_vel            
        real(8), dimension(:, :, :), pointer :: pressure             
        real(8), dimension(:, :, :), pointer :: density              
        real(8), dimension(:, :, :), pointer :: a_visc               
        real(8), dimension(:, :, :), pointer :: vof                  
        real(8), dimension(:, :, :), pointer :: x                    
        real(8), dimension(:, :, :), pointer :: y                    
        real(8), dimension(:, :, :), pointer :: z                    

        real(8), dimension(:, :, :), pointer :: dvel_x_dx   
        real(8), dimension(:, :, :), pointer :: dvel_x_dy   
        real(8), dimension(:, :, :), pointer :: dvel_x_dz   
        real(8), dimension(:, :, :), pointer :: dvel_y_dx   
        real(8), dimension(:, :, :), pointer :: dvel_y_dy   
        real(8), dimension(:, :, :), pointer :: dvel_y_dz   
        real(8), dimension(:, :, :), pointer :: dvel_z_dx   
        real(8), dimension(:, :, :), pointer :: dvel_z_dy   
        real(8), dimension(:, :, :), pointer :: dvel_z_dz   

        real(8) :: x1, x2, x3, x4, x5, x6, x7, x8  
        real(8) :: y1, y2, y3, y4, y5, y6, y7, y8  
        real(8) :: z1, z2, z3, z4, z5, z6, z7, z8  
        real(8) :: r1, r2, r3, r4, r5, r6, r7, r8  

        real(8) :: avg_acc_x             
        real(8) :: avg_acc_y             
        real(8) :: avg_acc_z             
        real(8) :: tot_acc               
        real(8) :: cos_acc, sin_acc      
        real(8) :: x_avg, y_avg, z_avg   
        real(8) :: eff_length            
        real(8) :: eff_vel_div           
        real(8) :: dt_cour_temp          
        real(8) :: linear_visc_fac_temp  
        real(8) :: quad_visc_fac_temp    
        real(8) :: quad_visc_fac         
        integer :: i_sphere              
        integer :: i_sphere2             
        real(8), dimension(:,:,:), pointer :: linear_visc_fac       
        integer :: i, j, k  
        integer :: ip, jp, kp
        integer :: mesh_type

        call this%a_visc%Update_visc_factors()
        linear_visc_fac => this%a_visc%xl_visc_mat
        quad_visc_fac   = this%a_visc%Get_quad_visc_factor()
        i_sphere        = this%a_visc%to_radial_index_sphere
        i_sphere2       = this%a_visc%from_radial_index_sphere
        mesh_type       = this%mesh%mesh_type

        call this%velocity%dvelocity_x_dx%Point_to_data(dvel_x_dx)
        call this%velocity%dvelocity_x_dy%Point_to_data(dvel_x_dy)
        call this%velocity%dvelocity_x_dz%Point_to_data(dvel_x_dz)
        call this%velocity%dvelocity_y_dx%Point_to_data(dvel_y_dx)
        call this%velocity%dvelocity_y_dy%Point_to_data(dvel_y_dy)
        call this%velocity%dvelocity_y_dz%Point_to_data(dvel_y_dz)
        call this%velocity%dvelocity_z_dx%Point_to_data(dvel_z_dx)
        call this%velocity%dvelocity_z_dy%Point_to_data(dvel_z_dy)
        call this%velocity%dvelocity_z_dz%Point_to_data(dvel_z_dz)

        call this%acceleration   %Point_to_data(acceleration_x, acceleration_y, acceleration_z)
        call this%velocity       %Point_to_data(velocity_x, velocity_y, velocity_z)
        call this%mesh           %Point_to_data(x, y, z)

        call this%total_pressure %Point_to_data(pressure)
        call this%total_density  %Point_to_data(density)
        call this%total_sound_vel%Point_to_data(sound_vel)
        call this%total_vof      %Point_to_data(vof)
        call this%a_visc         %Point_to_data(a_visc)

        call this%inversed_vertex_mass%Point_to_data(inversed_vertex_mass)





        time%dt_cour = 1d20
        call this%acceleration%Exchange_end()

        do k = 1, this%nz
            do j = 1, this%ny
                do i = 1, this%nx
                    if (vof(i, j, k) < this%emf .or. pressure(i, j, k) == 0d0) then  
                        a_visc(i, j, k) = 0d0
                    else
                        ip = i + 1
                        jp = j + 1
                        kp = k + 1

                        avg_acc_x = (acceleration_x(i , j , k ) + acceleration_x(ip , j , k  ) + acceleration_x(i , j+1, k  ) &
                            +acceleration_x(ip, jp, k ) + acceleration_x(i  , j , kp)  + acceleration_x(ip, j  , k+1) &
                            +acceleration_x(i , jp, kp) + acceleration_x(ip , jp, kp)) * 0.125d0

                        avg_acc_y = (acceleration_y(i , j , k ) + acceleration_y(ip, j , k )  + acceleration_y(i , jp,k ) &
                            +acceleration_y(ip, jp, k ) + acceleration_y(i , j , kp)  + acceleration_y(ip, j ,kp) &
                            +acceleration_y(i , jp, kp) + acceleration_y(ip, jp, kp)) * 0.125d0

                        avg_acc_z = (acceleration_z(i , j , k ) + acceleration_z(ip, j , k  )  + acceleration_z(i , jp, k  ) &
                            +acceleration_z(ip, jp, k ) + acceleration_z(i , j , kp )  + acceleration_z(ip, j , kp ) &
                            +acceleration_z(i , jp, kp) + acceleration_z(ip, jp, kp )) * 0.125d0

                        tot_acc = sqrt(avg_acc_x * avg_acc_x + avg_acc_y * avg_acc_y + avg_acc_z * avg_acc_z)


                        if (tot_acc < 1d-8 .or. (mesh_type /= 2 .and. k < i_sphere) ) then 


                            r1 = sqrt(x(i  , j  , kp) ** 2d0 + y(i  , j  , kp) ** 2d0 + z(i  , j  , kp) ** 2d0)
                            r2 = sqrt(x(ip, j  , kp) ** 2d0 + y(ip, j  , kp) ** 2d0 + z(ip, j  , kp) ** 2d0)
                            r3 = sqrt(x(i  , jp, kp) ** 2d0 + y(i  , jp, kp) ** 2d0 + z(i  , jp, kp) ** 2d0)
                            r4 = sqrt(x(ip, jp, kp) ** 2d0 + y(ip, jp, kp) ** 2d0 + z(ip, jp, kp) ** 2d0)
                            r5 = sqrt(x(i  , j  , k  ) ** 2d0 + y(i  , j  , k  ) ** 2d0 + z(i  , j  , k  ) ** 2d0)
                            r6 = sqrt(x(ip, j  , k  ) ** 2d0 + y(ip, j  , k  ) ** 2d0 + z(ip, j  , k  ) ** 2d0)
                            r7 = sqrt(x(i  , jp, k  ) ** 2d0 + y(i  , jp, k  ) ** 2d0 + z(i  , jp, k  ) ** 2d0)
                            r8 = sqrt(x(ip, jp, k  ) ** 2d0 + y(ip, jp, k  ) ** 2d0 + z(ip, jp, k  ) ** 2d0)
                            eff_length = (r1 + r2 + r3 + r4 - r5 - r6 - r7 - r8) / 4d0
                        else
                            x1 = x(i, j, k)
                            x2 = x(ip, j, k)
                            x3 = x(ip, jp, k)
                            x4 = x(i, jp, k)
                            x5 = x(i, j, kp)
                            x6 = x(ip, j, kp)
                            x7 = x(ip, jp, kp)
                            x8 = x(i, jp, kp)

                            y1 = y(i, j, k)
                            y2 = y(ip, j, k)
                            y3 = y(ip, jp, k)
                            y4 = y(i, jp, k)
                            y5 = y(i, j, kp)
                            y6 = y(ip, j, kp)
                            y7 = y(ip, jp, kp)
                            y8 = y(i, jp, kp)

                            z1 = z(i, j, k)
                            z2 = z(ip, j, k)
                            z3 = z(ip, jp, k)
                            z4 = z(i, jp, k)
                            z5 = z(i, j, kp)
                            z6 = z(ip, j, kp)
                            z7 = z(ip, jp, kp)
                            z8 = z(i, jp, kp)

                            x_avg = (x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8) * 0.125d0
                            y_avg = (y1 + y2 + y3 + y4 + y5 + y6 + y7 + y8) * 0.125d0
                            z_avg = (z1 + z2 + z3 + z4 + z5 + z6 + z7 + z8) * 0.125d0

                            eff_length = Line_length_in_cell(x1, y1, z1 ,x2, y2, z2 ,x3, y3, z3 ,x4, y4, z4 , &
                                x5, y5, z5 ,x6, y6, z6 ,x7, y7, z7 ,x8, y8, z8 , &
                                avg_acc_x, avg_acc_y, avg_acc_z, x_avg, y_avg, z_avg)

                        end if
                        if (tot_acc < 1d-8) then
                            eff_vel_div = 0d0
                        else
                            eff_vel_div = (dvel_x_dx(i, j, k) * avg_acc_x ** 2 + &
                                dvel_y_dy(i, j, k) * avg_acc_y ** 2 + &
                                dvel_z_dz(i, j, k) * avg_acc_z ** 2 + &
                                (dvel_x_dy(i, j, k) + dvel_y_dx(i, j, k)) * avg_acc_x * avg_acc_y + &
                                (dvel_x_dz(i, j, k) + dvel_z_dx(i, j, k)) * avg_acc_x * avg_acc_z + &
                                (dvel_y_dz(i, j, k) + dvel_z_dy(i, j, k)) * avg_acc_y * avg_acc_z) / &
                                (avg_acc_x ** 2 + avg_acc_y ** 2 + avg_acc_z ** 2 + 1d-20)
                        end if





                        quad_visc_fac_temp   = quad_visc_fac
                        linear_visc_fac_temp = linear_visc_fac(i,j,k)



                        a_visc(i, j, k) = density(i, j, k) * eff_length * &
                            (quad_visc_fac_temp * eff_length * eff_vel_div &
                            -linear_visc_fac_temp * sqrt(sound_vel(i ,j, k)))


                        a_visc(i, j, k) = min(0d0, eff_vel_div) * a_visc(i, j, k)
                        dt_cour_temp = eff_length * eff_length / sound_vel(i ,j, k)

                        if (dt_cour_temp < time%dt_cour) then
                            time%dt_cour = dt_cour_temp
                            time%i_cour = i
                            time%j_cour = j
                            time%k_cour = k
                        end if
                    end if
                end do
            end do
        end do

        call this%a_visc%Exchange_virtual_space_nonblocking()
        time%dt_cour = sqrt(time%dt_cour) / time%dt_cour_fac



    end subroutine Calculate_artificial_viscosity_3d


    subroutine Calculate_velocity_2d(this, dt_mid, wilkins_scheme)
        class (hydro_step_t), intent(in out) :: this
        real(8)             , intent(in)     :: dt_mid
        integer             , intent(in)     :: wilkins_scheme   

        real(8), dimension(:, :, :), pointer :: inversed_vertex_mass 
        real(8), dimension(:, :, :), pointer :: acceleration_x       
        real(8), dimension(:, :, :), pointer :: acceleration_y       
        real(8), dimension(:, :, :), pointer :: velocity_x           
        real(8), dimension(:, :, :), pointer :: velocity_y           
        real(8), dimension(:, :, :), pointer :: r_factor             
        real(8), dimension(:, :, :), pointer :: a_visc               
        real(8), dimension(:, :, :), pointer :: x                    
        real(8), dimension(:, :, :), pointer :: y                    


        real(8) :: r_i_j, r_i_jp, r_ip_j       
        real(8) :: r_i_jm, r_im_j              
        integer :: i, j                        


        call this%acceleration   %Point_to_data(acceleration_x, acceleration_y)
        call this%velocity       %Point_to_data(velocity_x, velocity_y)
        call this%mesh           %Point_to_data(x, y)
        call this%mesh           %Point_to_r_factor(r_factor)

        call this%inversed_vertex_mass%Point_to_data(inversed_vertex_mass)
        call this%a_visc         %Point_to_data(a_visc)


        do j = 1, this%nyp     
            do i = 1, this%nxp

                r_i_j  = (dble(1 - wilkins_scheme)) * r_factor(i    , j    , 1) + dble(wilkins_scheme)
                r_ip_j = (dble(1 - wilkins_scheme)) * r_factor(i + 1, j    , 1) + dble(wilkins_scheme)
                r_i_jp = (dble(1 - wilkins_scheme)) * r_factor(i    , j + 1, 1) + dble(wilkins_scheme)
                r_im_j = (dble(1 - wilkins_scheme)) * r_factor(i - 1, j    , 1) + dble(wilkins_scheme)
                r_i_jm = (dble(1 - wilkins_scheme)) * r_factor(i    , j - 1, 1) + dble(wilkins_scheme)



                velocity_x(i, j, 1) = velocity_x(i, j, 1) + acceleration_x(i, j, 1) +                  & 
                    0.5d0 * dt_mid * r_i_j * inversed_vertex_mass(i, j, 1) * (    & 
                    a_visc(i - 1, j    , 1) * (y(i    , j + 1, 1) - y(i - 1, j    , 1)) - &
                    a_visc(i    , j    , 1) * (y(i    , j + 1, 1) - y(i + 1, j    , 1)) - &
                    a_visc(i    , j - 1, 1) * (y(i + 1, j    , 1) - y(i    , j - 1, 1)) + &
                    a_visc(i - 1, j - 1, 1) * (y(i - 1, j    , 1) - y(i    , j - 1, 1)))

                velocity_y(i, j, 1) = velocity_y(i, j, 1) + acceleration_y(i, j, 1) +                                         &
                    0.25d0 * dt_mid        * inversed_vertex_mass(i, j, 1) * (  &
                    -a_visc(i - 1, j    , 1) * (x(i    , j + 1, 1) - x(i - 1, j    , 1)) * (r_i_jp + r_im_j) + &
                    a_visc(i    , j    , 1) * (x(i    , j + 1, 1) - x(i + 1, j    , 1)) * (r_ip_j + r_i_jp) + &
                    a_visc(i    , j - 1, 1) * (x(i + 1, j    , 1) - x(i    , j - 1, 1)) * (r_i_jm + r_ip_j) - &
                    a_visc(i - 1, j - 1, 1) * (x(i - 1, j    , 1) - x(i    , j - 1, 1)) * (r_im_j + r_i_jm))
            end do
        end do










        call this%velocity%Apply_boundary(this%mesh%coordinates%data)









    end subroutine Calculate_velocity_2d


    subroutine Calculate_velocity_3d(this, dt_mid)
        class (hydro_step_t), intent(in out) :: this
        real(8)             , intent(in) :: dt_mid
        real(8), dimension(:, :, :), pointer :: inversed_vertex_mass 
        real(8), dimension(:, :, :), pointer :: acceleration_x       
        real(8), dimension(:, :, :), pointer :: acceleration_y       
        real(8), dimension(:, :, :), pointer :: acceleration_z       
        real(8), dimension(:, :, :), pointer :: velocity_x           
        real(8), dimension(:, :, :), pointer :: velocity_y           
        real(8), dimension(:, :, :), pointer :: velocity_z           
        real(8), dimension(:, :, :), pointer :: a_visc               
        real(8), dimension(:, :, :), pointer :: x                    
        real(8), dimension(:, :, :), pointer :: y                    
        real(8), dimension(:, :, :), pointer :: z                    
        integer :: i, j, k  
        integer :: ip,im,jp,jm,kp,km  
        real(8) :: x1, x2, x3, x4, x5, x6, x7, x8  
        real(8) :: y1, y2, y3, y4, y5, y6, y7, y8  
        real(8) :: z1, z2, z3, z4, z5, z6, z7, z8  
        real(8) :: r1, r2, r3, r4, r5, r6, r7, r8  

        call this%acceleration   %Point_to_data(acceleration_x, acceleration_y, acceleration_z)
        call this%velocity       %Point_to_data(velocity_x, velocity_y, velocity_z)
        call this%mesh           %Point_to_data(x, y, z)

        call this%inversed_vertex_mass%Point_to_data(inversed_vertex_mass)

        call this%a_visc         %Point_to_data(a_visc)

        call this%a_visc%Exchange_end()
        do k = 1, this%nzp
            do j = 1, this%nyp
                do i = 1, this%nxp
                    ip = i + 1
                    im = i - 1
                    jp = j + 1
                    jm = j - 1
                    kp = k + 1
                    km = k - 1

                    x1 = x(ip, j, k)
                    x2 = x(i, jp, k)
                    x3 = x(i, j, kp)
                    x4 = x(im, j, k)
                    x5 = x(i, jm, k)
                    x6 = x(i, j, km)

                    y1 = y(ip, j, k)
                    y2 = y(i, jp, k)
                    y3 = y(i, j, kp)
                    y4 = y(im, j, k)
                    y5 = y(i, jm, k)
                    y6 = y(i, j, km)

                    z1 = z(ip, j, k)
                    z2 = z(i, jp, k)
                    z3 = z(i, j, kp)
                    z4 = z(im, j, k)
                    z5 = z(i, jm, k)
                    z6 = z(i, j, km)

                    velocity_x(i, j, k) = velocity_x(i, j, k) + acceleration_x(i, j, k) + &
                        dt_mid * 0.25d0 * inversed_vertex_mass(i, j, k) * &
                        (-a_visc(i , j , k ) * ((y2 - y1) * (z3 - z1) - (z2 - z1) * (y3-y1)) &
                        -a_visc(i , j , km) * ((y6 - y1) * (z2 - z1) - (z6 - z1) * (y2-y1)) &
                        -a_visc(im, j , k ) * ((y3 - y4) * (z2 - z4) - (z3 - z4) * (y2-y4)) &
                        -a_visc(im, j , km) * ((y2 - y4) * (z6 - z4) - (z2 - z4) * (y6-y4)) &
                        -a_visc(i , jm, k ) * ((y1 - y5) * (z3 - z5) - (z1 - z5) * (y3-y5)) &
                        -a_visc(i , jm, km) * ((y6 - y5) * (z1 - z5) - (z6 - z5) * (y1-y5)) &
                        -a_visc(im, jm, k ) * ((y5 - y4) * (z3 - z4) - (z5 - z4) * (y3-y4)) &
                        -a_visc(im, jm, km) * ((y6 - y4) * (z5 - z4) - (z6 - z4) * (y5-y4)))



                    velocity_y(i, j, k) = velocity_y(i, j, k) + acceleration_y(i, j, k) + &
                        dt_mid * 0.25d0 * inversed_vertex_mass(i, j, k) * &
                        (-a_visc(i , j , k ) * ((z2 - z1) * (x3 - x1) - (x2 - x1) * (z3 - z1)) &
                        -a_visc(i , j , km) * ((z6 - z1) * (x2 - x1) - (x6 - x1) * (z2 - z1)) &
                        -a_visc(im, j , k ) * ((z3 - z4) * (x2 - x4) - (x3 - x4) * (z2 - z4)) &
                        -a_visc(im, j , km) * ((z2 - z4) * (x6 - x4) - (x2 - x4) * (z6 - z4)) &
                        -a_visc(i , jm, k ) * ((z1 - z5) * (x3 - x5) - (x1 - x5) * (z3 - z5)) &
                        -a_visc(i , jm, km) * ((z6 - z5) * (x1 - x5) - (x6 - x5) * (z1 - z5)) &
                        -a_visc(im, jm, k ) * ((z5 - z4) * (x3 - x4) - (x5 - x4) * (z3 - z4)) &
                        -a_visc(im, jm, km) * ((z6 - z4) * (x5 - x4) - (x6 - x4) * (z5 - z4)))

                    velocity_z(i, j, k) = velocity_z(i, j, k) + (acceleration_z(i, j, k) + &
                        dt_mid * 0.25d0 * inversed_vertex_mass(i, j, k) * &
                        (-a_visc(i , j , k ) * ((x2 - x1) * (y3 - y1) - (y2 - y1) * (x3 - x1)) &
                        -a_visc(i , j , km) * ((x6 - x1) * (y2 - y1) - (y6 - y1) * (x2 - x1)) &
                        -a_visc(im, j , k ) * ((x3 - x4) * (y2 - y4) - (y3 - y4) * (x2 - x4)) &
                        -a_visc(im, j , km) * ((x2 - x4) * (y6 - y4) - (y2 - y4) * (x6 - x4)) &
                        -a_visc(i , jm, k ) * ((x1 - x5) * (y3 - y5) - (y1 - y5) * (x3 - x5)) &
                        -a_visc(i , jm, km) * ((x6 - x5) * (y1 - y5) - (y6 - y5) * (x1 - x5)) &
                        -a_visc(im, jm, k ) * ((x5 - x4) * (y3 - y4) - (y5 - y4) * (x3 - x4)) &
                        -a_visc(im, jm, km) * ((x6 - x4) * (y5 - y4) - (y6 - y4) * (x5 - x4))))




                end do
            end do
        end do

        if (this%mesh%mesh_type /= 2) then
            call this%velocity%Impose_spherical_symmetry(this%mesh%coordinates)
        end if

        call this%velocity%Impose_no_move_3d(this%mesh%coordinates)

        call this%velocity%Apply_boundary(this%mesh%coordinates%data, is_blocking=.false.)

        call this%velocity%Calculate_derivatives(this%mesh%coordinates, this%total_vof, this%total_volume&
            , this%nx, this%ny, this%nz, this%emf)
            call this%velocity%Exchange_end()
    end subroutine Calculate_velocity_3d


    subroutine Calculate_energy_2d(this, teps, dt)
        use constants_module, only : THIRD
        implicit none
        class (hydro_step_t), intent(in out) :: this

        real(8), intent(in) :: teps        
        real(8), intent(in) :: dt          

        real(8), dimension(:, :, :), pointer :: velocity_x   
        real(8), dimension(:, :, :), pointer :: velocity_y   
        real(8), dimension(:, :, :), pointer :: r_factor     
        real(8), dimension(:, :, :), pointer :: x            
        real(8), dimension(:, :, :), pointer :: y            

        real(8), dimension(:, :, :), pointer :: nmats_in_cell

        real(8), dimension(:, :, :), pointer :: temperature  
        real(8), dimension(:, :, :), pointer :: cell_mass    
        real(8), dimension(:, :, :), pointer :: pressure     
        real(8), dimension(:, :, :), pointer :: a_visc       
        real(8), dimension(:, :, :), pointer :: dp_drho      
        real(8), dimension(:, :, :), pointer :: dp_de        
        real(8), dimension(:, :, :), pointer :: vof          
        real(8), dimension(:, :, :), pointer :: sie          
        real(8), dimension(:, :, :), pointer :: vol          

        real(8), dimension(:, :, :, :), pointer :: temperature_vof
        real(8), dimension(:, :, :, :), pointer :: cell_mass_vof1
        real(8), dimension(:, :, :, :), pointer :: pressure_vof
        real(8), dimension(:, :, :, :), pointer :: density_vof
        real(8), dimension(:, :, :, :), pointer :: dp_drho_vof
        real(8), dimension(:, :, :, :), pointer :: dp_de_vof
        real(8), dimension(:, :, :, :), pointer :: sie_vof
        real(8), dimension(:, :, :, :), pointer :: mat_vof

        integer :: i, j, tmp_mat        
        real(8) :: x1, x2, x3, x4       
        real(8) :: y1, y2, y3, y4       
        real(8) :: r1, r2, r3, r4       
        real(8) :: u1, u2, u3, u4       
        real(8) :: v1, v2, v3, v4       
        real(8) :: x1n, x2n, x3n, x4n   
        real(8) :: y1n, y2n, y3n, y4n   
        real(8) :: top_right_area       
        real(8) :: bottom_left_area     
        real(8) :: new_volume           
        real(8) :: old_volume           
        real(8) :: vol_diff             
        real(8) :: x24_diff, x31_diff   
        real(8) :: y24_diff, y31_diff   
        real(8) :: r24_avg , r31_avg    
        real(8) :: a_visc_temp          
        real(8) :: vol_diff_stress      
        real(8) :: sie_vof_temp         
        real(8) :: pressure_temp        
        real(8) :: dp_drho_temp         
        real(8) :: dp_de_temp           
        real(8) :: mass_vof_temp        
        real(8) :: vol_vof_temp         
        real(8) :: vol_diff_vof_temp    
        real(8) :: vol_diff_vof_stress_temp 
        real(8) :: sie_diff             
        real(8) :: stress_fac           

        call this%velocity         %Point_to_data(velocity_x, velocity_y)
        call this%mesh             %Point_to_data(x, y)
        call this%mesh             %Point_to_r_factor(r_factor)

        call this%total_temperature%Point_to_data(temperature)
        call this%total_pressure   %Point_to_data(pressure)
        call this%total_dp_drho    %Point_to_data(dp_drho)
        call this%total_dp_de      %Point_to_data(dp_de)
        call this%total_volume     %Point_to_data(vol)
        call this%total_cell_mass  %Point_to_data(cell_mass)
        call this%total_sie        %Point_to_data(sie)
        call this%total_vof        %Point_to_data(vof)
        call this%a_visc           %Point_to_data(a_visc)
        call this%materials%temperature%Point_to_data(temperature_vof)
        call this%materials%cell_mass  %Point_to_data(cell_mass_vof1)
        call this%materials%pressure   %Point_to_data(pressure_vof)
        call this%materials%density    %Point_to_data(density_vof)
        call this%materials%dp_drho    %Point_to_data(dp_drho_vof)
        call this%materials%dp_de      %Point_to_data(dp_de_vof)
        call this%materials%sie        %Point_to_data(sie_vof)
        call this%materials%vof        %Point_to_data(mat_vof)

        call this%num_mat_cells    %Point_to_data(nmats_in_cell)
        do j = 1, this%ny
            do i = 1, this%nx
                if (vof(i, j, 1) < this%emf) then   
                    sie(i, j, 1) = 0d0
                    temperature (i, j, 1) = teps
                    do tmp_mat = 1, this%nmats
                        temperature_vof(tmp_mat, i, j, 1) = 0d0
                    end do


                else
                    r1 = r_factor(i + 1, j, 1)
                    x1 = x(i + 1, j, 1)
                    y1 = y(i + 1, j, 1)
                    u1 = velocity_x(i + 1, j, 1) + velocity_x(i + 1, j, 1)  
                    v1 = velocity_y(i + 1, j, 1) + velocity_y(i + 1, j, 1)

                    r2 = r_factor(i + 1, j + 1, 1)
                    x2 = x(i + 1, j + 1, 1)
                    y2 = y(i + 1, j + 1, 1)
                    u2 = velocity_x(i + 1, j + 1, 1) + velocity_x(i + 1, j + 1, 1)
                    v2 = velocity_y(i + 1, j + 1, 1) + velocity_y(i + 1, j + 1, 1)

                    r3 = r_factor(i, j + 1, 1)
                    x3 = x(i, j + 1, 1)
                    y3 = y(i, j + 1, 1)
                    u3 = velocity_x(i, j + 1, 1) + velocity_x(i, j + 1, 1)
                    v3 = velocity_y(i, j + 1, 1) + velocity_y(i, j + 1, 1)

                    r4 = r_factor(i, j, 1)
                    x4 = x(i, j, 1)
                    y4 = y(i, j, 1)
                    u4 = velocity_x(i, j, 1) + velocity_x(i, j, 1)
                    v4 = velocity_y(i, j, 1) + velocity_y(i, j, 1)

                    y24_diff  = y2 - y4
                    y31_diff  = y3 - y1
                    x24_diff  = x2 - x4
                    x31_diff  = x3 - x1
                    r24_avg = 0.5d0 * (r2 + r4)
                    r31_avg = 0.5d0 * (r1 + r3)

                    if ((this%wilkins_scheme == 1) .and. (this%mesh%cyl == 1d0)) then

                        x1n = x1 + velocity_x(i + 1, j, 1) * dt
                        y1n = y1 + velocity_y(i + 1, j, 1) * dt
                        x2n = x2 + velocity_x(i + 1, j + 1, 1) * dt
                        y2n = y2 + velocity_y(i + 1, j + 1, 1) * dt
                        x3n = x3 + velocity_x(i, j + 1, 1) * dt
                        y3n = y3 + velocity_y(i, j + 1, 1) * dt
                        x4n = x4 + velocity_x(i, j, 1) * dt
                        y4n = y4 + velocity_y(i, j, 1) * dt
                        top_right_area   = (x3n - x2n) * (y1n - y2n) - (x1n - x2n) * (y3n - y2n)
                        bottom_left_area = (x1n - x4n) * (y3n - y4n) - (x3n - x4n) * (y1n - y4n)
                        new_volume       = ((x1n + x2n + x3n) * top_right_area + (x3n + x4n + x1n) * bottom_left_area) * THIRD
                        top_right_area   = (x3 - x2) * (y1 - y2) - (x1 - x2) * (y3 - y2)
                        bottom_left_area = (x1 - x4) * (y3 - y4) - (x3 - x4) * (y1 - y4)
                        old_volume       = ((x1 + x2 + x3) * top_right_area + (x3 + x4 + x1) * bottom_left_area) * THIRD
                        vol_diff         = (new_volume - old_volume) / dt * 2.d0

                    else
                        vol_diff = y24_diff * (u1 * r1 - u3 * r3) + y31_diff * (u2 * r2 - u4 * r4) - &   
                            r24_avg * x24_diff * (v1 - v3) - r31_avg * x31_diff * (v2 - v4)
                    end if

                    vol_diff_stress = 0d0
                    vol_diff        = vol_diff        * 0.25d0 * dt    
                    vol_diff_stress = vol_diff_stress * 0.25d0 * dt

                    if (vol_diff < 0d0) then      
                        a_visc_temp = a_visc(i, j, 1)
                    else
                        a_visc_temp = 0d0
                    end if

                    sie(i, j, 1) = 0d0     

                    do tmp_mat = 1, this%nmats



                        if (mat_vof(tmp_mat, i, j, 1) < this%emf) then
                            !                            sie_vof(i, j, 1) = 0d0
                            temperature_vof(tmp_mat, i, j, 1) = 0d0
                        else
                            !                            sie_vof_temp = sie_vof(i, j, 1)
                            if (nmats_in_cell(i, j, 1) == 1d0) then
                                pressure_temp = pressure(i, j, 1)



                                dp_de_temp   = dp_de  (i, j, 1)
                                dp_drho_temp = dp_drho(i, j, 1)
                            else                                         
                                pressure_temp = pressure_vof(tmp_mat, i, j, 1) * vof(i, j, 1)




                                dp_de_temp   = dp_de_vof  (tmp_mat, i, j, 1)
                                dp_drho_temp = dp_drho_vof(tmp_mat, i, j, 1)
                            end if
                            mass_vof_temp = cell_mass_vof1(tmp_mat, i, j, 1)
                            vol_vof_temp  = vol(i, j, 1) * mat_vof(tmp_mat, i, j, 1)



                            vol_diff_vof_temp = vol_diff * mat_vof(tmp_mat, i, j, 1)

                            stress_fac = 0d0
                            vol_diff_vof_stress_temp = vol_diff_stress * stress_fac * vof(i, j, 1)




                            sie_diff = - ((pressure_temp - &
                                (0.5d0 * dp_drho_temp * mass_vof_temp / vol_vof_temp ** 2 * vol_diff_vof_temp) + a_visc_temp) *&
                                vol_diff_vof_temp + vol_diff_vof_stress_temp) / &
                                (1d0 + 0.5d0 * dp_de_temp * vol_diff_vof_temp / (mass_vof_temp + 1d-30)) / (mass_vof_temp+1d-30)

                                                    sie_vof(tmp_mat,i, j, 1) = sie_vof(tmp_mat,i, j, 1) + sie_diff



                                                    sie(i, j, 1) = sie(i, j, 1) + sie_vof(tmp_mat,i, j, 1) * mass_vof_temp
                        end if   
                    end do





                    sie(i, j, 1) = sie(i, j, 1) / (cell_mass(i, j, 1) + 1d-30)     
                end if
            end do
        end do





        call this%total_sie%Exchange_virtual_space_blocking()
        call this%materials%sie%Exchange_virtual_space_blocking()
    end subroutine Calculate_energy_2d

    subroutine Calculate_energy_3d(this, teps, dt)
        implicit none
        class (hydro_step_t), intent(in out) :: this


        real(8), intent(in) :: teps        
        real(8), intent(in) :: dt          

        real(8), dimension(:, :, :), pointer :: velocity_x   
        real(8), dimension(:, :, :), pointer :: velocity_y   
        real(8), dimension(:, :, :), pointer :: velocity_z   
        real(8), dimension(:, :, :), pointer :: x            
        real(8), dimension(:, :, :), pointer :: y            
        real(8), dimension(:, :, :), pointer :: z            

        real(8), dimension(:, :, :), pointer :: nmats_in_cell

        real(8), dimension(:, :, :), pointer :: temperature  
        real(8), dimension(:, :, :), pointer :: cell_mass    
        real(8), dimension(:, :, :), pointer :: pressure     
        real(8), dimension(:, :, :), pointer :: a_visc       
        real(8), dimension(:, :, :), pointer :: dp_drho      
        real(8), dimension(:, :, :), pointer :: dp_de        
        real(8), dimension(:, :, :), pointer :: vof          
        real(8), dimension(:, :, :), pointer :: sie          
        real(8), dimension(:, :, :), pointer :: vol          

        real(8), dimension(:, :, :, :), pointer :: temperature_vof
        real(8), dimension(:, :, :, :), pointer :: cell_mass_vof
        real(8), dimension(:, :, :, :), pointer :: pressure_vof
        real(8), dimension(:, :, :, :), pointer :: density_vof
        real(8), dimension(:, :, :, :), pointer :: dp_drho_vof
        real(8), dimension(:, :, :, :), pointer :: dp_de_vof
        real(8), dimension(:, :, :, :), pointer :: sie_vof
        real(8), dimension(:, :, :, :), pointer :: mat_vof

        real(8), dimension(:, :, :), pointer :: mat_vol          
        real(8), dimension(:, :, :), pointer :: material_x       
        real(8), dimension(:, :, :), pointer :: material_y       
        real(8), dimension(:, :, :), pointer :: material_z       

        real(8) :: vol_diff                  
        real(8) :: a_visc_temp               
        real(8) :: vol_diff_stress           
        real(8) :: sie_vof_temp              
        real(8) :: pressure_temp             
        real(8) :: dp_drho_temp              
        real(8) :: dp_de_temp                
        real(8) :: mass_vof_temp             
        real(8) :: vol_vof_temp              
        real(8) :: vol_diff_vof_temp         
        real(8) :: vol_diff_vof_stress_temp  
        real(8) :: sie_diff                  
        real(8) :: stress_fac                

        integer :: i, j, k, tmp_mat          

        call this%velocity         %Point_to_data(velocity_x, velocity_y, velocity_z)
        call this%mesh             %Point_to_data(x, y, z)
        call this%total_temperature%Point_to_data(temperature)
        call this%total_pressure   %Point_to_data(pressure)
        call this%total_dp_drho    %Point_to_data(dp_drho)
        call this%total_dp_de      %Point_to_data(dp_de)
        call this%total_volume     %Point_to_data(vol)
        call this%total_cell_mass  %Point_to_data(cell_mass)
        call this%total_sie        %Point_to_data(sie)
        call this%total_vof        %Point_to_data(vof)
        call this%a_visc           %Point_to_data(a_visc)
        call this%num_mat_cells    %Point_to_data(nmats_in_cell)
        call this%rezone%Point_to_data (material_x, material_y, material_z)
        call this%rezone%material_volume%Point_to_data(mat_vol)



        call this%materials%temperature%Point_to_data(temperature_vof)
        call this%materials%cell_mass  %Point_to_data(cell_mass_vof)
        call this%materials%pressure   %Point_to_data(pressure_vof)
        call this%materials%dp_drho    %Point_to_data(dp_drho_vof)
        call this%materials%dp_de      %Point_to_data(dp_de_vof)
        call this%materials%sie        %Point_to_data(sie_vof)
        call this%materials%vof        %Point_to_data(mat_vof)

        do k = 1, this%nz
            do j = 1, this%ny
                do i = 1, this%nx
                    do tmp_mat = 1, this%nmats

                        if (vof(i, j, k) < this%emf) then
                            sie(i, j, k) = 0d0
                            temperature(i, j, k) = teps
                            temperature_vof(tmp_mat, i, j, k) = 0d0
                            cycle
                        end if

                        vol_diff = mat_vol(i, j, k) - vol(i, j, k)
                        vol_diff_stress = 0d0

                        if (vol_diff < 0d0) then
                            a_visc_temp = a_visc(i, j, k)
                        else
                            a_visc_temp = 0d0
                        end if

                        if (mat_vof(tmp_mat, i, j, k) < this%emf) then
                            sie_vof(tmp_mat, i, j, k) = 0d0
                            temperature_vof(tmp_mat, i, j, k) = 0d0
                        else
                            sie_vof_temp = sie_vof(tmp_mat, i, j, k)
                            if (nmats_in_cell(i, j, k) == 1d0) then
                                pressure_temp = pressure(i, j, k) 
                                dp_de_temp    = dp_de   (i, j, k)
                                dp_drho_temp  = dp_drho (i, j, k)
                            else
                                pressure_temp = pressure_vof(tmp_mat, i, j, k)
                                dp_de_temp    = dp_de_vof   (tmp_mat, i, j, k)
                                dp_drho_temp  = dp_drho_vof (tmp_mat, i, j, k)
                            end if
                            pressure_temp  = pressure_temp * vof(i, j, k)   
                            mass_vof_temp  = cell_mass_vof(tmp_mat, i, j, k)


                            vol_diff_vof_temp = vol_diff * mat_vof(tmp_mat, i, j, k)
                            vol_vof_temp  = vol(i, j, k) * mat_vof(tmp_mat, i, j, k)

                            stress_fac = 0d0

                            vol_diff_vof_stress_temp = vol_diff_stress * stress_fac * vof(i, j, k)

                            sie_diff = - ((pressure_temp - &
                                (0.5d0 * dp_drho_temp * mass_vof_temp / vol_vof_temp ** 2 * vol_diff_vof_temp) + a_visc_temp) *&
                                vol_diff_vof_temp + vol_diff_vof_stress_temp) / &
                                (1d0 + 0.5d0 * dp_de_temp * vol_diff_vof_temp / (mass_vof_temp + 1d-30)) / (mass_vof_temp+1d-30)

                            if (sie_diff == 0d0) cycle
                            sie_vof(tmp_mat, i, j, k) = sie_vof(tmp_mat, i, j, k) + sie_diff
                            sie(i, j, k) = sie(i, j, k) + sie_vof(tmp_mat, i, j, k) * mass_vof_temp
                        end if
                    end do
                end do
            end do
        end do

        call this%total_sie%Exchange_virtual_space_nonblocking()

        call this%materials%sie%Exchange_virtual_space_nonblocking()


        return
    end subroutine Calculate_energy_3d


    subroutine Fix_vof_2d(this, dt)
        use geometry_module, only : Quad_volume
        implicit none


        class (hydro_step_t), intent(in out) :: this
        real(8)             , intent(in)     :: dt

        real(8), dimension(:, :, :), pointer :: velocity_x   
        real(8), dimension(:, :, :), pointer :: velocity_y   
        real(8), dimension(:, :, :), pointer :: x            
        real(8), dimension(:, :, :), pointer :: y            

        real(8), dimension(:, :, :), pointer :: nmats_in_cell
        real(8), dimension(:, :, :), pointer :: cell_mass           
        real(8), dimension(:, :, :), pointer :: pressure            
        real(8), dimension(:, :, :), pointer :: a_visc              
        real(8), dimension(:, :, :), pointer :: vol                 
        real(8), dimension(:, :, :), pointer :: vof                 
        real(8), dimension(:, :, :), pointer :: sie                 

        real(8), dimension(:, :, :, :), pointer :: cell_mass_vof
        real(8), dimension(:, :, :, :), pointer :: mat_vof
        real(8), dimension(:, :, :, :), pointer :: sie_vof

        real(8) :: x1, x2, x3, x4       
        real(8) :: y1, y2, y3, y4       
        real(8) :: vol_new              
        real(8) :: vof_old              
        real(8) :: vol_ratio            
        real(8) :: vof_sum              
        real(8) :: vof_max              
        integer :: mat_vof_max          
        real(8) :: cell_mass_vof_sum    
        real(8) :: sie_vof_sum          
        real(8) :: sie_tot              

        integer :: i, j, tmp_mat     
        real(8) :: emf1                 

        call this%velocity         %Point_to_data(velocity_x, velocity_y)
        call this%mesh             %Point_to_data(x, y)

        call this%total_cell_mass  %Point_to_data(cell_mass)
        call this%total_pressure   %Point_to_data(pressure)
        call this%num_mat_cells    %Point_to_data(nmats_in_cell)
        call this%total_volume     %Point_to_data(vol)
        call this%total_sie        %Point_to_data(sie)
        call this%total_vof        %Point_to_data(vof)
        call this%a_visc           %Point_to_data(a_visc)
                    call this%materials%sie      %Point_to_data(sie_vof)
                    call this%materials%vof      %Point_to_data(mat_vof)
                    call this%materials%cell_mass%Point_to_data(cell_mass_vof)
        emf1 = 1 - this%emf

        do j = 1, this%ny
            do i = 1, this%nx



                if (vof(i, j, 1) > emf1) then  
                    vof(i, j, 1) = 1d0
                else if (vof(i, j, 1) < this%emf) then 
                    nmats_in_cell(i, j, 1) = 0
                    vof(i, j, 1) = 0d0
                else


                    x1 = x(i    , j    , 1) + velocity_x(i    , j    , 1) * dt
                    x2 = x(i + 1, j    , 1) + velocity_x(i + 1, j    , 1) * dt
                    x3 = x(i + 1, j + 1, 1) + velocity_x(i + 1, j + 1, 1) * dt
                    x4 = x(i    , j + 1, 1) + velocity_x(i    , j + 1, 1) * dt

                    y1 = y(i    , j    , 1) + velocity_y(i    , j    , 1) * dt
                    y2 = y(i + 1, j    , 1) + velocity_y(i + 1, j    , 1) * dt
                    y3 = y(i + 1, j + 1, 1) + velocity_y(i + 1, j + 1, 1) * dt
                    y4 = y(i    , j + 1, 1) + velocity_y(i    , j + 1, 1) * dt

                    vol_new   = Quad_volume(x1, y1, x2, y2, x3, y3, x4, y4, this%mesh%cyl)
                    vol_ratio = vol(i, j, 1) / vol_new
                    vof_old   = vof(i, j, 1)

                    if ((1d0 - 1d0 / vol_ratio) * (pressure(i, j, 1) + a_visc(i, j, 1)) > 0d0) then
                        if (vof(i, j, 1) < emf1) vof(i, j, 1) = vof(i, j, 1) * vol_ratio
                        if (vof(i, j, 1) > emf1) vof(i, j, 1) = 1d0
                        vol_ratio = vof(i, j, 1) / vof_old
                        do tmp_mat = 1, this%nmats

                            if (mat_vof(tmp_mat, i, j, 1) < emf1) mat_vof(tmp_mat, i, j, 1) = mat_vof(tmp_mat, i, j, 1) * vol_ratio
                        end do
                    end if
                end if

                vof_sum   = 0d0
                vof_max   = 0d0
                sie_vof_sum = 0d0
                cell_mass_vof_sum = 0d0

                do tmp_mat = 1, this%nmats

                    if (mat_vof(tmp_mat, i, j, 1) > emf1) then
                        mat_vof      (tmp_mat, i, j, 1) = 1d0
                    else if (mat_vof(tmp_mat, i, j, 1) < this%emf) then
                        mat_vof      (tmp_mat, i, j, 1) = 0d0
                        cell_mass_vof(tmp_mat, i, j, 1) = 0d0
                        sie_vof      (tmp_mat, i, j, 1) = 0d0
                    end if
                    vof_sum            = vof_sum            + mat_vof(tmp_mat, i, j, 1)
                    sie_vof_sum        = sie_vof_sum        + sie_vof(tmp_mat, i, j, 1) * cell_mass_vof(tmp_mat, i, j, 1)
                    cell_mass_vof_sum  = cell_mass_vof_sum  + cell_mass_vof (tmp_mat, i, j, 1)
                    if (vof_max <= mat_vof(tmp_mat, i, j, 1)) then
                        vof_max     = mat_vof(tmp_mat, i, j, 1)
                        mat_vof_max = tmp_mat
                    end if
                end do

                if (vof_sum /= vof(i, j, 1)) then

                    mat_vof(mat_vof_max, i, j, 1) = mat_vof(mat_vof_max, i, j, 1) + vof(i, j, 1) - vof_sum
                end if


                sie_tot = sie(i, j, 1) * cell_mass(i, j, 1)
                if ((abs(sie_vof_sum - sie_tot) > 1d-4 * abs(sie_tot)) .and. (abs(sie_tot) > 1d-40)) then
                    sie(i, j, 1) = sie_vof_sum / (cell_mass(i, j, 1) + 1d-20)
                end if
            end do
        end do
        return
    end subroutine Fix_vof_2d


    subroutine Fix_vof_3d(this, dt, from_interp_mesh)
        implicit none

        class (hydro_step_t), intent(in out) :: this
        real(8)             , intent(in)     :: dt
        integer             , intent(in)     :: from_interp_mesh  

        real(8), dimension(:, :, :), pointer :: mat_vol  

        real(8), dimension(:, :, :), pointer :: nmats_in_cell
        real(8), dimension(:, :, :), pointer :: cell_mass           
        real(8), dimension(:, :, :), pointer :: pressure            
        real(8), dimension(:, :, :), pointer :: a_visc              
        real(8), dimension(:, :, :), pointer :: vol                 
        real(8), dimension(:, :, :), pointer :: vof                 
        real(8), dimension(:, :, :), pointer :: sie                 

        real(8), dimension(:, :, :, :), pointer :: cell_mass_vof
        real(8), dimension(:, :, :, :), pointer :: mat_vof
        real(8), dimension(:, :, :, :), pointer :: sie_vof

        real(8), dimension(:, :, :), allocatable :: vof_sum_arr          
        real(8), dimension(:, :, :), allocatable :: vof_max_arr          
        real(8), dimension(:, :, :), allocatable :: sie_vof_sum_arr          
        real(8), dimension(:, :, :), allocatable :: cell_mass_vof_sum_arr          
        integer, dimension(:, :, :), allocatable :: mat_vof_max_arr          

        real(8) :: vol_new              
        real(8) :: vof_old              
        real(8) :: vol_ratio            
        real(8) :: vof_sum              
        real(8) :: vof_max              
        integer :: mat_vof_max          
        real(8) :: cell_mass_vof_sum    
        real(8) :: sie_vof_sum          
        real(8) :: sie_tot              

        integer :: i, j, k, tmp_mat     
        real(8) :: emf1                 

        call this%total_cell_mass  %Point_to_data(cell_mass)
        call this%total_pressure   %Point_to_data(pressure)
        call this%num_mat_cells    %Point_to_data(nmats_in_cell)
        call this%total_volume     %Point_to_data(vol)
        call this%total_sie        %Point_to_data(sie)
        call this%total_vof        %Point_to_data(vof)
        call this%a_visc           %Point_to_data(a_visc)
        call this%rezone%material_volume%Point_to_data(mat_vol)

        emf1 = 1 - this%emf

        allocate(vof_sum_arr(1:this%nx, 1:this%ny, 1:this%nz))
        allocate(vof_max_arr(1:this%nx, 1:this%ny, 1:this%nz))
        allocate(sie_vof_sum_arr(1:this%nx, 1:this%ny, 1:this%nz))
        allocate(cell_mass_vof_sum_arr(1:this%nx, 1:this%ny, 1:this%nz))
        allocate(mat_vof_max_arr(1:this%nx, 1:this%ny, 1:this%nz))

        vof_sum_arr = 0d0
        vof_max_arr = 0d0
        sie_vof_sum_arr = 0d0
        cell_mass_vof_sum_arr = 0d0
        mat_vof_max_arr = 0

            call this%materials%sie      %Point_to_data(sie_vof)
            call this%materials%vof      %Point_to_data(mat_vof)
            call this%materials%cell_mass%Point_to_data(cell_mass_vof)
!        do tmp_mat = 1, this%nmats
!
!
!            do k = 1, this%nz
!                do j = 1, this%ny
!                    do i = 1, this%nx
!
!                        if (vof(i, j, k) > emf1) then
!                            vof(i, j, k) = 1d0
!                        else if (vof(i, j, k) < this%emf) then
!                            nmats_in_cell(i, j, k) = 0
!                            vof(i, j, k) = 0d0
!                        else if (from_interp_mesh /= 1) then
!
!                            vol_new   = mat_vol(i, j, k)
!                            vol_ratio = vol(i, j, k) / vol_new
!                            vof_old   = vof(i, j, k)
!
!                            if ((1d0 - 1d0 / vol_ratio) * (pressure(i, j, k) + a_visc(i, j, k)) > 0d0) then
!                                if (vof(i, j, k) < emf1) vof(i, j, k) = vof(i, j, k) * vol_ratio
!                                if (vof(i, j, k) > emf1) vof(i, j, k) = 1d0
!                                vol_ratio = vof(i, j, k) / vof_old
!                                if (mat_vof(i, j, k) < emf1) mat_vof(i, j, k) = mat_vof(i, j, k) * vol_ratio
!                            end if
!                        end if
!
!
!                        if (mat_vof(i, j, k) > emf1) then
!                            mat_vof      (i, j, k) = 1d0
!                        else if (mat_vof(i, j, k) < this%emf) then
!                            mat_vof      (i, j, k) = 0d0
!                            cell_mass_vof(i, j, k) = 0d0
!                            sie_vof      (i, j, k) = 0d0
!                        end if
!                        vof_sum_arr(i, j, k)            = vof_sum_arr(i, j, k)            + mat_vof(i, j, k)
!                        sie_vof_sum_arr(i, j, k)        = sie_vof_sum_arr(i, j, k)        + sie_vof(i, j, k) * cell_mass_vof(i, j, k)
!                        cell_mass_vof_sum_arr(i, j, k)  = cell_mass_vof_sum_arr(i, j, k)  + cell_mass_vof (i, j, k)
!                        if (vof_max_arr(i, j, k) <= mat_vof(i, j, k)) then
!                            vof_max_arr(i, j, k) = mat_vof(i, j, k)
!                            mat_vof_max_arr(i, j, k) = tmp_mat
!                        end if
!
!                    end do
!                end do
!            end do
!        end do
!
!        do k = 1, this%nz
!            do j = 1, this%ny
!                do i = 1, this%nx
!                    if (vof_sum_arr(i,j,k) /= vof(i,j,k)) then
!                        call this%materials(mat_vof_max_arr(i,j,k))%vof%Point_to_data(mat_vof)
!                        mat_vof(i, j, k) = mat_vof(i, j, k) + vof(i, j, k) - vof_sum_arr(i,j,k)
!                    end if
!
!                    sie_tot = sie(i, j, k) * cell_mass(i, j, k)
!                    if ((abs(sie_vof_sum_arr(i,j,k) - sie_tot) > 1d-4 * abs(sie_tot)) .and. (abs(sie_tot) > 1d-40)) then
!                        sie(i, j, k) = sie_vof_sum_arr(i,j,k) / (cell_mass(i, j, k) + 1d-20)
!                    end if
!                end do
!            end do
!        end do
call this%velocity%Calculate_derivatives(this%mesh%coordinates, &
            this%total_vof, this%total_volume, this%nx, this%ny, this%nz, this%emf)
call this%total_sie%exchange_end()
call this%materials%sie%exchange_end()

        deallocate(vof_sum_arr)
        deallocate(vof_max_arr)
        deallocate(sie_vof_sum_arr)
        deallocate(cell_mass_vof_sum_arr)
        deallocate(mat_vof_max_arr)
        call this%total_vof%Apply_boundarY(.false.)
            call this%materials%vof%Apply_boundarY(is_blocking=.false.)

        return
    end subroutine Fix_vof_3d


    subroutine Calculate_mesh_2d(this, time)
        implicit none
        class (hydro_step_t), intent(in out) :: this
        class (time_t)      , intent(in out) :: time

        real(8), dimension(:, :, :), pointer :: x               
        real(8), dimension(:, :, :), pointer :: y               
        real(8), dimension(:, :, :), pointer :: material_x      
        real(8), dimension(:, :, :), pointer :: material_y      
        real(8), dimension(:, :, :), pointer :: vol             
        real(8), dimension(:, :, :), pointer :: r_factor        
        real(8), dimension(:, :, :), pointer :: velocity_x      
        real(8), dimension(:, :, :), pointer :: velocity_y      
        real(8), dimension(:, :, :), pointer :: mesh_velocity_x 
        real(8), dimension(:, :, :), pointer :: mesh_velocity_y 
        real(8), dimension(:, :, :), pointer :: mat_vol         
        real(8), dimension(:, :, :), pointer :: momentum_x      
        real(8), dimension(:, :, :), pointer :: momentum_y      
        real(8), dimension(:, :, :), pointer :: adv_momentum_x  
        real(8), dimension(:, :, :), pointer :: adv_momentum_y  
        real(8), dimension(:, :, :), pointer :: vertex_mass           

        real(8), dimension(4) :: dx  
        real(8), dimension(4) :: dy  
        real(8), dimension(4) :: inv_dist_sq  

        real(8) :: vel_diff_x   
        real(8) :: vel_diff_y   
        real(8) :: vel_diff     
        real(8) :: vel_grad     
        real(8) :: vel_rez_max  

        integer :: i, j, ii 

        call this%velocity         %Point_to_data(velocity_x, velocity_y)
        call this%mesh             %Point_to_r_factor(r_factor)
        call this%mesh             %Point_to_data(x, y)
        call this%total_volume     %Point_to_data(vol)
        call this%vertex_mass%Point_to_data(vertex_mass)

        call this%rezone%Point_to_coordinates_2d(material_x, material_y)
        call this%rezone%mesh_velocity  %Point_to_data(mesh_velocity_x, mesh_velocity_y)
        call this%advect%momentum       %Point_to_data(momentum_x     , momentum_y)
        call this%advect%adv_momentum   %Point_to_data(adv_momentum_x , adv_momentum_y)
        call this%rezone%material_volume%Point_to_data(mat_vol)



        do j = 1, this%nyp
            do i = 1, this%nxp
                x(i, j, 1)        = x(i, j, 1) + time%dt * mesh_velocity_x(i, j, 1)
                y(i, j, 1)        = y(i, j, 1) + time%dt * mesh_velocity_y(i, j, 1)
                r_factor(i, j, 1) = x(i, j, 1) * this%mesh%cyl + this%mesh%omcyl
            end do
        end do



















        adv_momentum_x = 0d0
        adv_momentum_y = 0d0


        do j = 1, this%ny
            do i = 1, this%nx
                momentum_x(i, j, 1) = 0.25d0 * (velocity_x(i+1, j  , 1) + velocity_x(i+1, j+1, 1) + &
                    velocity_x(i  , j+1, 1) + velocity_x(i  , j  , 1))
                momentum_y(i, j, 1) = 0.25d0 * (velocity_y(i+1, j  , 1) + velocity_y(i+1, j+1, 1) + &
                    velocity_y(i  , j+1, 1) + velocity_y(i  , j  , 1))
            end do
        end do







        call this%rezone%material_volume%Calculate(this%rezone%material_coordinates, cyl_optional=this%mesh%cyl)

        call this%Calculate_density(this%rezone%material_volume)





        call this%total_volume%calculate(this%mesh%coordinates, cyl_optional=this%mesh%cyl)
        call this%total_volume%Apply_boundary()



    end subroutine Calculate_mesh_2d


    subroutine Calculate_mesh_3d(this, time)
        implicit none
        class (hydro_step_t), intent(in out) :: this
        class (time_t)      , intent(in out) :: time


        real(8), dimension(:, :, :), pointer :: x               
        real(8), dimension(:, :, :), pointer :: y               
        real(8), dimension(:, :, :), pointer :: z               
        real(8), dimension(:, :, :), pointer :: velocity_x      
        real(8), dimension(:, :, :), pointer :: velocity_y      
        real(8), dimension(:, :, :), pointer :: velocity_z      
        real(8), dimension(:, :, :), pointer :: mesh_velocity_x 
        real(8), dimension(:, :, :), pointer :: mesh_velocity_y 
        real(8), dimension(:, :, :), pointer :: mesh_velocity_z 
        real(8), dimension(:, :, :), pointer :: vertex_mass     
        real(8), dimension(:, :, :), pointer :: vol             
        real(8), dimension(:, :, :), pointer :: mat_vol             

        real(8) :: vel_diff_x   
        real(8) :: vel_diff_y   
        real(8) :: vel_diff_z   
        real(8) :: vel_grad     
        real(8) :: vel_rez_max  

        real(8), parameter :: eps = 1d-20  

        real(8) :: vel_grad_i_p  
        real(8) :: vel_grad_i_m  
        real(8) :: vel_grad_j_p  
        real(8) :: vel_grad_j_m  
        real(8) :: vel_grad_k_p  
        real(8) :: vel_grad_k_m  

        real(8) :: len_sq        

        integer :: i, j, k 
        logical :: is_wall_x_top, is_wall_x_bot,is_wall_y_top, is_wall_y_bot,is_wall_z_top,is_wall_z_bot

        is_wall_x_top = this%parallel_params%is_wall_x_top
        is_wall_x_bot = this%parallel_params%is_wall_x_bot

        is_wall_y_top = this%parallel_params%is_wall_y_top
        is_wall_y_bot = this%parallel_params%is_wall_y_bot

        is_wall_z_top = this%parallel_params%is_wall_z_top
        is_wall_z_bot = this%parallel_params%is_wall_z_bot

        call this%velocity         %Point_to_data(velocity_x, velocity_y, velocity_z)
        call this%mesh             %Point_to_data(x, y, z)
        call this%vertex_mass%Point_to_data(vertex_mass)
        call this%total_volume     %Point_to_data(vol)
        call this%rezone%mesh_velocity  %Point_to_data(mesh_velocity_x, mesh_velocity_y, mesh_velocity_z)









        call this%mesh%coordinates%Calculate(time%dt, this%rezone%mesh_velocity)
        call this%mesh%Exchange_virtual_space_blocking()









        call this%mesh%coordinates%Apply_boundary(this%mesh%coordinates%data)
        call this%rezone%material_coordinates%Apply_boundary(this%mesh%coordinates%data)







        call this%total_vof%Exchange_end()
            call this%materials%vof%Exchange_end()
        call this%Calculate_density(this%rezone%material_volume)
        call this%total_density%Apply_boundary()
!            call this%materials%density  %Exchange_virtual_space_blocking()






        call this%total_volume%Calculate(this%mesh%coordinates)
        call this%total_volume%Apply_boundary()


        return
    end subroutine Calculate_mesh_3d


    subroutine Boundary_mirror_image(this)
        use geometry_module, only: Mirror_image_3d
        implicit none
        class (hydro_step_t), intent(in out) :: this

        return
    end subroutine Boundary_mirror_image





    subroutine Apply_wilkins_mass(this)
        implicit none

        class (hydro_step_t), intent(in out) :: this

        real(8), dimension(:, :, :), pointer :: vertex_mass, prev_vertex_mass, cell_mass, prev_cell_mass


        call this%vertex_mass%Point_to_data(vertex_mass)
        call this%previous_vertex_mass%Point_to_data(prev_vertex_mass)
        call this%total_cell_mass%Point_to_data(cell_mass)
        call this%previous_cell_mass%Point_to_data(prev_cell_mass)

        prev_cell_mass  (1:this%nxp, 1:this%nyp, 1) = cell_mass  (1:this%nxp, 1:this%nyp, 1)
        prev_vertex_mass(1:this%nxp, 1:this%nyp, 1) = vertex_mass(1:this%nxp, 1:this%nyp, 1)


        call this%vertex_mass%Calculate_vertex_mass_2d(this%mesh%coordinates, this%total_density, this%total_cell_mass,&
            this%wilkins_scheme, 0d0) 
        call this%vertex_mass%Exchange_virtual_space_blocking()
        call this%Calculate_inversed_vertex_mass()
        call this%inversed_vertex_mass%Exchange_virtual_space_blocking()

        return


    end subroutine Apply_wilkins_mass

    subroutine Restore_wilkins_mass(this)
        implicit none

        class (hydro_step_t), intent(in out) :: this

        real(8), dimension(:, :, :), pointer :: vertex_mass, prev_vertex_mass, cell_mass, prev_cell_mass


        call this%vertex_mass%Point_to_data(vertex_mass)
        call this%previous_vertex_mass%Point_to_data(prev_vertex_mass)
        call this%total_cell_mass%Point_to_data(cell_mass)
        call this%previous_cell_mass%Point_to_data(prev_cell_mass)


        cell_mass  (1:this%nxp, 1:this%nyp, 1) = prev_cell_mass  (1:this%nxp, 1:this%nyp, 1)
        vertex_mass(1:this%nxp, 1:this%nyp, 1) = prev_vertex_mass(1:this%nxp, 1:this%nyp, 1)
        return

    end subroutine Restore_wilkins_mass



    subroutine Calculate_stresd(this)
        class (hydro_step_t), intent(in out) :: this

    end subroutine Calculate_stresd


    subroutine Set_communication(this, comm, comm_params_cell, comm_params_vertex, comm_material)
        class (hydro_step_t)            :: this 
        type(communication_t), pointer            :: comm
        type(communication_parameters_t), pointer :: comm_params_cell, comm_params_vertex, comm_material



        call this%advect%Set_communication(comm, comm_params_cell, comm_params_vertex, comm_material)
        call this%rezone%Set_communication(comm, comm_params_cell, comm_params_vertex)
    end subroutine Set_communication



    subroutine Point_to_velocity_data_2d (this, ptr_x, ptr_y)
        class (hydro_step_t)              , intent(in out) :: this         
        real(8), dimension(:,:,:), pointer, intent(out)    :: ptr_x, ptr_y 

        call this%velocity%Point_to_data(ptr_x, ptr_y)
    end subroutine Point_to_velocity_data_2d

    subroutine Point_to_velocity_data_3d (this, ptr_x, ptr_y, ptr_z)
        class (hydro_step_t)              , intent(in out) :: this                
        real(8), dimension(:,:,:), pointer, intent(out)    :: ptr_x, ptr_y, ptr_z 

        call this%velocity%Point_to_data(ptr_x, ptr_y, ptr_z)
    end subroutine Point_to_velocity_data_3d

    subroutine Point_to_mesh_data_2d (this, ptr_x, ptr_y)
        class (hydro_step_t)              , intent(in out) :: this         
        real(8), dimension(:,:,:), pointer, intent(out)    :: ptr_x, ptr_y 

        call this%mesh%Point_to_data(ptr_x, ptr_y)
    end subroutine Point_to_mesh_data_2d

    subroutine Point_to_mesh_data_3d (this, ptr_x, ptr_y, ptr_z)
        class (hydro_step_t)              , intent(in out) :: this         
        real(8), dimension(:,:,:), pointer, intent(out)    :: ptr_x, ptr_y, ptr_z 

        call this%mesh%Point_to_data(ptr_x, ptr_y,ptr_z)
    end subroutine Point_to_mesh_data_3d

    subroutine Point_to_r_factor_data (this, ptr)
        class (hydro_step_t)              , intent(in out) :: this  
        real(8), dimension(:,:,:), pointer, intent(out)    :: ptr   

        call this%mesh%Point_to_r_factor(ptr)
    end subroutine Point_to_r_factor_data

    subroutine Point_to_cell_mass_data (this, ptr)
        class (hydro_step_t)              , intent(in out) :: this  
        real(8), dimension(:,:,:), pointer, intent(out)    :: ptr   

        call this%total_cell_mass%Point_to_data(ptr)
    end subroutine Point_to_cell_mass_data

    subroutine Point_to_temperature_data (this, ptr)
        class (hydro_step_t)              , intent(in out) :: this  
        real(8), dimension(:,:,:), pointer, intent(out)    :: ptr   

        call this%total_temperature%Point_to_data(ptr)
    end subroutine Point_to_temperature_data

    subroutine Point_to_pressure_data (this, ptr)
        class (hydro_step_t)              , intent(in out) :: this  
        real(8), dimension(:,:,:), pointer, intent(out)    :: ptr   

        call this%total_pressure%Point_to_data(ptr)
    end subroutine Point_to_pressure_data

    subroutine Point_to_sie_data (this, ptr)
        class (hydro_step_t)              , intent(in out) :: this  
        real(8), dimension(:,:,:), pointer, intent(out)    :: ptr   

        call this%total_sie%Point_to_data(ptr)
    end subroutine Point_to_sie_data

    subroutine Point_to_sound_velocity_data (this, ptr)
        class (hydro_step_t)              , intent(in out) :: this  
        real(8), dimension(:,:,:), pointer, intent(out)    :: ptr   

        call this%total_sound_vel%Point_to_data(ptr)
    end subroutine Point_to_sound_velocity_data


    subroutine Point_to_vol_data (this, ptr)
        class (hydro_step_t)              , intent(in out) :: this  
        real(8), dimension(:,:,:), pointer, intent(out)    :: ptr   

        call this%total_volume%Point_to_data(ptr)
    end subroutine Point_to_vol_data

    subroutine Point_to_density_data (this, ptr)
        class (hydro_step_t)              , intent(in out) :: this  
        real(8), dimension(:,:,:), pointer, intent(out)    :: ptr   

        call this%total_density%Point_to_data(ptr)
    end subroutine Point_to_density_data


    subroutine Point_to_artificial_viscosity_data (this, ptr)
        class (hydro_step_t)              , intent(in out) :: this  
        real(8), dimension(:,:,:), pointer, intent(out)    :: ptr   

        call this%a_visc%Point_to_data(ptr)
    end subroutine Point_to_artificial_viscosity_data


    subroutine Point_to_vof_data (this, ptr)
        class (hydro_step_t)              , intent(in out) :: this  
        real(8), dimension(:,:,:), pointer, intent(out)    :: ptr   

        call this%total_vof%Point_to_data(ptr)
    end subroutine Point_to_vof_data


    subroutine Point_to_density_vof_data (this, ptr, material_num)
        class (hydro_step_t)              , intent(in out) :: this  
        real(8), dimension(:,:,:), pointer, intent(out)    :: ptr   
        integer                           , intent(in)     :: material_num

!        call this%materials(material_num)%density%Point_to_data(ptr)
    end subroutine Point_to_density_vof_data

    subroutine Point_to_cell_mass_vof_data (this, ptr, material_num)
        class (hydro_step_t)              , intent(in out) :: this  
        real(8), dimension(:,:,:), pointer, intent(out)    :: ptr   
        integer                           , intent(in)     :: material_num

!        call this%materials(material_num)%cell_mass%Point_to_data(ptr)
    end subroutine Point_to_cell_mass_vof_data

    subroutine Point_to_mat_vof_data (this, ptr, material_num)
        class (hydro_step_t)              , intent(in out) :: this  
        real(8), dimension(:,:,:), pointer, intent(out)    :: ptr   
        integer                           , intent(in)     :: material_num

!        call this%materials(material_num)%vof%Point_to_data(ptr)
    end subroutine Point_to_mat_vof_data

    subroutine Point_to_temperature_vof_data (this, ptr, material_num)
        class (hydro_step_t)              , intent(in out) :: this  
        real(8), dimension(:,:,:), pointer, intent(out)    :: ptr   
        integer                                            :: material_num

!        call this%materials(material_num)%temperature%Point_to_data(ptr)
    end subroutine Point_to_temperature_vof_data

    subroutine Point_to_temperature_old_vof_data (this, ptr, material_num)
        class (hydro_step_t)              , intent(in out) :: this  
        real(8), dimension(:,:,:), pointer, intent(out)    :: ptr   
        integer                           , intent(in)     :: material_num

!        call this%materials(material_num)%temperature_old%Point_to_data(ptr)
    end subroutine Point_to_temperature_old_vof_data

    subroutine Point_to_pressure_vof_data (this, ptr, material_num)
        class (hydro_step_t)              , intent(in out) :: this  
        real(8), dimension(:,:,:), pointer, intent(out)    :: ptr   
        integer                                            :: material_num

!        call this%materials(material_num)%pressure%Point_to_data(ptr)
    end subroutine Point_to_pressure_vof_data

    subroutine Point_to_sie_vof_data (this, ptr, material_num)
        class (hydro_step_t)              , intent(in out) :: this  
        real(8), dimension(:,:,:), pointer, intent(out)    :: ptr   
        integer                           , intent(in)     :: material_num

!        call this%materials(material_num)%sie%Point_to_data(ptr)
    end subroutine Point_to_sie_vof_data

    subroutine Point_to_sound_velocity_vof_data (this, ptr, material_num)
        class (hydro_step_t)              , intent(in out) :: this  
        real(8), dimension(:,:,:), pointer, intent(out)    :: ptr   
        integer                           , intent(in)     :: material_num

!        call this%materials(material_num)%sound_vel%Point_to_data(ptr)
    end subroutine Point_to_sound_velocity_vof_data

    subroutine Point_to_deriv_vof_data (this, dp_de_vof, dp_drho_vof, dt_de_vof, dt_drho_vof, material_num)
        class (hydro_step_t)              , intent(in out) :: this  
        real(8), dimension(:,:,:), pointer, intent(out)    :: dp_de_vof, dp_drho_vof 
        real(8), dimension(:,:,:), pointer, intent(out)    :: dt_de_vof, dt_drho_vof 
        integer                           , intent(in)     :: material_num

!        call this%materials(material_num)%dp_de%Point_to_data(dp_de_vof)
!        call this%materials(material_num)%dp_drho%Point_to_data(dp_drho_vof)
!        call this%materials(material_num)%dt_de%Point_to_data(dt_de_vof)
!        call this%materials(material_num)%dt_drho%Point_to_data(dt_drho_vof)
    end subroutine Point_to_deriv_vof_data

    subroutine Point_to_mesh_velocity_data (this, ptr_x, ptr_y)
        class (hydro_step_t)              , intent(in out) :: this         
        real(8), dimension(:,:,:), pointer, intent(out)    :: ptr_x, ptr_y 

        call this%rezone%Point_to_velocities(ptr_x, ptr_y)
    end subroutine Point_to_mesh_velocity_data

    subroutine Point_to_advect (this, advect_ptr)
        class (hydro_step_t), intent(in out) :: this  
        type(advect_t), pointer, intent(out) :: advect_ptr

        advect_ptr => this%advect

    end subroutine Point_to_advect


    subroutine Point_to_rezone (this, rezone_ptr)
        class (hydro_step_t), intent(in out) :: this  
        type(rezone_t), pointer, intent(out) :: rezone_ptr

        rezone_ptr => this%rezone

    end subroutine Point_to_rezone

    subroutine Point_to_mat_id_data (this, mat_id_ptr)
        class (hydro_step_t), intent(in out) :: this  
        real(8), dimension(:,:,:), pointer, intent(out)    :: mat_id_ptr     

        call this%mat_id%Point_to_data(mat_id_ptr)

    end subroutine Point_to_mat_id_data

    subroutine Point_to_num_mat_cells_data (this, ptr)
        class (hydro_step_t)              , intent(in out) :: this    
        real(8), dimension(:,:,:), pointer, intent(out)    :: ptr     

        call this%num_mat_cells%Point_to_data(ptr)
    end subroutine Point_to_num_mat_cells_data

    subroutine Point_to_vertex_mass_data (this, ptr)
        class (hydro_step_t)              , intent(in out) :: this    
        real(8), dimension(:,:,:), pointer, intent(out)    :: ptr     

        call this%vertex_mass%Point_to_data(ptr)
    end subroutine Point_to_vertex_mass_data

    subroutine Point_to_inversed_vertex_mass_data (this, ptr)
        class (hydro_step_t)              , intent(in out) :: this    
        real(8), dimension(:,:,:), pointer, intent(out)    :: ptr     

        call this%inversed_vertex_mass%Point_to_data(ptr)
    end subroutine Point_to_inversed_vertex_mass_data


    subroutine Point_to_acceleration_data_2d (this, ptr_x, ptr_y)
        class (hydro_step_t)              , intent(in out) :: this         
        real(8), dimension(:,:,:), pointer, intent(out)    :: ptr_x, ptr_y 

        call this%acceleration%Point_to_data(ptr_x, ptr_y)
    end subroutine Point_to_acceleration_data_2d


    subroutine debug_check_nan(this)
        class (hydro_step_t)              , intent(in out) :: this         
        integer :: i,j, size_d
        character*80 :: tmp
        character*80 :: tmp2

        write(tmp, '(I1)') this%parallel_params%my_rank


!        call this%num_mat_cells%debug_check_nan("Hydro: Number material Cells. Rank: " //trim(tmp))
!        call this%a_visc%debug_check_nan("Hydro: Artificial Viscosity. Rank: " //trim(tmp))
!        call this%mat_id%debug_check_nan("Hydro: Material Index. Rank: " //trim(tmp))
!        call this%total_sound_vel%debug_check_nan("Hydro: Sound Velocity. Rank: " //trim(tmp))
!        call this%acceleration%debug_check_nan("Hydro: Acceleration. Rank: " //trim(tmp))
!        call this%velocity%debug_check_nan("Hydro: Velocity. Rank: " //trim(tmp))
!        call this%total_temperature%debug_check_nan("Hydro: Temperature. Rank: " //trim(tmp))
!        call this%vertex_mass%debug_check_nan("Hydro: Vertex Mass. Rank: " //trim(tmp))
!        call this%mesh%coordinates%debug_check_nan("Hydro: Mesh coordinates. Rank: " //trim(tmp))
!        call this%total_cell_mass%debug_check_nan("Hydro: Cell Mass. Rank: " //trim(tmp))
!        call this%total_pressure%debug_check_nan("Hydro: Pressure. Rank: " //trim(tmp))
!        call this%total_density%debug_check_nan("Hydro: Density. Rank: " //trim(tmp))
!        call this%total_volume%debug_check_nan("Hydro: Volume. Rank: " //trim(tmp))
!        call this%total_sie%debug_check_nan("Hydro: SIE. Rank: " //trim(tmp))
!        call this%inversed_vertex_mass%debug_check_nan("Hydro: Inversed vertex mass. Rank: " //trim(tmp))
!        call this%total_vof%debug_check_nan("Hydro: Vof. Rank: " //trim(tmp))
!
!        do i = 1, this%nmats
!            call this%materials(i)%cell_mass%debug_check_nan('Hydro: material cell mass. Rank: ' //trim(tmp))
!            call this%materials(i)%sie%debug_check_nan('Hydro: material sie. Rank: ' //trim(tmp))
!            call this%materials(i)%vof%debug_check_nan('Hydro: material vof. Rank: ' //trim(tmp))
!        end do

    end subroutine debug_check_nan
   
   
    subroutine Write_hydro_step(this, unit, iostat, iomsg)
        class (hydro_step_t), intent(in) :: this  
        integer,      intent(in)     :: unit
        integer,      intent(out)    :: iostat
        character(*), intent(in out)  :: iomsg
        integer :: i

#ifdef DEBUG
        write(*,*) "@@@ in Write_hydro_step @@@"
#endif
!        call this%num_mat_cells%Write_quantity_abstract(unit, iostat=iostat, iomsg=iomsg)
!        call this%mat_id%Write_quantity_abstract(unit, iostat=iostat, iomsg=iomsg)
!        call this%mesh%Write_mesh_abstract(unit, iostat=iostat, iomsg=iomsg)
!        call this%total_cell_mass%Write_quantity_abstract(unit, iostat=iostat, iomsg=iomsg)
!        call this%velocity%Write_quantity_abstract(unit, iostat=iostat, iomsg=iomsg)
!        call this%total_density%Write_quantity_abstract(unit, iostat=iostat, iomsg=iomsg)
!        call this%total_volume%Write_quantity_abstract(unit, iostat=iostat, iomsg=iomsg)
!        call this%total_sie%Write_quantity_abstract(unit, iostat=iostat, iomsg=iomsg)
!        call this%total_vof%Write_quantity_abstract(unit, iostat=iostat, iomsg=iomsg)
!        call this%total_pressure_sum%Write_quantity_abstract(unit, iostat=iostat, iomsg=iomsg)
!
!        do i=1, size(this%materials)
!            call  this%materials(i)%Write_material_abstract(unit, iostat=iostat, iomsg=iomsg)
!        end do

#ifdef DEBUG
        write(*,*) "@@@ end Write_hydro_step @@@"
#endif

    end subroutine Write_hydro_step

    subroutine Read_hydro_step(this, unit, iostat, iomsg)
        class (hydro_step_t), intent(in out) :: this  
        integer,      intent(in)     :: unit
        integer,      intent(out)    :: iostat
        character(*), intent(in out)  :: iomsg
        integer :: i

#ifdef DEBUG
        write(*,*) "@@@ in Read_hydro_step @@@"
#endif
!        call this%num_mat_cells%Read_quantity_abstract(unit, iostat=iostat, iomsg=iomsg)
!        call this%mat_id%Read_quantity_abstract(unit, iostat=iostat, iomsg=iomsg)
!        call this%mesh%Read_mesh_abstract(unit, iostat=iostat, iomsg=iomsg)
!        call this%total_cell_mass%Read_quantity_abstract(unit, iostat=iostat, iomsg=iomsg)
!        call this%velocity%Read_quantity_abstract(unit, iostat=iostat, iomsg=iomsg)
!        call this%total_density%Read_quantity_abstract(unit, iostat=iostat, iomsg=iomsg)
!        call this%total_volume%Read_quantity_abstract(unit, iostat=iostat, iomsg=iomsg)
!        call this%total_sie%Read_quantity_abstract(unit, iostat=iostat, iomsg=iomsg)
!        call this%total_vof%Read_quantity_abstract(unit, iostat=iostat, iomsg=iomsg)
!        call this%total_pressure_sum%Read_quantity_abstract(unit, iostat=iostat, iomsg=iomsg)
!
!        do i=1, size(this%materials)
!!            call  this%materials(i)%Read_material_abstract(unit, iostat=iostat, iomsg=iomsg)
!        end do



#ifdef DEBUG
        write(*,*) "@@@ end Read_hydro_step @@@"
#endif

    end subroutine Read_hydro_step
   
   
   
end module hydro_step_module
