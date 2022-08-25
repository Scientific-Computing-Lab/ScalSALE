

module problem_module
    use vertex_boundary_condition_module, only : vertex_boundary_condition_t, vertex_bc_wrapper_t
    use cell_boundary_condition_module  , only : cell_boundary_condition_t, cell_bc_wrapper_t
    use artificial_viscosity_module     , only : artificial_viscosity_t
    use equation_of_state_module        , only : equation_of_state_t, eos_wrapper_t
    use sound_velocity_module           , only : sound_velocity_t
    use acceleration_module             , only : acceleration_t
    use materials_in_cells_module       , only : materials_in_cells_t
    use temperature_module              , only : temperature_t
    use ideal_gas_module                , only : ideal_gas_t
    use vertex_mass_module              , only : vertex_mass_t
    use slip_vertex_module              , only : slip_vertex_t
    use slip_vertex_3d_module           , only : slip_vertex_3d_t
    use lagrange_surface_cell_module    , only : lagrange_surface_cell_t
    use no_slip_vertex_3d_module        , only : no_slip_vertex_3d_t
    use boundary_mirror_2d_module       , only : boundary_mirror_2d_t
    use boundary_mirror_3d_module       , only : boundary_mirror_3d_t
    use lagrange_surface_vertex_module  , only : lagrange_surface_vertex_t
    use vof_module                      , only : vof_t
    use hydro_step_module               , only : hydro_step_t
    use cell_mass_module                , only : cell_mass_t
    use slip_cell_module                , only : slip_cell_t
    use slip_cell_3d_module             , only : slip_cell_3d_t
    use num_materials_in_cells_module   , only : num_materials_in_cells_t
    use pressure_module                 , only : pressure_t
    use velocity_module                 , only : velocity_t
    use material_module                 , only : material_t
    use datafile_module                 , only : datafile_t
    use density_module                  , only : density_t
    use energy_module                   , only : energy_t
    use volume_module                   , only : volume_t
    use data_module                     , only : data_t
    use mesh_base_module                , only : mesh_base_t!, mesh_t, mesh_3d_t
    use mesh_3d_module                  , only : mesh_3d_t
    use mesh_module                     , only : mesh_t
    use time_module                     , only : time_t
    use textual_diagnostic_module       , only : textual_diagnostic_t
    !    use silo_diagnostic_module       , only : silo_diagnostic_t
    !    use diagnostic_module               , only : diagnostic_t, textual_diagnostic_t, plot_diagnostic_t!, silo_diagnostic_t, plot_diagnostic_t, &
    !        textual_diagnostic_hdf5_t, Static_hdf5_init, Static_hdf5_close
    use communication_parameters_module , only : communication_parameters_t
    use communication_module            , only : communication_t
    use parallel_parameters_module      , only : parallel_parameters_t
    use boundary_parameters_module      , only : boundary_parameters_t
    use data_4d_module, only: data_4d_t
    use material_quantity_module    , only : material_quantity_t
    !   use hdf5
    !   use cr_module                       , only : cr_t
    use mpi
    implicit none
    private

    type, public :: problem_t
        private
        character(len=1), dimension(:),allocatable      :: name
        type (hydro_step_t), public, pointer                  :: hydro  
        type (time_t), pointer                                :: time   
        integer                                               :: input_file 
        integer                                               :: nxp        
        integer                                               :: nyp        
        integer                                               :: nzp        
        integer                                               :: nx         
        integer                                               :: ny         
        integer                                               :: nz         
        integer                                               :: dimension
        integer                                               :: rezone_type
        integer :: global_nx
        integer :: num_proc

        integer                                               :: n_materials     

        integer                                               :: wilkins_scheme 

        !       integer(hid_t)                                        :: main_hdf5_diagnostics_file_id

        real(8)                                               :: emf, emfm 

        type (artificial_viscosity_t)  , pointer :: a_visc              
        type (sound_velocity_t)        , pointer :: total_sound_vel     
        type (acceleration_t)          , pointer :: acceleration        
        type (vertex_mass_t)           , pointer :: total_vertex_mass   
        type (vertex_mass_t)           , pointer :: previous_vertex_mass
        type (temperature_t)           , pointer :: total_temperature   
        type (cell_mass_t)             , pointer :: total_cell_mass     
        type (cell_mass_t)             , pointer :: previous_cell_mass  
        type (pressure_t)              , pointer :: total_pressure      
        type (pressure_t)              , pointer :: total_pressure_sum  
        type (velocity_t)              , pointer :: velocity            
        type (density_t)               , pointer :: total_density       
        type (energy_t)                , pointer :: total_sie           
        type (volume_t)                , pointer :: total_volume        
        class (mesh_base_t)   , public , pointer :: mesh                
        type (vof_t)                   , pointer :: total_vof
        type (parallel_parameters_t)   , pointer :: parallel_params     
        type (materials_in_cells_t)    , pointer :: mat_cells           
        type (num_materials_in_cells_t), pointer :: num_mat_cells       
        type (data_t)                  , pointer :: total_inverse_vertex_mass  
        type (data_t)                  , pointer :: total_dp_de_deriv   
        type (data_t)                  , pointer :: total_dp_drho_deriv 
        type (data_t)                  , pointer :: total_dt_de_deriv   
        type (data_t)                  , pointer :: total_dt_drho_deriv 

        type (material_t)   , pointer             :: materials
        !        type (eos_wrapper_t), dimension(:), allocatable         :: total_eos
        type (textual_diagnostic_t), dimension(:), allocatable  :: textual_diagnostics
        !        type (textual_diagnostic_hdf5_t), dimension(:), allocatable  :: textual_diagnostics_hdf5
        !                type (silo_diagnostic_t), pointer                       :: silo_diagnostic
        !        type (plot_diagnostic_t), pointer                       :: plot_diagnostic

        type(boundary_parameters_t), pointer :: boundary_params

        type(communication_t), pointer            :: communication
        type(communication_parameters_t), pointer :: communication_parameters_cell
        type(communication_parameters_t), pointer :: communication_parameters_vertex
        type(communication_parameters_t), pointer :: communication_parameters_material
        integer                                   :: communicator
        !type (cr_t), pointer                                    :: cr

    contains
        procedure, public  :: Write_to_files
        procedure, public  :: Close_files
        procedure, public  :: Start_calculation
        procedure, public  :: Get_parallel_parameters
        procedure, public  :: Get_time
        procedure, private :: Create_materials
        procedure, private :: Create_inverse_vertex_mass
        procedure, private :: Initialize_sie
        procedure, private :: Initialize_communication
        procedure, private :: Initialize_openmp

        procedure, private :: Set_communication
    end type problem_t


    interface problem_t
        module procedure Constructor
    end interface problem_t

contains

    type(problem_t) function Constructor(df)
        use general_utils_module, only : Compare_string, Str_eqv
        use geometry_module            , only : Create_start_index_layer

        type(datafile_t), intent(in)                      :: df        
        logical                                           :: saw_index
        integer                                           :: i, j, k, lay_j, lay_i, i_curr, j_curr, num_mat, ele_i, ele_j, nc, f_id
        integer, dimension(:)      , allocatable          :: witness_index
        type(cell_bc_wrapper_t  )  , allocatable, target  :: bc_c_wrap     
        type(vertex_bc_wrapper_t)  , allocatable, target  :: bc_v_wrap     
        type(vertex_bc_wrapper_t)  , allocatable, target  :: bc_v_wrap_coordinates     
        type(cell_bc_wrapper_t  ), dimension(:), pointer  :: bc_c_wrap_arr 
        type(vertex_bc_wrapper_t), dimension(:), pointer  :: bc_v_wrap_arr 
        type(vertex_bc_wrapper_t), dimension(:), pointer  :: bc_v_wrap_coordinates_arr 
        type(slip_cell_t        )               , target  :: s_c_bc        
        type(slip_cell_3d_t     )               , target  :: s_c_3d_bc     
        type(slip_vertex_t      )               , target  :: s_v_bc        
        type(slip_vertex_3d_t      )            , target  :: s_v_3d_bc     
        type(no_slip_vertex_3d_t      )         , target  :: no_s_v_3d_bc     
        type(lagrange_surface_cell_t      )          , target  :: lag_c_bc     
        type(lagrange_surface_vertex_t      )               , target  :: lag_v_bc     
        type(boundary_mirror_2d_t      )        , target  :: reflect_v_bc     
        type(boundary_mirror_3d_t      )        , target  :: reflect_v_bc_3d     
        type (mesh_t)                           , pointer :: mesh_2d                
        type (mesh_3d_t)                        , pointer :: mesh_3d                
        real(8)                                           :: emfm, emf     
        real(8), dimension(:, :, :), pointer              :: x, y, z, vol
        integer, dimension(:)      , allocatable          :: boundary_type
        real(8), dimension(:, :, :), pointer              :: ad1, ad2, ad3
        integer                                           :: counter, text_diag_counter, hdf5_diag_counter, tmp_mat
        character(len=80)       :: word
        integer, dimension(:, :), allocatable :: start_index
        real(8), dimension(:,:,:,:), pointer :: cell_mass_vof,density_vof

        integer :: myid, numprocs, ierr
        integer :: nxp, nyp, nzp, nx,ny,nz, m
        character(5) :: my_id
        integer :: rank 

        call Constructor%Initialize_communication(df)
        call Constructor%Initialize_openmp(df%threads)

        !Constructor%name = df%name
        !Constructor%name = "a"
        Constructor%emf  = df%emf
        Constructor%emfm = df%emfm

        nxp = Constructor%nxp
        nyp = Constructor%nyp
        nzp = Constructor%nzp

        nx = Constructor%nx
        ny = Constructor%ny
        nz = Constructor%nz


        Constructor%dimension  = df%dimension
        Constructor%n_materials = df%reduct_num_mat
        Constructor%wilkins_scheme = df%sw_wilkins
        Constructor%rezone_type = df%rezone_type
        nc = 1




        allocate(Constructor%time)
        allocate(Constructor%total_pressure_sum)
        allocate(Constructor%a_visc)

        allocate(Constructor%total_sound_vel)
        allocate(Constructor%acceleration)
        allocate(Constructor%total_vertex_mass)
        allocate(Constructor%total_temperature)
        allocate(Constructor%total_cell_mass)
        allocate(Constructor%total_pressure)
        allocate(Constructor%velocity)
        allocate(Constructor%total_density)
        allocate(Constructor%total_sie)
        allocate(Constructor%total_volume)
        if (df%dimension == 2) allocate(mesh_t:: Constructor%mesh)
        if (df%dimension == 3) allocate(mesh_3d_t:: Constructor%mesh)
        allocate(Constructor%total_vof)
        allocate(Constructor%mat_cells)
        allocate(Constructor%total_inverse_vertex_mass)
        allocate(Constructor%total_dp_de_deriv)
        allocate(Constructor%total_dp_drho_deriv)
        allocate(Constructor%total_dt_drho_deriv)
        allocate(Constructor%total_dt_de_deriv)
        allocate(Constructor%num_mat_cells)
        allocate(Constructor%hydro)
        !        allocate(Constructor%cr)
        allocate(textual_diagnostic_t :: Constructor%textual_diagnostics(0:df%num_diag_text))
         !       allocate(Constructor%textual_diagnostics_hdf5(0:df%num_diag_hdf5))

        allocate(material_t    :: Constructor%materials)

        !        allocate(eos_wrapper_t :: Constructor%total_eos (df%reduct_num_mat))
        allocate(Constructor%boundary_params)
        if (df%dimension == 2) allocate(bc_v_wrap_arr(4))
        if (df%dimension == 3) allocate(bc_v_wrap_arr(6))
        allocate(bc_v_wrap_coordinates_arr(1))
        allocate(bc_v_wrap_coordinates)
        allocate(bc_v_wrap)
        if (df%dimension == 2) allocate(bc_c_wrap_arr(4))
        if (df%dimension == 3) allocate(bc_c_wrap_arr(6))
        allocate(bc_c_wrap)



        if (Constructor%parallel_params%my_rank <= 9) then
            write(my_id,'(i1)') Constructor%parallel_params%my_rank
        elseif (Constructor%parallel_params%my_rank <= 99) then
            write(my_id,'(i2)') Constructor%parallel_params%my_rank
        else
            write(my_id,'(i3)') Constructor%parallel_params%my_rank
        end if

        write(*,*) trim(my_id)

        open(unit=69, file="diags" // trim(my_id), status="replace", action="write")

        if (Constructor%wilkins_scheme == 1) then 
            allocate(Constructor%previous_vertex_mass)
            allocate(Constructor%previous_cell_mass)
        end if

        Constructor%boundary_params = boundary_parameters_t(nxp, nyp, nzp, Constructor%dimension,Constructor%parallel_params&
            , df%mesh_type, df%boundary_conditions)



        if (df%dimension == 2) then
            if (df%mesh_type == 2) then
                bc_v_wrap%bc => s_v_bc
                bc_v_wrap_coordinates%bc => reflect_v_bc

                bc_v_wrap_arr(1) = bc_v_wrap
                bc_v_wrap_arr(2) = bc_v_wrap
                bc_v_wrap_arr(3) = bc_v_wrap
                bc_v_wrap_arr(4) = bc_v_wrap

                bc_c_wrap%bc => s_c_bc
                bc_c_wrap_arr(1) = bc_c_wrap
                bc_c_wrap_arr(2) = bc_c_wrap
                bc_c_wrap_arr(3) = bc_c_wrap
                bc_c_wrap_arr(4) = bc_c_wrap
            else if (df%mesh_type == 1) then

            end if

        else if (df%dimension == 3) then
            bc_v_wrap_coordinates%bc => reflect_v_bc_3d

            do i = 1, 6
                select case (df%boundary_conditions(i))
                    case (0) 
                        bc_c_wrap%bc => lag_c_bc
                        bc_v_wrap%bc => lag_v_bc
                        bc_c_wrap_arr(i) = bc_c_wrap
                        bc_v_wrap_arr(i) = bc_v_wrap
                    case (2) 
                        bc_v_wrap%bc => s_v_3d_bc
                        bc_c_wrap%bc => s_c_3d_bc
                        bc_c_wrap_arr(i) = bc_c_wrap
                        bc_v_wrap_arr(i) = bc_v_wrap
                    case (3) 
                        bc_v_wrap%bc => no_s_v_3d_bc
                        bc_c_wrap%bc => lag_c_bc
                        bc_v_wrap_arr(i) = bc_v_wrap
                        bc_c_wrap_arr(i) = bc_c_wrap
                    case default
                        write(*,*) "NO BC FOUND"
                        stop
                end select
            end do
        end if
        bc_v_wrap_coordinates_arr(1) = bc_v_wrap_coordinates

        if (df%dimension == 2) then
            Constructor%mat_cells = materials_in_cells_t(nxp, nyp, bc_c_wrap_arr,Constructor%boundary_params, 1&
                , df%number_layers_i, df%number_layers_j, df%number_cells_i&
                , df%number_cells_j, df%mat_index)

            allocate(mesh_2d)
            mesh_2d = mesh_t(df, bc_v_wrap_coordinates_arr, Constructor%boundary_params, Constructor%parallel_params)
            Constructor%mesh => mesh_2d

        else if (df%dimension == 3) then
            if (df%mesh_type == 1) then
                Constructor%mat_cells = materials_in_cells_t(nxp, nyp, nzp, bc_c_wrap_arr,Constructor%boundary_params,&
                    df%number_layers_i, df%number_cells_i, &
                    df%mat_index, Constructor%parallel_params)
            else
                if (df%mat_sie_0(1) == 0d0) then
                    Constructor%mat_cells = materials_in_cells_t(nxp, nyp, nzp, bc_c_wrap_arr,Constructor%boundary_params,&
                        1,df%mat_index, Constructor%parallel_params)
                else
                    Constructor%mat_cells = materials_in_cells_t(nxp, nyp, nzp, bc_c_wrap_arr,Constructor%boundary_params, Constructor%parallel_params)
                end if
            end if
            allocate(mesh_3d)
            mesh_3d = mesh_3d_t(df, bc_v_wrap_coordinates_arr, Constructor%boundary_params, Constructor%parallel_params)
            Constructor%mesh => mesh_3d
        end if
        write(*,*) "done mesh"

        Constructor%total_temperature   = temperature_t(df%init_temperature, nxp, nyp, nzp, bc_c_wrap_arr,&
            Constructor%boundary_params)
        Constructor%total_dp_de_deriv   = data_t  (nxp, nyp, nzp)
        Constructor%total_dp_drho_deriv = data_t  (nxp, nyp, nzp)
        Constructor%total_dt_de_deriv   = data_t  (nxp, nyp, nzp)
        Constructor%total_dt_drho_deriv = data_t  (nxp, nyp, nzp)
        Constructor%total_volume        = volume_t              (0d0, nxp, nyp, nzp, bc_c_wrap_arr,Constructor%boundary_params)
        Constructor%total_pressure      = pressure_t            (0d0, nxp, nyp, nzp, bc_c_wrap_arr,Constructor%boundary_params)
        Constructor%total_pressure_sum  = pressure_t            (0d0, nxp, nyp, nzp, bc_c_wrap_arr,Constructor%boundary_params)
        Constructor%total_sound_vel     = sound_velocity_t      (0d0, nxp, nyp, nzp, bc_c_wrap_arr,Constructor%boundary_params)

        Constructor%a_visc              = artificial_viscosity_t(0d0, df%quad_visc_fac, df%linear_visc_fac,&
            nxp, nyp, nzp, bc_c_wrap_arr,Constructor%boundary_params)
        Constructor%total_sie           = energy_t(0d0, nxp, nyp, nzp, bc_c_wrap_arr,Constructor%boundary_params)
        Constructor%total_vof           = vof_t    (nxp, nyp, nzp, bc_c_wrap_arr,Constructor%boundary_params, Constructor%mat_cells)
        Constructor%num_mat_cells       = num_materials_in_cells_t(1d0, nxp, nyp, nzp, bc_c_wrap_arr,Constructor%boundary_params)
        Constructor%total_inverse_vertex_mass = data_t  (0d0, nxp, nyp, nzp)

        Constructor%time = time_t(df, Constructor%communication, Constructor%parallel_params)

        call Constructor%total_volume%Point_to_data (vol)



        call Constructor%Create_materials (df, bc_c_wrap_arr, Constructor%mat_cells)


        Constructor%total_density = density_t    (df%mat_rho_0, Constructor%mat_cells , Constructor%nxp, Constructor%nyp&
            , Constructor%nzp , bc_c_wrap_arr,Constructor%boundary_params)

        Constructor%total_cell_mass = cell_mass_t (0d0, Constructor%nxp, Constructor%nyp, Constructor%nzp , bc_c_wrap_arr&
            ,Constructor%boundary_params)

        if (df%dimension == 3) then
            Constructor%velocity          = velocity_t    (0d0, nxp+1, nyp+1, nzp+1,&
                df%dimension, bc_v_wrap_arr,Constructor%boundary_params)
            Constructor%acceleration      = acceleration_t(0d0, nxp+1, nyp+1, nzp+1, bc_v_wrap_arr,Constructor%boundary_params)
            Constructor%total_vertex_mass = vertex_mass_t(0d0, nxp + 1, nyp + 1, &
                nzp + 1, bc_v_wrap_arr,Constructor%boundary_params )
        else
            Constructor%velocity          = velocity_t    (0d0, nxp+1, nyp+1, 1, df%dimension, bc_v_wrap_arr,&
                Constructor%boundary_params)
            Constructor%acceleration      = acceleration_t(0d0, nxp+1, nyp+1, 1, bc_v_wrap_arr,Constructor%boundary_params)
            Constructor%total_vertex_mass = vertex_mass_t(0d0, nxp + 1, nyp + 1, &
                1, bc_v_wrap_arr,Constructor%boundary_params )
        end if


        if (Constructor%wilkins_scheme == 1) then
            Constructor%previous_cell_mass = cell_mass_t (0d0, Constructor%nxp, Constructor%nyp, Constructor%nzp , bc_c_wrap_arr&
                ,Constructor%boundary_params)

            Constructor%previous_vertex_mass = vertex_mass_t(0d0,&
                Constructor%nxp + 1, Constructor%nyp + 1,Constructor%nzp + 1, bc_v_wrap_arr&
                ,Constructor%boundary_params)
        end if

        !df, nx, ny, nz, nxp, nyp, nzp, wilkins_scheme, mesh, velocity, acceleration,&
        !        total_volume, total_vof, total_sie, total_pressure, total_pressure_sum, total_density, &
        !        total_temperature, total_cell_mass, previous_cell_mass, vertex_mass, &
        !        previous_vertex_mass, inversed_vertex_mass, total_sound_vel,&
        !        a_visc, total_dp_de, total_dp_drho, total_dt_de, total_dt_drho, init_temperature, &
        !        nmats, materials, num_mat_cells, mat_id, emf, emfm, parallel_params, mat_ids

        Constructor%hydro = hydro_step_t(df, Constructor%nx , Constructor%ny , Constructor%nz, Constructor%nxp         ,&
            Constructor%nyp, Constructor%nzp, Constructor%wilkins_scheme, Constructor%mesh,&
            Constructor%velocity , Constructor%acceleration, Constructor%total_volume     ,&
            Constructor%total_vof, Constructor%total_sie   , Constructor%total_pressure   ,&
            Constructor%total_pressure_sum, Constructor%total_density                     ,&
            Constructor%total_temperature,  Constructor%total_cell_mass                   ,&
            Constructor%previous_cell_mass, Constructor%total_vertex_mass                 ,&
            Constructor%previous_vertex_mass, Constructor%total_inverse_vertex_mass       ,&
            Constructor%total_sound_vel,    Constructor%a_visc                            ,&
            Constructor%total_dp_de_deriv,  Constructor%total_dp_drho_deriv               ,&
            Constructor%total_dt_de_deriv,  Constructor%total_dt_drho_deriv               ,&
            df%init_temperature           ,&
            Constructor%n_materials, Constructor%materials, Constructor%num_mat_cells     ,&
            Constructor%mat_cells, Constructor%emf, Constructor%emfm, Constructor%parallel_params, df%mat_index)

        !Constructor%cr = cr_t(Constructor%hydro, Constructor%time,Constructor%boundary_params ,df%run_name, df%with_cr)

        !        if ( df%num_diag_hdf5 > 0 ) then
         !   Constructor%main_hdf5_diagnostics_file_id = Static_hdf5_init()
        !        end if

        text_diag_counter = 1
        hdf5_diag_counter = 1

        do i=1, size(df%diag_types(:))
            word = df%diag_types(i)
            if (Str_eqv(word, 'plot') .eqv. .TRUE.) then
            !                allocate(Constructor%plot_diagnostic)
            !                Constructor%plot_diagnostic = plot_diagnostic_t()
            !                call Constructor%plot_diagnostic%Init_diagnostic(Constructor%hydro, Constructor%time, 110 + i)
            else if (Str_eqv(word, 'silo') .eqv. .TRUE.) then
            !                allocate(Constructor%silo_diagnostic)
            !                Constructor%silo_diagnostic = silo_diagnostic_t()
            !                call Constructor%silo_diagnostic%Init_diagnostic(Constructor%hydro,Constructor%time)
            else if (Str_eqv(word, 'text') .eqv. .TRUE.) then
                word = df%diag_names(i)
                Constructor%textual_diagnostics(text_diag_counter) = textual_diagnostic_t(word, 110 + i, Constructor%parallel_params%my_rank)
                call Constructor%textual_diagnostics(text_diag_counter)%Init_diagnostic(Constructor%hydro,Constructor%time,counter&
                    , Constructor%parallel_params%my_rank)
                text_diag_counter = text_diag_counter + 1
            !                write(*,*), "problem:", word
            !else if (Str_eqv(word, 'hdf5') .eqv. .TRUE.) then
            !    word = df%diag_names(i)
            !    Constructor%textual_diagnostics_hdf5(hdf5_diag_counter) = textual_diagnostic_hdf5_t(word, 110 + i)
            !    call Constructor%textual_diagnostics_hdf5(hdf5_diag_counter)%Init_diagnostic(Constructor%hydro,Constructor%time&
            !    ,counter&
            !        , Constructor%main_hdf5_diagnostics_file_id)
            !    hdf5_diag_counter = hdf5_diag_counter + 1
            end if
        end do

        call Constructor%Set_communication()

        call Constructor%mesh%Exchange_virtual_space_blocking()
        call Constructor%total_vof%Exchange_virtual_space_blocking()

        write(*,*) "Done diagnostics"

        if(Constructor%dimension == 3) then

            call Constructor%boundary_params%Calculate_edge_vector_3d(Constructor%mesh%coordinates%data)
            call Constructor%boundary_params%Calculate_boundary_normal_3d(Constructor%mesh%coordinates%data)


            call Constructor%mesh%coordinates%Apply_boundary(Constructor%mesh%coordinates%data)
            call Constructor%mesh%coordinates%Exchange_virtual_space_blocking()

            call Constructor%mesh%Calculate_average_coordinates()

            do i = 1, 3
                call mesh_3d%avg_coordinates(i)%Exchange_virtual_space_blocking()
            end do

            call Constructor%mesh%Point_to_data(x,y,z)
            call Constructor%total_volume%Calculate(Constructor%mesh%coordinates)
            call Constructor%total_volume%Exchange_virtual_space_blocking()
            call Constructor%total_density%Exchange_virtual_space_blocking()


            call Constructor%total_cell_mass%Calculate_cell_mass(Constructor%total_volume, Constructor%total_density)
            call Constructor%total_cell_mass%Exchange_virtual_space_blocking()

            call Constructor%total_vertex_mass%Calculate_vertex_mass_3d(Constructor%mesh%coordinates, Constructor%total_density&
                , Constructor%total_cell_mass)

        else
            call Constructor%mesh%Point_to_data(x,y)
            call Constructor%total_volume%Calculate(Constructor%mesh%coordinates, cyl_optional=Constructor%mesh%cyl)
            call Constructor%total_volume%Exchange_virtual_space_blocking()
            call Constructor%total_cell_mass%Calculate_cell_mass(Constructor%total_volume, Constructor%total_density)
            call Constructor%total_cell_mass%Exchange_virtual_space_blocking()
            call Constructor%total_vertex_mass%Calculate_vertex_mass_2d(Constructor%mesh%coordinates, Constructor%total_density&
                , Constructor%total_cell_mass, Constructor%wilkins_scheme, Constructor%mesh%cyl)

        end if

        if (Constructor%wilkins_scheme == 1) then
            call Constructor%previous_cell_mass%Calculate_cell_mass(Constructor%total_volume, Constructor%total_density)

            if(Constructor%dimension == 3) then
                call Constructor%previous_vertex_mass%Calculate_vertex_mass_3d(Constructor%mesh%coordinates,&
                    Constructor%total_density&
                    , Constructor%total_cell_mass)
            else
                call Constructor%previous_vertex_mass%Calculate_vertex_mass_2d(Constructor%mesh%coordinates,&
                    Constructor%total_density&
                    , Constructor%total_cell_mass, Constructor%wilkins_scheme, Constructor%mesh%cyl)
            end if

        end if

        call Constructor%velocity%Apply_boundary (Constructor%mesh%coordinates%data)

        call Constructor%total_vertex_mass%Exchange_virtual_space_blocking()

        call Constructor%Create_inverse_vertex_mass(Constructor%emfm, Constructor%total_vertex_mass, &
            Constructor%total_inverse_vertex_mass)
        call Constructor%total_inverse_vertex_mass%Exchange_virtual_space_blocking()

        call Constructor%materials%density%point_to_data(density_vof)
        call Constructor%materials%cell_mass%point_to_data(cell_mass_vof)

        do k = 1, Constructor%nz
            do j = 1, Constructor%ny
                do i = 1, Constructor%nx
                    do tmp_mat=1, Constructor%n_materials
                        cell_mass_vof(tmp_mat, i,j,k) = vol(i,j,k) * density_vof(tmp_mat,i,j,k)
                    end do
                end do
            end do
        end do

        call Constructor%materials%sie%Exchange_virtual_space_blocking()
        call Constructor%materials%density%Exchange_virtual_space_blocking()
        call Constructor%materials%cell_mass%Exchange_virtual_space_blocking()
        call Constructor%materials%temperature%Exchange_virtual_space_blocking()


        call Constructor%Initialize_sie(Constructor%emf)
        call Constructor%a_visc%Initialize(df%quad_visc_fac, df%linear_visc_fac, df%to_radial_index_sphere,&
            df%from_radial_index_sphere, df%start_layer_index_r,df%no_xl_flag, df%start_no_xl_visc, df%end_no_xl_visc, df%offset_no_xl_visc)
        call Constructor%velocity%Initialize(df%to_radial_index_sphere,&
            df%from_radial_index_sphere, df%no_move_layer, df%start_layer_index_r)

        write (*,*) "finished building problem"

    end function


    subroutine Set_communication(this)
        class (problem_t), intent(in out) :: this
        type(communication_parameters_t), pointer       :: comm_params
        type(communication_t)    , pointer       :: communication
        integer, dimension(:), allocatable :: pdims
        logical, dimension(:), allocatable :: periods
        integer :: ierr
        integer :: i, num_dim, communicator

        call this%mat_cells%Set_communication(this%communication, this%communication_parameters_cell)

        call this%mesh%Set_communication(this%communication, this%communication_parameters_vertex)

        call this%total_temperature%Set_communication(this%communication, this%communication_parameters_cell)
        call this%velocity%Set_communication(this%communication, this%communication_parameters_vertex)
        call this%acceleration%Set_communication(this%communication, this%communication_parameters_vertex)
        call this%total_dp_de_deriv%Set_communication(this%communication, this%communication_parameters_cell)
        call this%total_dp_drho_deriv%Set_communication(this%communication, this%communication_parameters_cell)
        call this%total_dt_de_deriv%Set_communication(this%communication, this%communication_parameters_cell)
        call this%total_dt_drho_deriv%Set_communication(this%communication, this%communication_parameters_cell)
        call this%total_volume%Set_communication(this%communication, this%communication_parameters_cell)
        call this%total_pressure%Set_communication(this%communication, this%communication_parameters_cell)
        call this%total_pressure_sum%Set_communication(this%communication, this%communication_parameters_cell)
        call this%total_sound_vel%Set_communication(this%communication, this%communication_parameters_cell)
        call this%a_visc%Set_communication(this%communication, this%communication_parameters_cell)

        call this%total_sie%Set_communication(this%communication, this%communication_parameters_cell)
        call this%total_vof%Set_communication(this%communication, this%communication_parameters_cell)
        call this%num_mat_cells%Set_communication(this%communication, this%communication_parameters_cell)
        call this%total_inverse_vertex_mass%Set_communication(this%communication, this%communication_parameters_cell)

        !        do i = 1, this%n_materials
        call this%materials%Set_communication_material(this%communication, this%communication_parameters_material)
        !        end do

        call this%total_density%Set_communication(this%communication, this%communication_parameters_cell)
        call this%total_cell_mass%Set_communication(this%communication, this%communication_parameters_cell)
        call this%total_vertex_mass%Set_communication(this%communication, this%communication_parameters_vertex)

        call this%boundary_params%Set_communication(this%communication, this%communication_parameters_cell)

        if (this%wilkins_scheme == 1) then
            call this%previous_cell_mass%Set_communication(this%communication, this%communication_parameters_cell)
            call this%previous_vertex_mass%Set_communication(this%communication, this%communication_parameters_vertex)
        end if

        call this%hydro%Set_communication(this%communication, this%communication_parameters_cell, &
            this%communication_parameters_vertex, this%communication_parameters_material)
        return
    end subroutine Set_communication

    subroutine Write_to_files(this)
        class (problem_t), intent(in out) :: this
        integer :: i
        integer, save :: counter_diag = 0

        do i=1,size(this%textual_diagnostics(1:))
            call this%textual_diagnostics(i)%Apply()
        !!            write(*,*) "WRITING"
        end do
        !if ( associated(this%silo_diagnostic) ) call this%silo_diagnostic%Apply()
        !       do i=1,size(this%textual_diagnostics_hdf5(1:))
        !           call this%textual_diagnostics_hdf5(i)%Apply
        !       end do
        !if (associated(this%plot_diagnostic)) call this%plot_diagnostic%Apply
        !if (mod(this%hydro%cyc_delete, 100) == 0 .or. counter_diag == 0) then


        !end if
        counter_diag = counter_diag + 1
    end subroutine Write_to_files

    subroutine Close_files(this)
        class (problem_t), intent(in out) :: this
        integer :: i, error_hdf5

        close(69)
        do i=1,size(this%textual_diagnostics(1:))
            call this%textual_diagnostics(i)%Close_diagnostic
        end do

        !        do i=1,size(this%textual_diagnostics_hdf5(1:))
        !            call this%textual_diagnostics_hdf5(i)%Close_diagnostic
        !        end do

    !        if ( associated(this%plot_diagnostic) ) call this%plot_diagnostic%Close_diagnostic
    !        if ( associated(this%silo_diagnostic) ) call this%silo_diagnostic%Close_diagnostic

        !if ( size(this%textual_diagnostics_hdf5(1:)) > 0 ) error_hdf5 = Static_hdf5_close(this%main_hdf5_diagnostics_file_id)

    end subroutine Close_files


    subroutine Start_calculation(this)
        use omp_lib
        use general_utils_module, only : int2str
        class (problem_t) , intent(in out) :: this 
        real(8) :: reem_start, reem_total
        integer :: ierr
        !character(len=size(this%cr%run_name)+size(this%name)+2)  :: ckpt_name
        integer :: np
        integer :: first_comm, comm, rank
        integer :: left, right, rear, front, down, up
        integer :: npx, npy, npz 
        integer :: i, j, k, n_l, xs, ys, zs
      
        real(8), dimension(:,:,:), pointer  :: values
        real(8), dimension(:), allocatable  :: array
        real(8) :: val
        integer, save :: counter = 0
        integer :: comm_dist_graph
        integer                              :: my_rank, num_neighbors
        integer                              :: x, y, z   
        integer, dimension(26)   :: destinations, weights
        integer, dimension(1)   :: sources, degrees
        integer, dimension(3)   :: coords
        type(communication_parameters_t), pointer       :: comm_params
        type(communication_t)    , pointer       :: communication
        integer :: aaaai
        integer :: ncyc, max_ncyc
        !        call this%cr%Get_ckpt_name(this%name, ckpt_name)
#ifdef DEBUG
        write(*,*) '@@@@ ATTACHE TO PROCESS AND CHANGE i to 0 TO CONTINUE EXECUTION @@@@'
        aaaai = 1
        do while ( aaaai .eq. 1)
            aaaai = aaaai
        end do
        write(*,*) '@@@@@@@@@@@@@@@@@@@@@@@@ STARTED EXECUTION @@@@@@@@@@@@@@@@@@@@@@@@'
#endif

        !call this%cr%Restart(ckpt_name)

        reem_total = omp_get_wtime()
                        !call this%Write_to_files()
        ncyc = 1
        if (this%rezone_type == 0) then
            max_ncyc = 21
        else
            max_ncyc = 201
        end if

        if (this%mesh%dimension == 2) then
            do while (this%time%Should_continue() .and. ncyc < max_ncyc)
                call this%hydro%do_time_step_2d(this%time)
                call this%time%Update_time()
                !call this%Write_to_files()
                ncyc = ncyc + 1
            !       call this%cr%Checkpoint(ckpt_name)
            end do

        else if (this%mesh%dimension == 3) then
            do while (this%time%Should_continue() .and. ncyc < max_ncyc)
                reem_start = omp_get_wtime()
                call this%hydro%do_time_step_3d(this%time)
                call this%time%Update_time()
                !  call this%Write_to_files()
                counter = counter + 1
                ncyc = ncyc + 1
            !      call this%cr%Checkpoint(ckpt_name)
            end do
        end if



        call this%Close_files()
        if (this%parallel_params%my_rank == 0) then

!            write(str1, *) this%global_nx
!            write(str2, *) this%num_proc
!            write(*,*) str1,str2
!            open (71, status = 'replace', file = trim(tmp2(:len(trim(tmp2)) -4 ) // "_" &
!                // trim(tmp) // ".txt"))
            open (71, status = 'replace', file = "runtime_" // trim(int2str(this%global_nx)) // "_" // trim(int2str(this%num_proc)) // ".txt")
            write(71,*) "Total Time:", omp_get_wtime() - reem_total
        end if
        write(*,*) "Total Time:", omp_get_wtime() - reem_total
        write(*,*) "ncyc: ", ncyc-1
    end subroutine Start_calculation

    subroutine Create_materials(this, df, bc_c_wrap_arr, mat_cell)
        use datafile_module
        implicit none
        class(problem_t)                                , intent(inout) :: this          
        type(datafile_t)                                , intent(in)    :: df            
        type(cell_bc_wrapper_t  ), dimension(:), pointer, intent(inout) :: bc_c_wrap_arr 
        type(materials_in_cells_t), pointer                      , intent(inout) :: mat_cell
        !        type(eos_wrapper_t), allocatable                                :: eos_c_wrap
        !        type(ideal_gas_t), target                                       :: ig_eos_c
        integer :: nxp, nyp, nzp
        integer                                                         :: i, j, k, p, mat
        integer, dimension(:), allocatable                              :: witness_index
        logical                                                         :: saw_index
        type(material_t), pointer :: mat1
        nxp = this%parallel_params%nxp
        nyp = this%parallel_params%nyp
        nzp = this%parallel_params%nzp

        !        allocate(eos_c_wrap)
        !        eos_c_wrap%eos => ig_eos_c
        saw_index = .true.
        i = 1
        allocate(witness_index(0:df%n_materials))
        witness_index = -1
        !        do i = 1, this%n_materials
        !            this%total_eos(i) = eos_c_wrap
        !        end do

        !                    write(*,*) "MATS", i, nxp, nyp, nzp, mat, df%mat_gamma_gas(mat),&
        !                        df%mat_atomic_mass(mat), df%mat_z(mat), df%mat_z2(mat),&
        !                        df%mat_rho_0(mat), df%init_temperature, df%mat_sie_0(mat)
        this%materials = material_t(nxp, nyp, nzp, this%n_materials, df%mat_index, df%mat_gamma_gas,&
            df%mat_atomic_mass, df%mat_z, df%mat_z2,&
            df%mat_rho_0, df%init_temperature, df%mat_sie_0, mat_cell, bc_c_wrap_arr, this%boundary_params)

        !        nxp, nyp, nzp, nmats, mat_ids, gamma_gas, atomic_mass,&
        !        num_protons, num_protons_2, rho_0, temperature_init, sie_0, eos, mat_cells, bc_cell&
        !        , bc_params
        deallocate(witness_index)
    end subroutine Create_materials

    subroutine Create_inverse_vertex_mass(this, emfm, tar_total_vertex_mass, tar_total_inverse_vertex_mass)
        class(problem_t)           , intent(inout) :: this                          
        real(8)                    , intent(in)    :: emfm                          
        type(data_t)       , target, intent(inout) :: tar_total_inverse_vertex_mass 
        type(vertex_mass_t), target, intent(inout) :: tar_total_vertex_mass         

        integer :: i, j, k
        real(8), dimension(:, :, :), pointer    :: v_mass, i_v_mass

        call tar_total_vertex_mass%Point_to_data(v_mass)
        call tar_total_inverse_vertex_mass%Point_to_data(i_v_mass)
        do k = 1, this%nzp
            do j = 1, this%nyp
                do i = 1, this%nxp
                    if (v_mass(i, j, k) > emfm) then
                        i_v_mass(i, j, k) = 1.d0 / v_mass(i, j, k)
                    else
                        i_v_mass(i, j, k) = 0.d0
                    end if
                end do
            end do
        end do

    end subroutine Create_inverse_vertex_mass


    subroutine Initialize_sie(this, emf)
        implicit none
        class(problem_t), intent(inout)      :: this             
        real(8)         , intent(in)         :: emf              


        integer                              :: i, j, k, tmp_mat
        real(8), dimension(:, :,:, :), pointer :: sie_vof
        real(8), dimension(:, :,:, :), pointer :: cell_mass_vof
        real(8), dimension(:, :, :), pointer :: cell_mass    
        real(8), dimension(:, :, :), pointer :: sie
        real(8), dimension(:, :, :), pointer :: t

        call this%materials%sie      %Point_to_data(sie_vof)
        call this%materials%cell_mass%Point_to_data(cell_mass_vof)

        call this%total_sie%Point_to_data(sie)
        call this%total_cell_mass%Point_to_data(cell_mass)
        call this%materials%Apply_eos(this%nx, this%ny, this%nz,emf,.false.)



        do k = 1, this%nz
            do j = 1, this%ny
                do i = 1, this%nx
                    do tmp_mat = 1, this%n_materials
                        sie(i, j, k) = sie(i, j, k) + sie_vof(tmp_mat, i, j, k) * cell_mass_vof(tmp_mat, i, j, k)
                    end do
                end do
            end do
        end do

        do k = 1, this%nz
            do j = 1, this%ny
                do i = 1, this%nx
                    sie(i, j, k) = sie(i, j, k) / (cell_mass(i, j, k) + 1.d-30)

                end do
            end do
        end do

    end subroutine Initialize_sie


    subroutine Initialize_communication(this, df)
        use communication_utils_module      , only : Count_neighbors, Create_diag_topology
        implicit none
        class(problem_t), intent(inout)      :: this
        class(datafile_t), intent(in) :: df

        integer, dimension(3)   :: coords
        integer                              :: my_rank, my_rank1, num_neighbors
        integer                              :: local_nxp, local_nyp, local_nzp, local_nx, local_ny, local_nz
        integer, dimension(:), allocatable :: pdims
        logical, dimension(:), allocatable :: periods
        integer :: ierr
        integer :: i, communicator
        local_nx = df%virt_nx / df%npx
        local_ny = df%virt_ny / df%npy

        local_nz = 1
        local_nzp = 1
        if (df%dimension == 3) then
            local_nz = df%virt_nz / df%npz
            local_nzp = local_nz + 1
        end if
        local_nxp = local_nx + 1
        local_nyp = local_ny + 1
        this%global_nx = df%virt_nz
        this%num_proc = df%npx * df%npy * df%npz
        this%nxp = local_nxp
        this%nyp = local_nyp
        this%nzp = local_nzp
        this%nx = local_nx
        this%ny = local_ny
        this%nz = local_nz

        allocate(this%communication_parameters_cell)
        allocate(this%communication_parameters_vertex)
        allocate(this%communication_parameters_material)

        call MPI_comm_rank(MPI_COMM_WORLD, my_rank1, ierr)

        coords(1) = mod(my_rank1, df%npx) + 1
        coords(2) = (mod(my_rank1, (df%npx * df%npy))) / df%npx + 1
        coords(3) = my_rank1 / (df%npx * df%npy) + 1

        call Create_diag_topology(df%npx, df%npy, df%npz, coords, communicator, num_neighbors)

        call MPI_comm_rank(communicator, my_rank, ierr)

        coords(1) = mod(my_rank, df%npx) + 1
        coords(2) = (mod(my_rank, (df%npx * df%npy))) / df%npx + 1
        coords(3) = my_rank / (df%npx * df%npy) + 1
        this%communication_parameters_cell = communication_parameters_t(this%nx, this%ny, this%nz, num_neighbors, communicator,&
            df%npx, df%npy, df%npz, coords, 0)

        this%communication_parameters_vertex = communication_parameters_t(this%nxp, this%nyp, this%nzp, num_neighbors, communicator,&
            df%npx, df%npy, df%npz, coords, 1)


        this%communication_parameters_material = communication_parameters_t(df%reduct_num_mat, this%nxp, this%nyp, this%nzp, num_neighbors, communicator,&
            df%npx, df%npy, df%npz, coords, 0)


        allocate(this%parallel_params)
        this%parallel_params = parallel_parameters_t(my_rank, coords, df%np, df%npx, df%npy, df%npz, df%virt_nxp, df%virt_nyp&
            , df%virt_nzp&
            , local_nxp, local_nyp, local_nzp, local_nx&
            , local_ny, local_nz)
        allocate(this%communication)
        this%communication = communication_t(this%parallel_params)
    end subroutine Initialize_communication


    subroutine Initialize_openmp(this, max_threads)
        use omp_lib
        implicit none
        class(problem_t), intent(inout)   :: this
        integer         , intent(in)      :: max_threads

        integer :: threads

        call omp_set_num_threads(max_threads)

        threads = omp_get_num_threads()

        if (threads /= max_threads) then
            write(*,*) "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
            write(*,*) "Warning, number of threads doesn't match number of threads set in datafile"
            write(*,*) "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"
        end if
    end subroutine Initialize_openmp


    function Get_time (this)
        class (problem_t), intent(in out) :: this  
        type(time_t) :: Get_time

        Get_time = this%time

    end function Get_time

    subroutine Get_parallel_parameters(this, parallel_param)
        class (problem_t), intent(in out) :: this  
        type(parallel_parameters_t),pointer, intent(inout) :: parallel_param

        parallel_param => this%parallel_params

    end subroutine Get_parallel_parameters
end module problem_module
