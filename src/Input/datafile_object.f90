
module datafile_module
    use json_module
    use, intrinsic :: iso_fortran_env , only: error_unit, output_unit
    use general_utils_module, only: str_eqv, int2str, error_msg, warn_msg, lower
    implicit none
    private
    public :: datafile_t

    type :: datafile_t
        private


        integer, public :: virt_nxp
        integer, public :: virt_nyp
        integer, public :: virt_nzp
        integer, public :: virt_nx
        integer, public :: virt_ny
        integer, public :: virt_nz

        integer, public :: with_cr                

        character(len=255), public        :: name

        integer, public ::  sw_vert_mass      
        integer, public ::  sw_nraz           
        integer, public ::  sw_symmetry       
        integer, public ::  sw_wilkins        

        real(8), public :: finish_time        
        real(8), public :: cyl        
        real(8), public :: dt_max             
        real(8), public :: init_temperature   
        real(8), public :: dt_factor          
        real(8), public :: t0                 
        real(8), public :: dt0                
        real(8), public :: dt_cour_fac        
        real(8), public :: emf
        real(8), public :: emfm
        real(8), public :: quad_visc_fac
        real(8), public :: linear_visc_fac
        integer, public :: to_radial_index_sphere 
        integer, public :: from_radial_index_sphere 
        integer, public :: no_move_layer           

        integer, public :: start_no_xl_visc
        integer, public :: end_no_xl_visc
        logical, public :: no_xl_flag
        integer, public :: offset_no_xl_visc

        integer, public :: np  
        integer, public :: npx 
        integer, public :: npy 
        integer, public :: npz 
        integer, public :: threads 

        integer, public                              :: mesh_type           
        integer, dimension(:), allocatable, public   :: boundary_conditions

        integer, public                              :: number_layers_i  
        integer, public                              :: number_layers_j  
        integer, public                              :: number_layers_k  
        integer, dimension(:), allocatable, public   :: number_cells_i
        integer, dimension(:), allocatable, public   :: number_cells_j
        integer, dimension(:), allocatable, public   :: number_cells_k


        integer                           , public   :: reduct_num_mat   
        integer                           , public   :: n_materials   
        integer, dimension(:), allocatable, public   :: mat_index
        integer, dimension(:), allocatable, public   :: mat_ids
        integer, dimension(:), allocatable, public   :: mat_eos
        real(8), dimension(:), allocatable, public   :: mat_atomic_mass
        real(8), dimension(:), allocatable, public   :: mat_z
        real(8), dimension(:), allocatable, public   :: mat_z2
        real(8), dimension(:), allocatable, public   :: mat_rho_0
        real(8), dimension(:), allocatable, public   :: mat_sie_0
        real(8), dimension(:), allocatable, public   :: mat_gamma_gas

        integer, dimension(:, :), allocatable, public    :: start_layer_index_r
        integer, dimension(:, :), allocatable, public    :: start_layer_index_t
        integer, dimension(:, :), allocatable, public    :: start_layer_index_p

        integer, dimension(:)   , allocatable, public     :: contour_i_type
        real(8), dimension(:,:) , allocatable, public     :: contour_i_line
        real(8), dimension(:)   , allocatable, public     :: zone_i_d
        integer, dimension(:)   , allocatable, public     :: zone_i_type


        real(8), dimension(:, :), allocatable, public     :: theta0
        integer, dimension(:)   , allocatable, public     :: contour_j_type
        real(8), dimension(:)   , allocatable, public     :: zone_j_d
        integer, dimension(:)   , allocatable, public     :: zone_j_type

        real(8), dimension(:, :), allocatable, public     :: phi0
        integer, dimension(:)   , allocatable, public     :: contour_k_type
        real(8), dimension(:)   , allocatable, public     :: zone_k_d
        integer, dimension(:)   , allocatable, public     :: zone_k_type

        integer, public :: rezone_type
        logical, public :: shorter_advect
        logical, public :: fix_overflow
        logical, public :: line_calc
        integer, public :: dimension


        integer                                       , public      :: num_diag_text 
        integer                                       , public      :: num_diag_hdf5 
        character(len=80), dimension(:), allocatable, public        :: diag_types
        character(len=80), dimension(:), allocatable, public        :: diag_variables
        character(len=80), dimension(:), allocatable, public        :: diag_names

    contains
           procedure, public :: Parse_diagnostics
           procedure, public :: Parse_data
           procedure, public :: Parse_switches
           procedure, public :: Parse_simulation_parameters
           procedure, public :: Parse_rezone_advect
           procedure, public :: Parse_mesh
           procedure, public :: Parse_layers
           procedure, public :: Parse_materials
           procedure, public :: Parse_contours
           procedure, public :: Parse_zones
           procedure, public :: Parse_parallel
           procedure, public :: Parse_checkpoint_restart
    end type datafile_t

    interface datafile_t
        module procedure Constructor

    end interface datafile_t
contains
    type(datafile_t) function Constructor(datafile)
        use replace_words_module, only : Replace_words
        use defaults_module, only : Set_defaults
        use geometry_module, only : Create_start_index_layer

        implicit none

        character(len=*), intent(in)    :: datafile 
        type(json_file) :: json
        type(json_value), pointer :: p
        integer :: i,j,k
        logical :: found

        call json%initialize()

        call json%load(filename = datafile)

        call Set_defaults(json)

        call Replace_words(json)

        call Constructor%Parse_data(json)

        call Constructor%Parse_switches(json)
        call Constructor%Parse_simulation_parameters(json)
        call Constructor%Parse_rezone_advect(json)
        call Constructor%Parse_mesh(json)
        call Constructor%Parse_layers(json)
        call Constructor%Parse_materials(json)

        call Constructor%Parse_contours(json)
        call Constructor%Parse_zones(json)
        call Constructor%Parse_parallel(json)
        call Constructor%Parse_checkpoint_restart(json)
        call Constructor%Parse_diagnostics(json)

        if (json%failed()) then
            call error_msg("Json parser failed")

        end if


    end function Constructor

    subroutine Parse_diagnostics(this, json)
        implicit none
        class (datafile_t) , intent(inout) :: this
        type(json_file)    , intent(inout) :: json

        character(len=*), parameter    :: segment  = "diagnostics"
        character(len=255) :: tmp
        character(len=255) :: path_to_diagnostic
        character(len=255) :: root
        character(len=:),allocatable :: diag, diag_variable, diag_type
        logical :: found
        integer :: num_grps, i, j, num_diags,k,num_variables, total_files, diag_counter

        this%num_diag_text = 0
        tmp = segment // ".number_diagnostics"
        root = segment // ".group"
        call json%get(trim(tmp), num_grps, found)
        total_files = 0
        do i = 1, num_grps
            tmp = root
            path_to_diagnostic = trim(tmp) // "(" // int2str(i) // ").diagnostic"
            call json%info(trim(path_to_diagnostic), n_children=num_diags, found=found)
            if (found .eqv. .false. ) cycle

            do j = 1, num_diags
                tmp = trim(path_to_diagnostic) // "(" // int2str(j) // ").variables"
                call json%info(trim(tmp), n_children=num_variables, found=found)
                if (found .eqv. .false.) total_files = total_files + 1
                total_files = total_files + num_variables
            end do
        end do
        allocate(this%diag_names     (total_files))
        allocate(this%diag_types     (total_files))
        allocate(this%diag_variables (total_files))

        diag_counter = 0
        do i = 1, num_grps
            tmp = root
            path_to_diagnostic = trim(tmp) // "(" // int2str(i) // ").diagnostic"
            call json%info(trim(path_to_diagnostic), n_children=num_diags, found=found)
            if (found .eqv. .false. ) cycle
            do j = 1, num_diags
                tmp = trim(path_to_diagnostic) // "(" // int2str(j) // ").type"
                call json%get(trim(tmp), diag_type, found)

                root = trim(path_to_diagnostic) // "(" // int2str(j) // ")"
                select case(trim(diag_type))
                    case("text")
                        tmp = trim(root) // ".variables"
                        call json%info(trim(tmp), n_children=num_variables, found=found)
                        do k = 1, num_variables
                            diag_counter = diag_counter + 1
                            this%num_diag_text = this%num_diag_text + 1
                            tmp = trim(root) // ".variables(" // int2str(k) // ")"

                            call json%get(trim(tmp), diag, found)
                            if (found .eqv. .false.) call warn_msg("Warning no variable found in text or hdf5 segment")

                            this%diag_types     (diag_counter) = "text"
                            this%diag_names     (diag_counter) = lower(diag) // ".txt"
                            this%diag_variables (diag_counter) = lower(diag)
                        end do
                    case ("hdf5")
                        call json%info(trim(tmp), n_children=num_variables, found=found)
                        do k = 1, num_variables
                            diag_counter = diag_counter + 1
                            this%num_diag_hdf5 = this%num_diag_hdf5 + 1
                            tmp = trim(root) // ".variables(" // int2str(k) // ")"
                            call json%get(trim(tmp), diag, found)
                            if (found .eqv. .false.) call warn_msg("Warning no variable found in text or hdf5 segment")
                            this%diag_types     (diag_counter) = "hdf5"
                            this%diag_names     (diag_counter) = lower(diag)
                            this%diag_variables (diag_counter) = lower(diag)
                        end do
                    case ("silo")
                        diag_counter = diag_counter + 1
                        tmp = trim(root)
                        this%diag_types     (diag_counter) = "silo"
                        this%diag_names     (diag_counter) = ""
                        this%diag_variables (diag_counter) = ""
                    case ("plot")
                        diag_counter = diag_counter + 1
                        this%diag_types     (diag_counter) = "plot"
                        this%diag_names     (diag_counter) = "plot"
                        this%diag_variables (diag_counter) = ""
                end select
            end do
        end do
    end subroutine Parse_diagnostics

    subroutine Parse_parallel(this, json)
        implicit none
        class (datafile_t) , intent(inout) :: this
        type(json_file)    , intent(inout) :: json

        character(len=*), parameter    :: segment  = "parallel"
        character(len=255) :: tmp
        character(len=255) :: root
        character(len=:),allocatable :: parallel_param
        logical :: found
        integer :: var

        tmp = segment // ".np"
        call json%get(trim(tmp), var, found)
        if (found .eqv. .false.) call error_msg("Bad value in parallel segment")
        this%np = var

        tmp = segment // ".npx"
        call json%get(trim(tmp), var, found)
        if (found .eqv. .false.) call error_msg("Bad value in parallel segment")
        this%npx = var

        tmp = segment // ".npy"
        call json%get(trim(tmp), var, found)
        if (found .eqv. .false.) call error_msg("Bad value in parallel segment")
        this%npy = var

        tmp = segment // ".npz"
        call json%get(trim(tmp), var, found)
        if (found .eqv. .false.) call error_msg("Bad value in parallel segment")
        this%npz = var

        tmp = segment // ".threads"
        call json%get(trim(tmp), var, found)
        if (found .eqv. .false.) call error_msg("Bad value in parallel segment")
        this%threads = var
    end subroutine Parse_parallel

    subroutine Parse_data(this, json)
        implicit none
        class (datafile_t) , intent(inout) :: this
        type(json_file)    , intent(inout) :: json

        character(len=*), parameter    :: segment  = "data"
        character(len=255) :: tmp
        character(len=255) :: root
        character(len=:),allocatable :: parallel_param
        logical :: found
        integer :: var

        tmp = segment // ".name"


        call json%get(trim(tmp), parallel_param, found)

        if (found .eqv. .false.) call error_msg("Bad value in data segment")
        this%name = parallel_param

    end subroutine Parse_data

    subroutine Parse_switches(this, json)
        implicit none
        class (datafile_t) , intent(inout) :: this
        type(json_file)    , intent(inout) :: json

        character(len=*), parameter    :: segment  = "switches"
        character(len=255) :: tmp
        character(len=255) :: root
        logical :: found, var_l
        integer :: var

        tmp = segment // ".sw_nraz"
        call json%get(trim(tmp), var, found)
        if (found .eqv. .false.) call error_msg("Bad value in switch segment")
        this%sw_nraz = var

        tmp = segment // ".sw_symmetry"
        call json%get(trim(tmp), var, found)
        if (found .eqv. .false.) call error_msg("Bad value in switch segment")
        this%sw_symmetry = var

        tmp = segment // ".sw_vert_mass"
        call json%get(trim(tmp), var, found)
        if (found .eqv. .false.) call error_msg("Bad value in switch segment")
        this%sw_vert_mass = var

        tmp = segment // ".sw_wilkins"
        call json%get(trim(tmp), var, found)
        if (found .eqv. .false.) call error_msg("Bad value in switch segment")
        this%sw_wilkins = var

        tmp = segment // ".sw_cr"
        call json%get(trim(tmp), var, found)
        if (found .eqv. .false.) call error_msg("Bad value in switch segment")
        this%with_cr = var

    end subroutine Parse_switches

    subroutine Parse_simulation_parameters(this, json)
        implicit none
        class (datafile_t) , intent(inout) :: this
        type(json_file)    , intent(inout) :: json

        character(len=*), parameter    :: segment  = "simulation_parameters"
        character(len=255) :: tmp
        character(len=255) :: root
        logical :: found, var_l
        real(8) :: var, var_r
        integer :: var_i

        tmp = segment // ".time_final"
        call json%get(trim(tmp), var, found)
        if (found .eqv. .false.) call error_msg("Bad value in simulation_parameters segment")
        this%finish_time = var

        tmp = segment // ".dt_max"
        call json%get(trim(tmp), var, found)
        if (found .eqv. .false.) call error_msg("Bad value in simulation_parameters segment")
        this%dt_max = var

        tmp = segment // ".cyl"
        call json%get(trim(tmp), var, found)
        if (found .eqv. .false.) call error_msg("Bad value in simulation_parameters segment")
        this%cyl = var

        tmp = segment // ".init_temperature"
        call json%get(trim(tmp), var, found)
        if (found .eqv. .false.) call error_msg("Bad value in simulation_parameters segment")
        this%init_temperature = var

        tmp = segment // ".dt_factor"
        call json%get(trim(tmp), var, found)
        if (found .eqv. .false.) call error_msg("Bad value in simulation_parameters segment")
        this%dt_factor = var

        tmp = segment // ".dt0"
        call json%get(trim(tmp), var, found)
        if (found .eqv. .false.) call error_msg("Bad value in simulation_parameters segment")
        this%dt0 = var

        tmp = segment // ".time_final"
        call json%get(trim(tmp), var, found)
        if (found .eqv. .false.) call error_msg("Bad value in simulation_parameters segment")
        this%finish_time = var

        tmp = segment // ".dt_cour_fac"
        call json%get(trim(tmp), var, found)
        if (found .eqv. .false.) call error_msg("Bad value in simulation_parameters segment")
        this%dt_cour_fac = var

        tmp = segment // ".quad_visc_fac"
        call json%get(trim(tmp), var, found)
        if (found .eqv. .false.) call error_msg("Bad value in simulation_parameters segment")
        this%quad_visc_fac = var

        tmp = segment // ".linear_visc_fac"
        call json%get(trim(tmp), var, found)
        if (found .eqv. .false.) call error_msg("Bad value in simulation_parameters segment")
        this%linear_visc_fac = var

        tmp = segment // ".emf"
        call json%get(trim(tmp), var, found)
        if (found .eqv. .false.) call error_msg("Bad value in simulation_parameters segment")
        this%emf = var

        tmp = segment // ".emfm"
        call json%get(trim(tmp), var, found)
        if (found .eqv. .false.) call error_msg("Bad value in simulation_parameters segment")
        this%emfm = var

        tmp = segment // ".i_sphere"
        call json%get(trim(tmp), var_i, found)
        if (found .eqv. .false.) call error_msg("Bad value in simulation_parameters segment")
        this%to_radial_index_sphere = var_i

        tmp = segment // ".i_sphere_up"
        call json%get(trim(tmp), var_i, found)
        if (found .eqv. .false.) call error_msg("Bad value in simulation_parameters segment")
        this%from_radial_index_sphere = var_i

        tmp = segment // ".no_move_layer"
        call json%get(trim(tmp), var_i, found)
        if (found .eqv. .false.) call error_msg("Bad value in simulation_parameters segment")
        this%no_move_layer = var_i


        tmp = segment // ".no_linear_visc.flag"
        call json%get(trim(tmp), var_l, found)
        if (found .eqv. .false.) then
            call warn_msg("No flag found in NO_XL_VISC, we assume you do want it and act as it is true")
            var_l = .true.
        end if
        this%no_xl_flag = var_l
        if (var_l .eqv. .true.) then
            tmp = segment // ".no_linear_visc.start_layer"
            call json%get(trim(tmp), var_i, found)
            if (found .eqv. .false.) call error_msg("Parser: Bad value in simulation parameter no_xl_visc")
            this%start_no_xl_visc = var_i

            tmp = segment // ".no_linear_visc.end_layer"
            call json%get(trim(tmp), var_i, found)
            if (found .eqv. .false.) call error_msg("Parser: Bad value in simulation parameter no_xl_visc")
            this%end_no_xl_visc = var_i

            tmp = segment // ".no_linear_visc.offset"
            call json%get(trim(tmp), var_i, found)
            if (found .eqv. .false.) call error_msg("Parser: Bad value in simulation parameter no_xl_visc")
            this%offset_no_xl_visc = var_i
        end if

    end subroutine Parse_simulation_parameters

    subroutine Parse_rezone_advect(this, json)
        implicit none
        class (datafile_t) , intent(inout) :: this
        type(json_file)    , intent(inout) :: json

        character(len=*), parameter    :: segment  = "rezone_advect"
        character(len=255) :: tmp
        character(len=255) :: root
        logical :: found, var_l
        integer :: var

        tmp = segment // ".shorter_advect"
        call json%get(trim(tmp), var_l, found)
        if (found .eqv. .false.) call error_msg("Bad value in rezone_advect segment")
        this%shorter_advect = var_l

        tmp = segment // ".fix_overflow"
        call json%get(trim(tmp), var_l, found)
        if (found .eqv. .false.) call error_msg("Bad value in rezone_advect segment")
        this%fix_overflow = var_l

        tmp = segment // ".line_calc"
        call json%get(trim(tmp), var_l, found)
        if (found .eqv. .false.) call error_msg("Bad value in rezone_advect segment")
        this%line_calc = var_l

        tmp = segment // ".rezone_type"
        call json%get(trim(tmp), var, found)
        if (found .eqv. .false.) call error_msg("Bad value in rezone_advect segment")
        this%rezone_type = var

    end subroutine Parse_rezone_advect

    subroutine Parse_mesh(this, json)
        implicit none
        class (datafile_t) , intent(inout) :: this
        type(json_file)    , intent(inout) :: json

        character(len=*), parameter    :: segment  = "cell_set"
        character(len=255) :: tmp
        character(len=255) :: root
        logical :: found, var_l
        integer :: var

        tmp = segment // ".mesh_type"
        call json%get(trim(tmp), var, found)
        if (found .eqv. .false.) call error_msg("Bad value in mesh segment")
        this%mesh_type = var


        if (this%dimension == 2) then
            allocate(this%boundary_conditions(4))
            if (this%mesh_type == 2) this%boundary_conditions = (/2,2,2,2/)
        else
            allocate(this%boundary_conditions(6))
            if (this%mesh_type == 2) this%boundary_conditions = (/2,2,2,2,2,2/)
            if (this%mesh_type == 1) this%boundary_conditions = (/2,2,2,2,3,0/)
        end if

    end subroutine Parse_mesh

    subroutine Parse_layers(this, json)
        use geometry_module, only : Create_start_index_layer
        implicit none
        class (datafile_t) , intent(inout) :: this
        type(json_file)    , intent(inout) :: json

        character(len=*), parameter    :: segment  = "layers_materials"
        character(len=*), dimension(3), parameter :: layers = (/".number_cells_i", ".number_cells_j", ".number_cells_k"/)
        character(len=255) :: tmp
        character(len=255) :: root
        logical :: found, var_l
        integer :: var_i, num_cells,i,j
        real(8) :: var_r

        root = segment // layers(1)
        call json%info(trim(root), n_children=num_cells, found=found)
        if (found .eqv. .false.)  call error_msg("Parser: Problem in layers material datafile")
        this%number_layers_i = num_cells
        allocate(this%number_cells_i(num_cells))
        do i = 1, num_cells
            tmp = trim(root) // "(" // int2str(i) // ")"
            call json%get(trim(tmp), var_i, found)
            if (found .eqv. .false.) call error_msg("Parser: Bad value in layers datafile")
            this%number_cells_i(i) = var_i
        end do
        this%virt_nx = sum(this%number_cells_i)
        this%virt_nxp = this%virt_nx + 1
        allocate(this%start_layer_index_r(this%number_layers_i + 1 , 1))
        this%start_layer_index_r = Create_start_index_layer(this%number_cells_i, this%number_layers_i, 1)

        root = segment // layers(2)
        call json%info(trim(root), n_children=num_cells, found=found)
        if (found .eqv. .false.)   call error_msg("Parser: Problem in layers material datafile")
        this%number_layers_j = num_cells
        allocate(this%number_cells_j(num_cells))
        do i = 1, num_cells
            tmp = trim(root) // "(" // int2str(i) // ")"
            call json%get(trim(tmp), var_i, found)
            if (found .eqv. .false.) call error_msg("Parser: Bad value in layers datafile")
            this%number_cells_j(i) = var_i
        end do
        this%virt_ny = sum(this%number_cells_j)
        this%virt_nyp = this%virt_ny + 1
        allocate(this%start_layer_index_t(this%number_layers_j + 1 , 1))
        this%start_layer_index_t = Create_start_index_layer(this%number_cells_j, this%number_layers_j,1)

        root = segment // layers(3)
        call json%info(trim(root), n_children=num_cells, found=found)
        if (found .eqv. .false.) then
            this%dimension = 2
            allocate(this%number_cells_k(1))
            this%number_cells_k(1) = 1
            this%number_layers_k = 1
            this%virt_nz  = 1
            this%virt_nzp = 1
        else
            this%dimension = 3
            this%number_layers_k = num_cells
            allocate(this%number_cells_k(num_cells))
            do i = 1, num_cells
                tmp = trim(root) // "(" // int2str(i) // ")"
                call json%get(trim(tmp), var_i, found)
                if (found .eqv. .false.)  call error_msg("Bad value in layers datafile")
                this%number_cells_k(i) = var_i
            end do
            this%virt_nz  = this%virt_nx
            this%virt_nzp = this%virt_nxp
            this%virt_nx = sum(this%number_cells_k)
            this%virt_nxp = this%virt_nx + 1
            allocate(this%start_layer_index_p(this%number_layers_k + 1 , 1))
            this%start_layer_index_p = Create_start_index_layer(this%number_cells_k, this%number_layers_k,1)
        end if
    end subroutine Parse_layers

    subroutine Parse_materials(this, json)
        implicit none
        class (datafile_t) , intent(inout) :: this
        type(json_file)    , intent(inout) :: json

        character(len=*), parameter    :: segment  = "layers_materials"
        character(len=255) :: tmp
        character(len=255) :: root, mat_segment
        character(len=80), dimension(:), allocatable :: witness_array
        character(len=:),allocatable :: material_name
        logical :: found, already_added
        integer :: var_i,i, j,k, num_mats, mat_ind, last_index
        real(8) :: var_r

        root = segment // ".materials"
        call json%info(trim(root), n_children=num_mats, found=found)
        if (found .eqv. .false.) call error_msg("Problem in material datafile")
        allocate(witness_array(num_mats))
        witness_array = ""

        this%n_materials = num_mats
        this%reduct_num_mat = 0
        allocate(this%mat_index(num_mats))

        allocate(this%mat_eos(num_mats))
        allocate(this%mat_atomic_mass(num_mats))
        allocate(this%mat_z(num_mats))
        allocate(this%mat_z2(num_mats))
        allocate(this%mat_rho_0(num_mats))
        allocate(this%mat_sie_0(num_mats))
        allocate(this%mat_gamma_gas(num_mats))

        j = 1
        mat_ind = 1
        last_index = 1
        do i = 1, num_mats

            already_added = .false.
            tmp = trim(root) // "(" // int2str(i) // ")"
            call json%get(trim(tmp), material_name, found)

            do j = 1, num_mats
                if (str_eqv(trim(witness_array(j)), trim(material_name))) then
                    if (str_eqv(trim(witness_array(j)), "VACCUM")) then
                        mat_ind = 0
                    else
                        mat_ind = this%mat_index(j)
                    end if
                    already_added = .true.
                end if

            end do

            if (already_added .eqv. .false.) then
                this%mat_index(i) = mat_ind
                witness_array(i) = trim(material_name)
                if (Str_eqv(trim(material_name), "Vaccum") .eqv. .false.) then
                    mat_ind = mat_ind + 1
                    this%reduct_num_mat = this%reduct_num_mat + 1
                    last_index = mat_ind
                end if

            else
                this%mat_index(i) = mat_ind
                mat_ind = last_index
            end if

            if (Str_eqv(trim(material_name), "Vaccum") .eqv. .true.) then
                this%mat_index(i) = 0
                cycle
            end if
            mat_segment = material_name

            tmp = trim(mat_segment) // ".eos_type"
            call json%get(trim(tmp), var_i, found)
            if (found .eqv. .false.) call error_msg("Parse Bad value for material")
            this%mat_eos(this%mat_index(i)) = var_i

            tmp = trim(mat_segment) // ".A"
            call json%get(trim(tmp), var_r, found)
            if (found .eqv. .false.) call error_msg("Parse Bad value for material")
            this%mat_atomic_mass(this%mat_index(i)) = var_r

            tmp = trim(mat_segment) // ".Z"
            call json%get(trim(tmp), var_r, found)
            if (found .eqv. .false.) call error_msg("Parse Bad value for material")
            this%mat_z(this%mat_index(i)) = var_r

            tmp = trim(mat_segment) // ".Z_2"
            call json%get(trim(tmp), var_r, found)
            if (found .eqv. .false.) call error_msg("Parse Bad value for material")
            this%mat_z2(this%mat_index(i)) = var_r

            tmp = trim(mat_segment) // ".rho_0"
            call json%get(trim(tmp), var_r, found)
            if (found .eqv. .false.) call error_msg("Parse Bad value for material")
            this%mat_rho_0(this%mat_index(i)) = var_r

            tmp = trim(mat_segment) // ".sie_0"
            call json%get(trim(tmp), var_r, found)
            if (found .eqv. .false.) call error_msg("Parse Bad value for material")
            this%mat_sie_0(this%mat_index(i)) = var_r

            tmp = trim(mat_segment) // ".gamma_gas"
            call json%get(trim(tmp), var_r, found)
            if (found .eqv. .false.) call error_msg("Parse Bad value for material")
            this%mat_gamma_gas(this%mat_index(i)) = var_r

        end do

        deallocate(witness_array)

    end subroutine Parse_materials

    subroutine Parse_contours(this, json)
        implicit none
        class (datafile_t) , intent(inout) :: this
        type(json_file)    , intent(inout) :: json

        character(len=*), parameter    :: segment  = "contours"
        character(len=*), dimension(3), parameter :: contours = (/".contours_i", ".contours_j", ".contours_k"/)
        character(len=255) :: tmp
        character(len=255) :: root, cntr_root
        logical :: found, var_l
        integer :: var_i, num_cntrs,i,j, cntr_type
        real(8) :: var_r


        root = segment // contours(1)
        call json%info(trim(root), n_children=num_cntrs, found=found)
        if (found .eqv. .false.) call error_msg("Problem in contours datafile")
        allocate(this%contour_i_type(num_cntrs))
        allocate(this%contour_i_line(num_cntrs, 11))
        do i = 1, num_cntrs
            cntr_root = trim(root) // "(" // int2str(i) // ")"
            tmp = trim(cntr_root) // ".contour_type"
            call json%get(trim(tmp),cntr_type, found)
            if (found .eqv. .false.) call error_msg("Error in countour in datafile")
            this%contour_i_type(i) = cntr_type
            select case(cntr_type)
                case (0)
                    tmp = trim(cntr_root) // ".y1"
                    call json%get(trim(tmp), var_r, found)
                    if (found .eqv. .false.) call error_msg("Error in contour i in datafile")
                    this%contour_i_line(i, 1) = var_r

                    tmp = trim(cntr_root) // ".x1"
                    call json%get(trim(tmp), var_r, found)
                    if (found .eqv. .false.) call error_msg("Error in contour i in datafile")
                    this%contour_i_line(i, 2) = var_r

                    tmp = trim(cntr_root) // ".y2"
                    call json%get(trim(tmp), var_r, found)
                    if (found .eqv. .false.) call error_msg("Error in contour i in datafile")
                    this%contour_i_line(i, 3) = var_r

                    tmp = trim(cntr_root) // ".x2"
                    call json%get(trim(tmp), var_r, found)
                    if (found .eqv. .false.) call error_msg("Error in contour i in datafile")
                    this%contour_i_line(i, 4) = var_r

                case(1)
                    tmp = trim(cntr_root) // ".y1"
                    call json%get(trim(tmp), var_r, found)
                    if (found .eqv. .false.) call error_msg("Error in contour i in datafile")
                    this%contour_i_line(i, 1) = var_r

                    tmp = trim(cntr_root) // ".x1"
                    call json%get(trim(tmp), var_r, found)
                    if (found .eqv. .false.) call error_msg("Error in contour i in datafile")
                    this%contour_i_line(i, 2) = var_r

                    tmp = trim(cntr_root) // ".y2"
                    call json%get(trim(tmp), var_r, found)
                    if (found .eqv. .false.) call error_msg("Error in contour i in datafile")
                    this%contour_i_line(i, 3) = var_r

                    tmp = trim(cntr_root) // ".x2"
                    call json%get(trim(tmp), var_r, found)
                    if (found .eqv. .false.) call error_msg("Error in contour i in datafile")
                    this%contour_i_line(i, 4) = var_r

                    tmp = trim(cntr_root) // ".a"
                    call json%get(trim(tmp), var_r, found)
                    if (found .eqv. .false.) call error_msg("Error in contour i in datafile")
                    this%contour_i_line(i, 5) = var_r

                    tmp = trim(cntr_root) // ".b"
                    call json%get(trim(tmp), var_r, found)
                    if (found .eqv. .false.) call error_msg("Error in contour i in datafile")
                    this%contour_i_line(i, 6) = var_r

                    tmp = trim(cntr_root) // ".c"
                    call json%get(trim(tmp), var_r, found)
                    if (found .eqv. .false.) call error_msg("Error in contour i in datafile")
                    this%contour_i_line(i, 7) = var_r

                    tmp = trim(cntr_root) // ".d"
                    call json%get(trim(tmp), var_r, found)
                    if (found .eqv. .false.) call error_msg("Error in contour i in datafile")
                    this%contour_i_line(i, 8) = var_r

                     tmp = trim(cntr_root) // ".e"
                    call json%get(trim(tmp), var_r, found)
                    if (found .eqv. .false.) call error_msg("Error in contour i in datafile")
                    this%contour_i_line(i, 9) = var_r

                    tmp = trim(cntr_root) // ".f"
                    call json%get(trim(tmp), var_r, found)
                    if (found .eqv. .false.) call error_msg("Error in contour i in datafile")
                    this%contour_i_line(i, 10) = var_r
                case default
                    call error_msg("Unrecognized contour type")
            end select
        end do

        root = segment // contours(2)
        call json%info(trim(root), n_children=num_cntrs, found=found)
        if (found .eqv. .false.) call error_msg("Problem in contours datafile")
        allocate(this%contour_j_type(num_cntrs))
        allocate(this%theta0(num_cntrs, 1))
        do i = 1, num_cntrs
            cntr_root = trim(root) // "(" // int2str(i) // ")"
            tmp = trim(cntr_root) // ".units"
            call json%get(trim(tmp), cntr_type, found)
            if (found .eqv. .false.) call error_msg("Error in countour j in datafile")
            this%contour_j_type = cntr_type

            tmp = trim(cntr_root) // ".theta0"
            call json%get(trim(tmp), var_r, found)
            if (found .eqv. .false.) call error_msg("Error in countour j in datafile")
            this%theta0(i, 1) = var_r
        end do

        root = segment // contours(3)
        call json%info(trim(root), n_children=num_cntrs, found=found)
        if (found .eqv. .true. .and. this%dimension == 2) call error_msg("Problem should be 3d, but it appears to have no contour k")

        if (found .eqv. .false.) return

        allocate(this%contour_k_type(num_cntrs))
        allocate(this%phi0(num_cntrs, 1))
        do i = 1, num_cntrs
            cntr_root = trim(root) // "(" // int2str(i) // ")"
            tmp = trim(cntr_root) // ".units"
            call json%get(trim(tmp), cntr_type, found)
            if (found .eqv. .false.) call error_msg("Error in countour in datafile")
            this%contour_j_type = cntr_type

            tmp = trim(cntr_root) // ".phi0"
            call json%get(trim(tmp), var_r, found)
            if (found .eqv. .false.) call error_msg("Error in countour in datafile")
            this%phi0(i, 1) = var_r
        end do
    end subroutine Parse_contours

    subroutine Parse_zones(this, json)
        implicit none
        class (datafile_t) , intent(inout) :: this
        type(json_file)    , intent(inout) :: json

        character(len=*), parameter    :: segment  = "zone"
        character(len=*), dimension(3), parameter :: zones = (/".zone_i", ".zone_j", ".zone_k"/)
        character(len=255) :: tmp
        character(len=255) :: root, zone_root
        logical :: found, var_l
        integer :: var_i, num_zones,i,j, zone_type
        real(8) :: var_r


        root = segment // zones(1)
        call json%info(trim(root), n_children=num_zones, found=found)
        if (found .eqv. .false.) call error_msg("Problem in zones datafile")
        allocate(this%zone_i_type(num_zones))
        allocate(this%zone_i_d(num_zones))

        do i = 1, num_zones
            zone_root = trim(root) // "(" // int2str(i) // ")"
            tmp = trim(zone_root) // ".type"
            call json%get(trim(tmp), zone_type, found)
            if (found .eqv. .false.) call error_msg("Error in zones in datafile")
            this%zone_i_type(i) = zone_type

            tmp = trim(zone_root) // ".dr"
            call json%get(trim(tmp), var_r, found)
            if (found .eqv. .false.) call error_msg("Error in zones in datafile")
            this%zone_i_d(i) = var_r
        end do

        root = segment // zones(2)
        call json%info(trim(root), n_children=num_zones, found=found)
        if (found .eqv. .false.) call error_msg("Problem in zones datafile")
        allocate(this%zone_j_type(num_zones))
        allocate(this%zone_j_d(num_zones))

        do i = 1, num_zones
            zone_root = trim(root) // "(" // int2str(i) // ")"
            tmp = trim(zone_root) // ".type"
            call json%get(trim(tmp), zone_type, found)
            if (found .eqv. .false.) call error_msg("Error in zones in datafile")
            this%zone_j_type(i) = zone_type

            tmp = trim(zone_root) // ".d_theta"
            call json%get(trim(tmp), var_r, found)
            if (found .eqv. .false.) call error_msg("Error in zones in datafile")
            this%zone_j_d(i) = var_r
        end do

        root = segment // zones(3)
        call json%info(trim(root), n_children=num_zones, found=found)
        if (found .eqv. .false.) return
        allocate(this%zone_k_type(num_zones))
        allocate(this%zone_k_d(num_zones))

        do i = 1, num_zones
            zone_root = trim(root) // "(" // int2str(i) // ")"
            tmp = trim(zone_root) // ".type"
            call json%get(trim(tmp), zone_type, found)
            if (found .eqv. .false.) call error_msg("Error in zones in datafile")
            this%zone_k_type(i) = zone_type

            tmp = trim(zone_root) // ".d_phi"
            call json%get(trim(tmp), var_r, found)
            if (found .eqv. .false.) call error_msg("Error in zones in datafile")
            this%zone_k_d(i) = var_r
        end do

    end subroutine Parse_zones


    subroutine Parse_checkpoint_restart(this, json)
        implicit none
        class (datafile_t) , intent(inout) :: this
        type(json_file)    , intent(inout) :: json

        character(len=*), parameter    :: segment  = "checkpoint_restart"
        character(len=255) :: tmp
        character(len=255) :: root
        character(len=:),allocatable :: parallel_param
        logical :: found
        integer :: var


    end subroutine Parse_checkpoint_restart


end module datafile_module
