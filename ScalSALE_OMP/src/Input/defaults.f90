
module defaults_module
    use json_module
    use, intrinsic :: iso_fortran_env , only: error_unit, output_unit
    use general_utils_module, only: str_eqv, int2str, error_msg

    implicit none

contains

    subroutine Set_defaults(json)
        implicit none
        type(json_file)    , intent(inout) :: json
        call Default_switches(json)
        call Default_simulation_parameters(json)
        call Default_rezone_advect(json)
        call Default_parallel(json)
        call Default_checkpoint_restart(json)
        call Default_material(json)
        call Default_mesh(json)
    end subroutine Set_defaults


    subroutine Default_switches(json)
        implicit none
        type(json_file)    , intent(inout) :: json

        character(len=*), parameter    :: segment  = "switches"
        character(len=255) :: tmp
        character(len=255) :: root
        character(len=:),allocatable :: value_str
        type(json_value) , pointer :: p, inp
        type(json_core)            :: json_creator
        real(8) :: value_r
        integer :: value_i
        logical :: found

        tmp = segment
        call json%info(trim(tmp), found=found)
        if (found .eqv. .false.) then
            call json_creator%create_object(inp, 'switches')
            call json%add("$switches", inp)
        end if

        tmp = segment // ".sw_wilkins"
        call json%info(trim(tmp), found=found)
        if (found .eqv. .false.) then
            call json%add(trim(tmp), 1, found=found)
        end if

        tmp = segment // ".sw_cr"
        call json%info(trim(tmp), found=found)
        if (found .eqv. .false.) call json%add(trim(tmp), 0, found=found)

        return
    end subroutine Default_switches

    subroutine Default_simulation_parameters(json)
        use general_utils_module  ,only : Str2dble
        implicit none

        type(json_file)    , intent(inout) :: json

        character(len=*), parameter    :: segment  = "simulation_parameters"
        character(len=255) :: tmp
        character(len=255) :: root
        character(len=:),allocatable :: value_str
        type(json_value) , pointer :: p, inp
        type(json_core)            :: json_creator
        logical :: found
        real(8) :: val


        tmp = segment
        call json%info(trim(tmp), found=found)
        if (found .eqv. .false.) then
            call json_creator%create_object(inp, segment)
            call json%add("$simulation_parameters", inp)
        end if

        tmp = segment // ".time_final"
        val = 30e-6
        call json%info(trim(tmp), found=found)
        if (found .eqv. .false.) call json%add(trim(tmp), val, found=found)

        tmp = segment // ".dt0"
        val = 1e-14
        call json%info(trim(tmp), found=found)
        if (found .eqv. .false.) call json%add(trim(tmp), Str2dble("1e-14"), found=found)

        tmp = segment // ".dt"
        call json%info(trim(tmp), found=found)
        if (found .eqv. .false.) call json%add(trim(tmp), Str2dble("1e-14"), found=found)

        tmp = segment // ".dt_factor"
        call json%info(trim(tmp), found=found)
        if (found .eqv. .false.) call json%add(trim(tmp), Str2dble("1e-1"), found=found)

        tmp = segment // ".emf"
        call json%info(trim(tmp), found=found)
        if (found .eqv. .false.) call json%add(trim(tmp), Str2dble("1e-5"), found=found)

        tmp = segment // ".cyl"
        call json%info(trim(tmp), found=found)
        if (found .eqv. .false.) call json%add(trim(tmp), 1, found=found)

        tmp = segment // ".emfm"
        call json%info(trim(tmp), found=found)
        if (found .eqv. .false.) call json%add(trim(tmp), Str2dble("1e-20"), found=found)

        tmp = segment // ".linear_visc_fac"
        call json%info(trim(tmp), found=found)
        if (found .eqv. .false.) call json%add(trim(tmp), Str2dble("0.15"), found=found)

        tmp = segment // ".quad_visc_fac"
        call json%info(trim(tmp), found=found)
        if (found .eqv. .false.) call json%add(trim(tmp), Str2dble("4.0"), found=found)

        tmp = segment // ".dt_cour_fac"
        call json%info(trim(tmp), found=found)
        if (found .eqv. .false.) call json%add(trim(tmp), Str2dble("3.0"), found=found)

        tmp = segment // ".i_sphere"
        call json%info(trim(tmp), found=found)
        if (found .eqv. .false.) call json%add(trim(tmp), -1, found=found)

        tmp = segment // ".i_sphere_up"
        call json%info(trim(tmp), found=found)
        if (found .eqv. .false.) call json%add(trim(tmp), -1, found=found)

        tmp = segment // ".no_move_layer"
        call json%info(trim(tmp), found=found)
        if (found .eqv. .false.) call json%add(trim(tmp), -1, found=found)

        tmp = segment // ".no_linear_visc"
        call json%info(trim(tmp), found)
        if (found .eqv. .false.) then
            call json_creator%create_object(inp, segment)
            call json%add("$" // trim(tmp), inp)
            call json%add(trim(tmp) // ".flag"       , .false., found=found)
            call json%add(trim(tmp) // ".start_layer", 1, found=found)
            call json%add(trim(tmp) // ".end_layer"  , 1, found=found)
            call json%add(trim(tmp) // ".offset"     , 0, found=found)
        end if

        tmp = segment // ".min_pressure_shock"
        call json%info(trim(tmp), found)
        if (found .eqv. .false.) call json%add(trim(tmp), Str2dble("1e9"), found=found)

        return
    end subroutine Default_simulation_parameters

    subroutine Default_parallel(json)
        implicit none
        type(json_file)    , intent(inout) :: json

        character(len=*), parameter    :: segment  = "parallel"
        character(len=255) :: tmp
        character(len=255) :: root
        character(len=:),allocatable :: value_str
        type(json_value) , pointer :: p, inp
        type(json_core)            :: json_creator
        logical :: found_np, found_npx, found_npy, found_npz, found
        integer :: value_i


        tmp = segment
        call json%info(trim(tmp), found=found)
        if (found .eqv. .false.) then
            call json_creator%create_object(inp, segment)
            call json%add("$parallel", inp)
        end if

        tmp = segment // ".np"
        call json%info(trim(tmp), found=found)
        if (found .eqv. .false.) call json%add(trim(tmp), 1, found=found)

        tmp = segment // ".npx"
        call json%info(trim(tmp), found=found)
        if (found .eqv. .false.) call json%add(trim(tmp), 1, found=found)

        tmp = segment // ".npy"
        call json%info(trim(tmp), found=found)
        if (found .eqv. .false.) call json%add(trim(tmp), 1, found=found)

        tmp = segment // ".npz"
        call json%info(trim(tmp), found=found)
        if (found .eqv. .false.) call json%add(trim(tmp), 1, found=found)

        tmp = segment // ".threads"
        call json%info(trim(tmp), found=found)
        if (found .eqv. .false.) call json%add(trim(tmp), 1, found=found)

        return
    end subroutine Default_parallel

    subroutine Default_rezone_advect(json)
        implicit none
        type(json_file)    , intent(inout) :: json

        character(len=*), parameter    :: segment  = "rezone_advect"
        character(len=255) :: tmp
        character(len=255) :: root
        character(len=:),allocatable :: value_str
        type(json_value) , pointer :: p, inp
        type(json_core)            :: json_creator
        logical :: found


        tmp = segment
        call json%info(trim(tmp), found=found)
        if (found .eqv. .false.) then
            call json_creator%create_object(inp, segment)
            call json%add("$rezone_advect", inp)
        end if

        tmp = segment // ".shorter_advect"
        call json%info(trim(tmp), found=found)
        if (found .eqv. .false.) call json%add(trim(tmp), .true., found=found)

        tmp = segment // ".fix_overflow"
        call json%info(trim(tmp), found=found)
        if (found .eqv. .false.) call json%add(trim(tmp), .true., found=found)

        tmp = segment // ".rezone_type"
        call json%info(trim(tmp), found=found)
        if (found .eqv. .false.) call json%add(trim(tmp), 0, found=found)

        tmp = segment // ".line_calc"
        call json%info(trim(tmp), found=found)
        if (found .eqv. .false.) call json%add(trim(tmp), .false., found=found)

        return
    end subroutine Default_rezone_advect

    subroutine Default_checkpoint_restart(json)
        implicit none
        type(json_file)    , intent(inout) :: json

        character(len=*), parameter    :: segment  = "checkpoint_restart"
        character(len=255) :: tmp
        character(len=255) :: root
        character(len=:),allocatable :: value_str
        type(json_value) , pointer :: p, inp
        type(json_core)            :: json_creator
        logical :: found


        tmp = segment
        call json%info(trim(tmp), found=found)
        if (found .eqv. .false.) then
            call json_creator%create_object(inp, segment)
            call json%add("$checkpoint_restart", inp)
        end if

        tmp = segment // ".run_name"
        call json%info(trim(tmp), found=found)
        if (found .eqv. .false.) call json%add(trim(tmp), "", found=found)

        tmp = segment // ".scr_prefix"
        call json%info(trim(tmp), found=found)
        if (found .eqv. .false.) call json%add(trim(tmp), "/home/", found=found)

        tmp = segment // ".checkpoint_seconds"
        call json%info(trim(tmp), found=found)
        if (found .eqv. .false.) call json%add(trim(tmp), 60.0, found=found)

        tmp = segment // ".checkpoint_overhead"
        call json%info(trim(tmp), found=found)
        if (found .eqv. .false.) call json%add(trim(tmp), 0.0, found=found)

        return
    end subroutine Default_checkpoint_restart

    subroutine Default_material(json)
        implicit none
        type(json_file), intent(inout) :: json

        character(len=*), parameter    :: segment  = "layers_materials"
        character(len=255) :: tmp
        character(len=255) :: root
        character(len=:),allocatable :: value_str, material
        logical :: found
        integer :: num_mats,i
        real(8) :: var_r

        tmp = segment // ".materials"
        call json%info(trim(tmp), n_children=num_mats, found=found)
        root = tmp
        do i = 1, num_mats

            tmp = root
            tmp = trim(tmp) // "(" // int2str(i) // ")"
            call json%get(trim(tmp) , material, found)
            if (found .eqv. .false.) cycle
            call json%get("$"//trim(material) // ".sie_0",var_r, found)

            if ( (found .eqv. .false.) .and. (Str_eqv(material, "Vaccum") .eqv. .false.) ) then
                call json%update("$"//trim(material) // ".sie_0", 0, found)
            end if
        end do
    end subroutine Default_material


    subroutine Default_mesh(json)
        implicit none
        type(json_file), intent(inout) :: json

        character(len=*), parameter    :: segment  = "cell_set"
        character(len=255) :: tmp
        character(len=255) :: root
        character(len=:),allocatable :: value_str, material
        logical :: found
        integer :: num_mats,i, var
        real(8) :: var_r



    end subroutine Default_mesh
end module defaults_module
