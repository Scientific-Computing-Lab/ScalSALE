
module replace_words_module
    use json_module
    use, intrinsic :: iso_fortran_env , only: error_unit, output_unit
    use general_utils_module, only: str_eqv, int2str, error_msg, lower
    implicit none

contains

    subroutine Replace_words(json)
        implicit none
        type(json_file), intent(inout) :: json

        call Replace_zone(json)
        call Replace_contour(json)
        call Replace_advect(json)
        call Replace_mesh(json)
        call Replace_material(json)

    end subroutine Replace_words

    subroutine Replace_zone(json)
        implicit none
        type(json_file), intent(inout) :: json
        character(len=*), parameter    :: zone  = "zone"
        character(len=*), dimension(3), parameter :: zones = (/".zone_i", ".zone_j", ".zone_k"/)
        character(len=255) :: tmp
        character(len=255) :: root
        character(len=:),allocatable :: value_str
        integer :: num_zones
        logical :: found
        integer :: i, j

        do i = 1, 3
            tmp = zone // zones(i)
            root = tmp
            call json%info(trim(tmp), n_children=num_zones, found=found)
            if (found .eqv. .false.) cycle
            do j = 1, num_zones
                tmp = root
                tmp = trim(tmp) // "(" // int2str(j) // ").type"
                call json%get(trim(tmp) , value_str, found)
                if (found .eqv. .false.) then
                    call error_msg("Error no type for zone or zone type without string")
                end if
                if (str_eqv(value_str, "constant")) then
                    call json%update(trim(tmp), 0, found)
                elseif (str_eqv(value_str, "geom_first") .or. str_eqv(value_str, "geometry_first")) then
                    call json%update(trim(tmp), 2, found)
                elseif (str_eqv(value_str, "geom_last") .or.     str_eqv(value_str, "geometry_last")) then
                    call json%update(trim(tmp), 3, found)
                end if
            end do
        end do
    end subroutine Replace_zone


    subroutine Replace_advect(json)
        implicit none
        type(json_file), intent(inout) :: json

        character(len=*), parameter    :: segment  = "rezone_advect"
        character(len=255) :: tmp
        character(len=:),allocatable :: value_str
        integer :: num_zones
        logical :: found
        integer :: i, j

        tmp = segment // ".rezone_type"
        call json%get(trim(tmp) , value_str, found)
        if (str_eqv(value_str, "lagrange")) then
            call json%update(trim(tmp), 0, found)
        elseif (str_eqv(value_str, "euler")) then
            call json%update(trim(tmp), 1, found)
        else if (str_eqv(value_str, "ale")) then
            call json%update(trim(tmp), 2, found)
        end if
    end subroutine Replace_advect

    subroutine Replace_mesh(json)
        implicit none
        type(json_file), intent(inout) :: json

        character(len=*), parameter    :: segment  = "cell_set"
        character(len=255) :: tmp
        character(len=:),allocatable :: value_str
        logical :: found

        tmp = segment // ".mesh_type"
        call json%get(trim(tmp) , value_str, found)
        if (str_eqv(value_str, "x_y") .or. str_eqv(value_str, "x_y_z")) then
            call json%update(trim(tmp), 2, found)
        elseif (str_eqv(value_str, "x y") .or. str_eqv(value_str, "x y z")) then
            call json%update(trim(tmp), 2, found)
        elseif (str_eqv(value_str, "pyramid")) then
            call json%update(trim(tmp), 1, found)
        end if

    end subroutine Replace_mesh

    subroutine Replace_contour(json)
        implicit none
        type(json_file), intent(inout) :: json

        character(len=*), parameter    :: segment  = "contours"
        character(len=*), dimension(3), parameter :: contours = (/".contours_i", ".contours_j", ".contours_k"/)
        character(len=255) :: tmp
        character(len=255) :: root, cntr_root
        character(len=:),allocatable :: value_str
        integer :: num_cntrs
        logical :: found

        integer :: i, j
        do i = 1, 3
            tmp = segment // contours(i)
            root = tmp
            call json%info(trim(tmp), n_children=num_cntrs, found=found)
            if (found .eqv. .false.) cycle
            do j = 1, num_cntrs
                tmp = root
                cntr_root = trim(root) // "(" // int2str(j) // ")"
                tmp = trim(cntr_root) // ".contour_type"
                call json%get(trim(tmp), value_str, found)
                if (found .eqv. .true.) then
                    if (str_eqv(value_str, "line") .or. str_eqv(value_str, "kav")) then
                        call json%update(trim(tmp), 0, found)
                    elseif (str_eqv(value_str, "super") .or. str_eqv(value_str, "elipse")) then
                        call json%update(trim(tmp), 1, found)
                        tmp = trim(cntr_root) // ".x2"
                        call json%get(trim(tmp), value_str, found)
                        select case(lower(value_str))
                            case ("radian")
                                call json%update(trim(tmp), 0, found)
                            case ("pi")
                                call json%update(trim(tmp), 1, found)
                            case ("degree")
                                call json%update(trim(tmp), 2, found)
                            case default
                                call error_msg("Replace words: Unrecognized option in x2 in conoutrs")
                        end select
                        tmp = trim(cntr_root) // ".y2"
                        call json%get(trim(tmp), value_str, found)
                        select case(lower(value_str))
                            case ("teta")
                                call json%update(trim(tmp), 0, found)
                            case default
                                call error_msg("Replace words: No other option than teta in y2 in conoutrs")
                        end select
                    end if
                end if

                tmp = root
                tmp = trim(tmp) // "(" // int2str(j) // ").units"
                call json%get(trim(tmp), value_str, found)
                if (found .eqv. .true.) then
                    select case(lower(value_str))
                        case ("regular")
                            call json%update(trim(tmp), 0, found)
                        case ("pi")
                            call json%update(trim(tmp), 1, found)
                        case ("degree")
                            call json%update(trim(tmp), 2, found)
                        case default
                            call error_msg("Replace words: bad value in units j")
                    end select
                end if
            end do
        end do

        if (json%failed()) then
            call json%print_error_message(error_unit)
        end if
    end subroutine Replace_contour

    subroutine Replace_material(json)
        implicit none
        type(json_file), intent(inout) :: json

        character(len=*), parameter    :: segment  = "layers_materials"
        character(len=255) :: tmp
        character(len=255) :: root
        character(len=:),allocatable :: value_str, material
        logical :: found
        integer :: num_mats,i

        tmp = segment // ".materials"
        call json%info(trim(tmp), n_children=num_mats, found=found)
        root = tmp
        do i = 1, num_mats
            tmp = root
            tmp = trim(tmp) // "(" // int2str(i) // ")"
            call json%get(trim(tmp) , material, found)
            if (found .eqv. .false. .or. Str_eqv(material, "Vaccum")) cycle
            call json%get("$"//trim(material) // ".eos_type",value_str, found)

            if (str_eqv(value_str, "0")) cycle
            select case(value_str)
                case ('ideal')
                    call json%update("$"//trim(material) // ".eos_type", 0, found)
                case ('ideal_gas')
                    call json%update("$"//trim(material) // ".eos_type", 0, found)
                case ('ideal gas')
                    call json%update("$"//trim(material) // ".eos_type", 0, found)
                case default
                    call error_msg("Replace words: Unrecognized EoS in datafile")
                    stop
            end select
        end do
    end subroutine Replace_material
end module replace_words_module
