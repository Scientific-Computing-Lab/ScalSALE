
module boundary_mirror_3d_module
    use vertex_boundary_condition_module, only : vertex_boundary_condition_t
    use data_module, only : data_t
    use quantity_module, only : quantity_t
    implicit none
    private

    type, extends(vertex_boundary_condition_t), public :: boundary_mirror_3d_t
        private
    contains

        procedure, public :: Calculate => Boundary_mirror_3d_calculate
        procedure, public :: Calculate_face
        procedure, public :: Calculate_corner

    end type boundary_mirror_3d_t



contains
    subroutine Boundary_mirror_3d_calculate (this, v_quantity, coordinates, angle, edge_num)
        use geometry_module, only: Mirror_image_3d
        implicit none
        class (boundary_mirror_3d_t)           , intent (in out) :: this      

        class(quantity_t), intent (in out) :: v_quantity  
        type(data_t), dimension(:), pointer,intent(inout)        :: coordinates
        real (8)                           , intent (in)     :: angle    
        integer                            , intent (in)     :: edge_num 
        integer, save :: cyc = 0
        real(8), dimension(:, :, :), pointer :: vol1, vol2

        real(8), dimension(:, :, :), pointer :: x               
        real(8), dimension(:, :, :), pointer :: y               
        real(8), dimension(:, :, :), pointer :: z               
        real(8), dimension(:, :, :), pointer :: material_x      
        real(8), dimension(:, :, :), pointer :: material_y      
        real(8), dimension(:, :, :), pointer :: material_z      
        real(8), dimension(:, :, :), pointer :: mat_vol      

        real(8), dimension(3) :: normal_i1, normal_inxp,normal_j1, normal_jnyp,normal_k1, normal_knzp

        logical :: wall_x_top, wall_x_bot, wall_y_top, wall_y_bot, wall_z_top, wall_z_bot
        integer :: ii, jj, kk, nxp, nyp, nzp   

        call v_quantity%Point_to_data(x, y, z)

        wall_x_top = v_quantity%boundary_params%parallel_params%is_wall_x_top
        wall_x_bot = v_quantity%boundary_params%parallel_params%is_wall_x_bot
        wall_y_top = v_quantity%boundary_params%parallel_params%is_wall_y_top
        wall_y_bot = v_quantity%boundary_params%parallel_params%is_wall_y_bot
        wall_z_top = v_quantity%boundary_params%parallel_params%is_wall_z_top
        wall_z_bot = v_quantity%boundary_params%parallel_params%is_wall_z_bot

        normal_i1   = v_quantity%boundary_params%normal_i1
        normal_inxp = v_quantity%boundary_params%normal_inxp
        normal_j1   = v_quantity%boundary_params%normal_j1
        normal_jnyp = v_quantity%boundary_params%normal_jnyp
        normal_k1   = v_quantity%boundary_params%normal_k1
        normal_knzp = v_quantity%boundary_params%normal_knzp

        nxp = v_quantity%d1
        nyp = v_quantity%d2
        nzp = v_quantity%d3
        if (wall_x_bot .eqv. .true.) then
            ii = 0
            do kk = 0, nzp + 1
                do jj = 0, nyp + 1

                    call Mirror_image_3d(normal_i1(1), normal_i1(2)  , &
                        normal_i1(3), v_quantity%boundary_params%plane_const_i1, &
                        x(2, jj, kk), y(2, jj, kk), z(2, jj, kk), &
                        x(0, jj, kk), y(0, jj, kk), z(0, jj, kk))

                end do
            end do
        end if

        if (wall_x_top .eqv. .true.) then
            ii = nxp
            do kk = 0, nzp + 1
                do jj = 0, nyp + 1
                    call Mirror_image_3d(normal_inxp(1), normal_inxp(2)  , &
                        normal_inxp(3), v_quantity%boundary_params%plane_const_inxp, &
                        x(nxp - 1, jj, kk), y(nxp-1, jj, kk), z(nxp-1, jj, kk), &
                        x(nxp + 1, jj, kk), y(nxp+1, jj, kk), z(nxp+1, jj, kk))
                end do
            end do
        end if

        if (wall_y_bot .eqv. .true.) then
            jj = 0
            do kk = 0, nzp + 1
                do ii = 0, nxp + 1
                    call Mirror_image_3d(normal_j1(1), normal_j1(2)  , &
                        normal_j1(3), v_quantity%boundary_params%plane_const_j1, &
                        x(ii, 2, kk), y(ii, 2, kk), z(ii, 2, kk), &
                        x(ii, 0, kk), y(ii, 0, kk), z(ii, 0, kk))
                end do
            end do
        end if
        if (wall_y_top .eqv. .true.) then
            jj = nyp
            do kk = 0,nzp + 1
                do ii = 0, nxp + 1
                    call Mirror_image_3d(normal_jnyp(1), normal_jnyp(2)  , &
                        normal_jnyp(3), v_quantity%boundary_params%plane_const_jnyp, &
                        x(ii, nyp-1, kk), y(ii, nyp-1, kk), z(ii, nyp-1, kk), &
                        x(ii, nyp+1, kk), y(ii, nyp+1, kk), z(ii, nyp+1, kk))

                end do
            end do
        end if

        if (wall_z_bot .eqv. .true.) then
            kk = 0
            do jj = 0, nyp + 1
                do ii = 0, nxp + 1
                    call Mirror_image_3d(normal_k1(1), normal_k1(2)  , &
                        normal_k1(3), v_quantity%boundary_params%plane_const_k1, &
                        x(ii, jj, 2), y(ii, jj, 2), z(ii, jj, 2), &
                        x(ii, jj, 0), y(ii, jj, 0), z(ii, jj, 0))

                end do
            end do
        end if


        if (wall_z_top .eqv. .true.) then
            kk = nzp
            do jj = 0, nyp + 1
                do ii = 0, nxp + 1
                    call Mirror_image_3d(normal_knzp(1), normal_knzp(2)  , &
                        normal_knzp(3), v_quantity%boundary_params%plane_const_knzp, &
                        x(ii, jj, nzp-1), y(ii, jj, nzp-1), z(ii, jj, nzp-1), &
                        x(ii, jj, nzp+1), y(ii, jj, nzp+1), z(ii, jj, nzp+1))
                end do
            end do
        end if
        cyc = cyc + 1
    end subroutine Boundary_mirror_3d_calculate

    subroutine Calculate_face (this, v_quantity, angle, edge_num)
        class (boundary_mirror_3d_t), intent (in out)     :: this      

        class(quantity_t), intent (in out) :: v_quantity  
        real (8)                           , intent (in)     :: angle    
        integer                            , intent (in)     :: edge_num 

    end subroutine

    subroutine Calculate_corner (this, v_quantity, angle, edge_num)
        implicit none
        class (boundary_mirror_3d_t), intent (in out)     :: this      

        class(quantity_t), intent (in out) :: v_quantity  
        real (8)                           , intent (in)     :: angle    
        integer                            , intent (in)     :: edge_num 

    end subroutine
end module Boundary_mirror_3d_module
