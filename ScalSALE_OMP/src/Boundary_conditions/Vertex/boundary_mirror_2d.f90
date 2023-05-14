
module boundary_mirror_2d_module
    use vertex_boundary_condition_module, only : vertex_boundary_condition_t
    use data_module, only : data_t
    use quantity_module, only : quantity_t
    implicit none
    private

    type, extends(vertex_boundary_condition_t), public :: boundary_mirror_2d_t
        private
    contains

        procedure, public :: Calculate => Boundary_mirror_2d_calculate
        procedure, public :: Calculate_face
        procedure, public :: Calculate_corner

    end type boundary_mirror_2d_t



contains
    subroutine Boundary_mirror_2d_calculate (this, v_quantity, coordinates, angle, edge_num)
        use geometry_module, only : Mirror_image
        class (boundary_mirror_2d_t)           , intent (in out) :: this      

        class(quantity_t), intent (in out) :: v_quantity  
        type(data_t), dimension(:), pointer,intent(inout)        :: coordinates
        real (8)                           , intent (in)     :: angle    
        integer                            , intent (in)     :: edge_num 

        real(8), dimension(:,:,:), pointer :: x, y
        integer :: nxp, nyp, ii,jj,i1,i2,j1,j2

        nxp = v_quantity%d1
        nyp = v_quantity%d2
        call coordinates(1)%point_to_data (x)
        call coordinates(2)%point_to_data (y)

        do ii = 1, nxp
            if (ii < nxp) then
                i1 = ii
                i2 = ii + 1
            else
                i1 = ii - 1
                i2 = ii
            end if

            call Mirror_image(x(i1, 1, 1), y(i1, 1, 1), &
                x(i2, 1, 1), y(i2, 1, 1), &
                x(ii, 2, 1), y(ii, 2, 1),&
                x(ii, 0, 1), y(ii, 0, 1))
        end do
        do ii = 1, nxp
            if (ii < nxp) then
                i1 = ii
                i2 = ii + 1
            else
                i1 = ii - 1
                i2 = ii
            end if
            call Mirror_image(x(i1, nyp, 1), y(i1, nyp, 1), &
                x(i2, nyp, 1), y(i2, nyp, 1), &
                x(ii, nyp-1, 1), y(ii, nyp-1, 1), &
                x(ii, nyp+1, 1), y(ii, nyp+1, 1))
        end do
        do jj = 0, nyp + 1
            if (jj < nyp) then
                j1 = jj
                j2 = jj + 1
            else
                j1 = jj - 1
                j2 = jj
            end if
            call Mirror_image(x(nxp, j1, 1), y(nxp, j1, 1), &
                x(nxp, j2, 1), y(nxp, j2, 1), &
                x(nxp-1, jj, 1), y(nxp-1, jj, 1), &
                x(nxp+1, jj, 1), y(nxp+1, jj, 1))
        end do


        do jj = 0, nyp + 1
            if (jj < nyp) then
                j1 = jj
                j2 = jj + 1
            else
                j1 = jj - 1
                j2 = jj
            end if
            call Mirror_image(x(1, j1, 1), y(1, j1, 1), &
                x(1, j2, 1), y(1, j2, 1), &
                x(2, jj, 1), y(2, jj, 1), &
                x(0, jj, 1), y(0, jj, 1))
        end do

    end subroutine Boundary_mirror_2d_calculate

    subroutine Calculate_face (this, v_quantity, angle, edge_num)
        class (boundary_mirror_2d_t), intent (in out)     :: this      

        class(quantity_t), intent (in out) :: v_quantity  
        real (8)                           , intent (in)     :: angle    
        integer                            , intent (in)     :: edge_num 

    end subroutine

    subroutine Calculate_corner (this, v_quantity, angle, edge_num)
        implicit none
        class (boundary_mirror_2d_t), intent (in out)     :: this      

        class(quantity_t), intent (in out) :: v_quantity  
        real (8)                           , intent (in)     :: angle    
        integer                            , intent (in)     :: edge_num 

    end subroutine
end module Boundary_mirror_2d_module
