
module slip_vertex_module
    use vertex_boundary_condition_module, only : vertex_boundary_condition_t
    use data_module, only : data_t
    use quantity_module, only : quantity_t
    implicit none
    private

    type, extends(vertex_boundary_condition_t), public :: slip_vertex_t
        private
    contains

        procedure, public  :: Calculate => Slip_vertex_calculate_2d
        procedure, public :: Calculate_face
        procedure, public :: Calculate_corner

    end type



contains
    subroutine Slip_vertex_calculate_2d (this, v_quantity, coordinates, angle, edge_num)
        class (slip_vertex_t)              , intent (in out) :: this      

        class(quantity_t), intent (in out) :: v_quantity  
        type(data_t), dimension(:), pointer,intent(inout)        :: coordinates
        real (8)                           , intent (in)     :: angle    
        integer                            , intent (in)     :: edge_num 

        real(8), dimension (:,:,:), pointer :: x 
        real(8), dimension (:,:,:), pointer :: y 

        real(8) :: slope        
        real(8) :: delta_x      
        real(8) :: delta_y      
        real(8) :: slope_factor 
        real(8) :: temp

        real(8), dimension (:,:,:), pointer :: values_dim1

        real(8), dimension (:,:,:), pointer :: values_dim2
        integer :: nxp, nyp, i,j

        nxp = v_quantity%d1
        nyp = v_quantity%d2

        call v_quantity%Point_to_data(values_dim1, values_dim2)

        call coordinates(1)%Point_to_data(x)
        call coordinates(2)%Point_to_data(y)

        select case (edge_num)

            case(1,2) 
                if (edge_num == 1) i = 1
                if (edge_num == 2) i = nxp

                do j = 1, nyp
                    slope = 1d20
                    delta_x = x(i, nyp, 1) - x(i, 1, 1)
                    if (abs(delta_x) > 1d-20) then
                        slope = (y(i, nyp, 1) - y(i, 1, 1)) / delta_x
                        slope_factor = 1d0 / (1d0 + slope * slope)

                        temp = values_dim1(i, j, 1)
                        values_dim1(i, j, 1) = (slope * values_dim2(i, j, 1) + temp) * slope_factor
                        values_dim2(i, j, 1) = (slope * slope * values_dim2(i, j, 1) + slope * temp) * slope_factor

                    else
                        values_dim1(i, j, 1) = 0d0
                    end if
                end do

            case(3, 4) 
                if (edge_num == 3) j = 1
                if (edge_num == 4) j = nyp

                do i = 1, nxp
                    slope = 1d20
                    delta_y = y(nxp, j, 1) - y(1, j, 1)
                    if (abs(delta_y) > 1d-20) then
                        slope = delta_y / (x(nxp, j, 1) - x(1, j, 1))
                        slope_factor = 1d0 / (1d0 + slope * slope)

                        temp = values_dim1(i, j, 1)
                        values_dim1(i, j, 1) = (slope * values_dim2(i, j, 1) + temp) * slope_factor
                        values_dim2(i, j, 1) = (slope * slope * values_dim2(i, j, 1) + slope * temp) * slope_factor

                    else
                        values_dim2(i, j, 1) = 0d0
                    end if
                end do
        end select
    end subroutine


    subroutine Calculate_face (this, v_quantity, angle, edge_num)
        class (slip_vertex_t), intent (in out)     :: this      

        class(quantity_t), intent (in out) :: v_quantity  
        real (8)                           , intent (in)     :: angle    
        integer                            , intent (in)     :: edge_num 

    end subroutine

    subroutine Calculate_corner (this, v_quantity, angle, edge_num)
        implicit none
        class (slip_vertex_t), intent (in out)     :: this      

        class(quantity_t), intent (in out) :: v_quantity  
        real (8)                           , intent (in)     :: angle    
        integer                            , intent (in)     :: edge_num 

    end subroutine
end module slip_vertex_module
