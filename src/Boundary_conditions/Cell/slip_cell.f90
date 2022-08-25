
module slip_cell_module
    use cell_boundary_condition_module, only : cell_boundary_condition_t
    use data_module    , only : data_t
    use data_4d_module    , only : data_4d_t
    use quantity_module, only : quantity_t
    implicit none
    private

    type, extends(cell_boundary_condition_t), public :: slip_cell_t
        private
    contains

        procedure, public :: Calculate => Slip_cell_calculate
    end type Slip_cell_t

contains
    subroutine Slip_cell_calculate (this, c_quantity, edge_num)
        class (slip_cell_t) , intent (in out)     :: this      
        class(quantity_t)   , intent (in out)     :: c_quantity  
        integer             , intent (in)         :: edge_num  

        real(8), dimension(:,:,:), pointer :: values
        real(8), dimension(:, :,:,:), pointer :: values_4d
        integer :: i, j, nxp, nyp, nmats, m
        nxp = c_quantity%d1 + 1
        nyp = c_quantity%d2 + 1

        if (associated(c_quantity%data_4d)) then
            nmats = c_quantity%data_4d%nmats
            call c_quantity%Point_to_data(values_4d)
            select case(edge_num)
                case(1)
                    i = 1
                    do j = 0, nyp
                        do m = 1, nmats
                            values_4d(m, i - 1, j, 1) = values_4d(m, i, j, 1)
                        end do
                    end do

                case(2)
                    i = nxp
                    do j = 0, nyp
                        do m = 1, nmats
                            values_4d(m, i, j, 1) = values_4d(m, i - 1, j, 1)
                        end do
                    end do

                case(3)
                    j = 1
                    do i = 0, nxp
                        do m = 1, nmats
                            values_4d(m, i, j - 1, 1) = values_4d(m, i, j, 1)
                        end do
                    end do

                case(4)
                    j = nyp
                    do i = 0, nxp
                        do m = 1, nmats
                            values_4d(m, i, j, 1) = values_4d(m, i, j - 1, 1)
                        end do
                    end do

            end select
        else
            call c_quantity%Point_to_data(values)
            select case(edge_num)
                case(1)
                    i = 1
                    do j = 0, nyp
                        values(i - 1, j, 1) = values(i, j, 1)
                    end do

                case(2)
                    i = nxp
                    do j = 0, nyp
                        values(i, j, 1) = values(i - 1, j, 1)
                    end do

                case(3)
                    j = 1
                    do i = 0, nxp
                        values(i, j - 1, 1) = values(i, j, 1)
                    end do

                case(4)
                    j = nyp
                    do i = 0, nxp
                        values(i, j, 1) = values(i, j - 1, 1)
                    end do
            end select
        end if






    end subroutine

end module slip_cell_module
