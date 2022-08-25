
module slip_cell_3d_module
    use cell_boundary_condition_module, only : cell_boundary_condition_t
    use data_module    , only : data_t
    use quantity_module, only : quantity_t

    implicit none
    private

    type, extends(cell_boundary_condition_t), public :: slip_cell_3d_t
        private
    contains

        procedure, public :: Calculate => Slip_cell_3d_calculate
    end type slip_cell_3d_t

contains
    subroutine Slip_cell_3d_calculate (this, c_quantity, edge_num)
        class (slip_cell_3d_t) , intent (in out)     :: this      
        class(quantity_t)   , intent (in out)     :: c_quantity  
        integer             , intent (in)         :: edge_num  

        real(8), dimension(:,:,:), pointer :: values
        real(8), dimension(:,:,:,:), pointer :: values_4d
        integer :: i, j, k, nxp, nyp, nzp,m, nmats

        nxp = c_quantity%d1 + 1
        nyp = c_quantity%d2 + 1
        nzp = c_quantity%d3 + 1
        if (associated(c_quantity%data_4d)) then
            call c_quantity%Point_to_data(values_4d)
            nmats = c_quantity%data_4d%nmats
            select case(edge_num)
                case(1)
                    i = 1
                    do k = 0, nzp
                        do j = 0, nyp
                            do m = 1, nmats
                                values_4d(m, i-1, j, k) = values_4d(m, i, j, k)
                            end do
                        end do
                    end do


                case(2)
                    i = nxp
                    do k = 0, nzp
                        do j = 0, nyp
                            do m = 1, nmats

                                values_4d(m, i, j, k) = values_4d(m, i-1, j, k)
                            end do
                        end do
                    end do

                case(3)
                    j = 1
                    do k = 0, nzp
                        do i = 0, nxp
                            do m = 1, nmats

                                values_4d(m, i, j-1, k) = values_4d(m, i, j, k)
                            end do
                        end do
                    end do

                case(4)
                    j = nyp
                    do k = 0, nzp
                        do i = 0, nxp
                            do m = 1, nmats

                                values_4d(m, i, j, k) = values_4d(m, i, j-1, k)
                            end do
                        end do
                    end do


                case(5)
                    k = 1
                    do j = 0, nyp
                        do i = 0, nxp
                            do m = 1, nmats

                                values_4d(m, i, j, k-1) = values_4d(m, i, j, k)
                            end do
                        end do
                    end do

                case(6)
                    k = nzp
                    do j = 0, nyp
                        do i = 0, nxp
                            do m = 1, nmats

                                values_4d(m, i, j, k) = values_4d(m, i, j, k-1)
                            end do
                        end do
                    end do

            end select
        else
            call c_quantity%Point_to_data(values)

            select case(edge_num)
                case(1)
                    i = 1
                    do k = 0, nzp
                        do j = 0, nyp
                            values(i-1, j, k) = values(i, j, k)
                        end do
                    end do


                case(2)
                    i = nxp
                    do k = 0, nzp
                        do j = 0, nyp
                            values(i, j, k) = values(i-1, j, k)
                        end do
                    end do

                case(3)
                    j = 1
                    do k = 0, nzp
                        do i = 0, nxp
                            values(i, j-1, k) = values(i, j, k)
                        end do
                    end do

                case(4)
                    j = nyp
                    do k = 0, nzp
                        do i = 0, nxp
                            values(i, j, k) = values(i, j-1, k)
                        end do
                    end do

                case(5)
                    k = 1
                    do j = 0, nyp
                        do i = 0, nxp
                            values(i, j, k-1) = values(i, j, k)
                        end do
                    end do

                case(6)
                    k = nzp
                    do j = 0, nyp
                        do i = 0, nxp
                            values(i, j, k) = values(i, j, k-1)
                        end do
                    end do
            end select
        end if
        return
    end subroutine

end module slip_cell_3d_module
