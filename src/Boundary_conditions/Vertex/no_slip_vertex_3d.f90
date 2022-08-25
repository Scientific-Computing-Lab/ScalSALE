
module no_slip_vertex_3d_module
    use vertex_boundary_condition_module, only : vertex_boundary_condition_t
    use data_module, only : data_t
    use quantity_module, only : quantity_t
    implicit none
    private

    type, extends(vertex_boundary_condition_t), public :: no_slip_vertex_3d_t
        private
    contains

        procedure, public :: Calculate => No_slip_vertex_3d_calculate
        procedure, public :: Calculate_face
        procedure, public :: Calculate_corner
    end type no_slip_vertex_3d_t



contains
    subroutine No_slip_vertex_3d_calculate (this, v_quantity, coordinates, angle, edge_num)
        class (no_slip_vertex_3d_t)           , intent (in out) :: this      
        class(quantity_t), intent (in out) :: v_quantity  
        type(data_t), dimension(:), pointer,intent(inout)        :: coordinates
        real (8)                           , intent (in)     :: angle    
        integer                            , intent (in)     :: edge_num 

        real(8) :: normal_vel      

        real(8), dimension (:,:,:), pointer :: x 
        real(8), dimension (:,:,:), pointer :: y 

        real(8) :: slope        
        real(8) :: delta_x      
        real(8) :: delta_y      
        real(8) :: slope_factor 
        integer :: nxp, nyp, nzp
        integer :: virt_nxp, virt_nyp, virt_nzp
        integer :: vec_e_num, ivirt,jvirt,kvirt

        logical :: is_wall_x_top, is_wall_x_bot,is_wall_y_top, is_wall_y_bot,is_wall_z_top,is_wall_z_bot


        real(8), dimension (:,:,:), pointer :: values_dim1

        real(8), dimension (:,:,:), pointer :: values_dim2

        real(8), dimension (:,:,:), pointer :: values_dim3
        real(8), dimension (:), pointer :: edge_vector_x
        real(8), dimension (:), pointer :: edge_vector_y
        real(8), dimension (:), pointer :: edge_vector_z
        integer :: i, j, k
        call v_quantity%Point_to_data(values_dim1, values_dim2, values_dim3)
        call v_quantity%boundary_params%Point_to_edge_vectors(edge_vector_x, edge_vector_y, edge_vector_z)

        nxp = v_quantity%d1
        nyp = v_quantity%d2
        nzp = v_quantity%d3

        virt_nxp = v_quantity%boundary_params%parallel_params%virt_nxp
        virt_nyp = v_quantity%boundary_params%parallel_params%virt_nyp
        virt_nzp = v_quantity%boundary_params%parallel_params%virt_nzp

        select case (edge_num)
            case(1) 
                i = 1
                do k = 1, nzp
                    do j = 1, nyp
                        values_dim1(i, j, k) = 0d0
                        values_dim2(i, j, k) = 0d0
                        values_dim3(i, j, k) = 0d0
                    end do
                end do

            case(2) 
                i = nxp
                do k = 1, nzp
                    do j = 1, nyp
                        values_dim1(i, j, k) = 0d0
                        values_dim2(i, j, k) = 0d0
                        values_dim3(i, j, k) = 0d0
                    end do
                end do
            case(3) 
                j = 1
                do k = 1, nzp
                    do i = 1, nxp
                        values_dim1(i, j, k) = 0d0
                        values_dim2(i, j, k) = 0d0
                        values_dim3(i, j, k) = 0d0
                    end do
                end do
            case(4) 
                j = nyp
                do k = 1, nzp
                    do i = 1, nxp
                        values_dim1(i, j, k) = 0d0
                        values_dim2(i, j, k) = 0d0
                        values_dim3(i, j, k) = 0d0
                    end do
                end do
            case(5) 
                k = 1
                do j = 1, nyp
                    do i = 1, nxp
                        values_dim1(i, j, k) = 0d0
                        values_dim2(i, j, k) = 0d0
                        values_dim3(i, j, k) = 0d0
                    end do
                end do
            case(6) 
                k = nzp
                do j = 1, nyp
                    do i = 1, nxp
                        values_dim1(i, j, k) = 0d0
                        values_dim2(i, j, k) = 0d0
                        values_dim3(i, j, k) = 0d0
                    end do
                end do
        end select
    end subroutine No_slip_vertex_3d_calculate

    subroutine Calculate_face (this, v_quantity, angle, edge_num)
        class (no_slip_vertex_3d_t), intent (in out)     :: this      

        class(quantity_t), intent (in out) :: v_quantity  
        real (8)                           , intent (in)     :: angle    
        integer                            , intent (in)     :: edge_num 

    end subroutine

    subroutine Calculate_corner (this, v_quantity, angle, edge_num)
        implicit none
        class (no_slip_vertex_3d_t), intent (in out)     :: this      

        class(quantity_t), intent (in out) :: v_quantity  
        real (8)                           , intent (in)     :: angle    
        integer                            , intent (in)     :: edge_num 

    end subroutine
end module no_slip_vertex_3d_module
