
module slip_vertex_3d_module
    use vertex_boundary_condition_module, only : vertex_boundary_condition_t
    use data_module, only : data_t
    use quantity_module, only : quantity_t
    implicit none
    private

    type, extends(vertex_boundary_condition_t), public :: slip_vertex_3d_t
        private
    contains

        procedure, public  :: Calculate => Slip_vertex_calculate_3d


        procedure, private  :: Slip_vertex_calculate_3d
        procedure, private :: Calculate_on_edge
        procedure, private :: Calculate_egde

        procedure, public :: Calculate_face
        procedure, public :: Calculate_corner

    end type

contains


    subroutine Slip_vertex_calculate_3d (this, v_quantity, coordinates, angle, edge_num)
        implicit none
        class (slip_vertex_3d_t)              , intent (in out) :: this      

        class(quantity_t), intent (in out)                   :: v_quantity  
        type(data_t), dimension(:), pointer,intent(inout)    :: coordinates
        real (8)                           , intent (in)     :: angle    
        integer                            , intent (in)     :: edge_num 

        real(8) :: normal_vel      

        integer :: nxp, nyp, nzp
        integer :: virt_nxp, virt_nyp, virt_nzp
        integer :: vec_e_num
        integer, dimension(:),allocatable :: ivirt,jvirt,kvirt

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
        ivirt = v_quantity%boundary_params%parallel_params%i_virt
        jvirt = v_quantity%boundary_params%parallel_params%j_virt
        kvirt = v_quantity%boundary_params%parallel_params%k_virt

        select case (edge_num)
            case(1) 
                i = 1
                do k = 1, nzp
                    do j = 1, nyp
                        normal_vel = values_dim1(i, j, k) * v_quantity%boundary_params%normal_i1(1)&
                            + values_dim2(i, j, k) * v_quantity%boundary_params%normal_i1(2)&
                            + values_dim3(i, j, k) * v_quantity%boundary_params%normal_i1(3)
                        values_dim1(i, j, k) = values_dim1(i, j, k) - normal_vel*v_quantity%boundary_params%normal_i1(1)
                        values_dim2(i, j, k) = values_dim2(i, j, k) - normal_vel*v_quantity%boundary_params%normal_i1(2)
                        values_dim3(i, j, k) = values_dim3(i, j, k) - normal_vel*v_quantity%boundary_params%normal_i1(3)
                    end do
                end do

            case(2) 
                i = nxp
                do k = 1, nzp
                    do j = 1, nyp
                        normal_vel = values_dim1(i, j, k)*v_quantity%boundary_params%normal_inxp(1)&
                            + values_dim2(i, j, k)*v_quantity%boundary_params%normal_inxp(2)&
                            + values_dim3(i, j, k)*v_quantity%boundary_params%normal_inxp(3)
                        values_dim1(i, j, k) = values_dim1(i, j, k) - normal_vel*v_quantity%boundary_params%normal_inxp(1)
                        values_dim2(i, j, k) = values_dim2(i, j, k) - normal_vel*v_quantity%boundary_params%normal_inxp(2)
                        values_dim3(i, j, k) = values_dim3(i, j, k) - normal_vel*v_quantity%boundary_params%normal_inxp(3)
                    end do
                end do

            case(3) 
                j = 1
                do k = 1, nzp
                    do i = 1, nxp
                        normal_vel = values_dim1(i, j, k)*v_quantity%boundary_params%normal_j1(1) &
                            + values_dim2(i, j, k)*v_quantity%boundary_params%normal_j1(2) &
                            + values_dim3(i, j, k)*v_quantity%boundary_params%normal_j1(3)
                        values_dim1(i, j, k) = values_dim1(i, j, k) - normal_vel*v_quantity%boundary_params%normal_j1(1)
                        values_dim2(i, j, k) = values_dim2(i, j, k) - normal_vel*v_quantity%boundary_params%normal_j1(2)
                        values_dim3(i, j, k) = values_dim3(i, j, k) - normal_vel*v_quantity%boundary_params%normal_j1(3)
                    end do
                end do

            case(4) 
                j = nyp
                do k = 1, nzp
                    do i = 1, nxp
                        normal_vel = values_dim1(i, j, k) * v_quantity%boundary_params%normal_jnyp(1)&
                            + values_dim2(i, j, k) * v_quantity%boundary_params%normal_jnyp(2)&
                            + values_dim3(i, j, k) * v_quantity%boundary_params%normal_jnyp(3)
                        values_dim1(i, j, k) = values_dim1(i, j, k) - normal_vel * v_quantity%boundary_params%normal_jnyp(1)
                        values_dim2(i, j, k) = values_dim2(i, j, k) - normal_vel * v_quantity%boundary_params%normal_jnyp(2)
                        values_dim3(i, j, k) = values_dim3(i, j, k) - normal_vel * v_quantity%boundary_params%normal_jnyp(3)
                    end do
                end do

            case(5) 
                k = 1
                do j = 1, nyp
                    do i = 1, nxp
                        normal_vel = values_dim1(i, j, k)*v_quantity%boundary_params%normal_k1(1)&
                            + values_dim2(i, j, k)*v_quantity%boundary_params%normal_k1(2)&
                            + values_dim3(i, j, k)*v_quantity%boundary_params%normal_k1(3)
                        values_dim1(i, j, k) = values_dim1(i, j, k) - normal_vel*v_quantity%boundary_params%normal_k1(1)
                        values_dim2(i, j, k) = values_dim2(i, j, k) - normal_vel*v_quantity%boundary_params%normal_k1(2)
                        values_dim3(i, j, k) = values_dim3(i, j, k) - normal_vel*v_quantity%boundary_params%normal_k1(3)
                    end do
                end do
            case(6) 
                k = nzp
                do j = 1, nyp
                    do i = 1, nxp
                        normal_vel = values_dim1(i, j, k)*v_quantity%boundary_params%normal_knzp(1)&
                            + values_dim2(i, j, k)*v_quantity%boundary_params%normal_knzp(2)&
                            + values_dim3(i, j, k)*v_quantity%boundary_params%normal_knzp(3)
                        values_dim1(i, j, k) = values_dim1(i, j, k) - normal_vel*v_quantity%boundary_params%normal_knzp(1)
                        values_dim2(i, j, k) = values_dim2(i, j, k) - normal_vel*v_quantity%boundary_params%normal_knzp(2)
                        values_dim3(i, j, k) = values_dim3(i, j, k) - normal_vel*v_quantity%boundary_params%normal_knzp(3)
                    end do
                end do
        end select




    end subroutine Slip_vertex_calculate_3d

    subroutine Calculate_on_edge (this, values_dim1, values_dim2, values_dim3, edge_vector_x,edge_vector_y,edge_vector_z)
        implicit none
        class (slip_vertex_3d_t)              , intent (in out) :: this      

        real(8), intent (in out) :: values_dim1

        real(8), intent (in out) :: values_dim2

        real(8), intent (in out) :: values_dim3
        real(8), intent (in out) :: edge_vector_x
        real(8), intent (in out) :: edge_vector_y
        real(8), intent (in out) :: edge_vector_z
        real(8) :: scalar_mult         
        scalar_mult = values_dim1 * edge_vector_x + values_dim2 * edge_vector_y +&
            values_dim3 * edge_vector_z
        values_dim1 = scalar_mult*edge_vector_x
        values_dim2 = scalar_mult*edge_vector_y
        values_dim3 = scalar_mult*edge_vector_z

    end subroutine Calculate_on_edge

    subroutine Calculate_egde(this,values_dim1, values_dim2, values_dim3, edge_vector_x,edge_vector_y,edge_vector_z,&
        ivirt, jvirt, kvirt, virt_nxp, virt_nyp, virt_nzp  )
        class (slip_vertex_3d_t)              , intent (in out) :: this      

        real(8), intent (in out) :: values_dim1

        real(8), intent (in out) :: values_dim2

        real(8), intent (in out) :: values_dim3
        real(8), dimension(:), intent (in out) :: edge_vector_x
        real(8), dimension(:), intent (in out) :: edge_vector_y
        real(8), dimension(:), intent (in out) :: edge_vector_z
        integer, intent (in)     :: ivirt,jvirt,kvirt,virt_nxp, virt_nyp, virt_nzp

        integer :: vec_e_num
        vec_e_num = -1

        if (jvirt == 1        .and. kvirt == 1   ) then
            vec_e_num = 1
            call this%Calculate_on_edge(values_dim1, values_dim2, values_dim3, &
                edge_vector_x(vec_e_num), edge_vector_y(vec_e_num), edge_vector_z(vec_e_num))
        end if
        if (ivirt == virt_nxp .and. kvirt == 1   ) then
            vec_e_num = 2
            call this%Calculate_on_edge(values_dim1, values_dim2, values_dim3, &
                edge_vector_x(vec_e_num), edge_vector_y(vec_e_num), edge_vector_z(vec_e_num))
        end if
        if (jvirt == virt_nyp .and. kvirt == 1   ) then
            vec_e_num = 3
            call this%Calculate_on_edge(values_dim1, values_dim2, values_dim3, &
                edge_vector_x(vec_e_num), edge_vector_y(vec_e_num), edge_vector_z(vec_e_num))
        end if
        if (ivirt == 1        .and. kvirt == 1       ) then
            vec_e_num = 4
            call this%Calculate_on_edge(values_dim1, values_dim2, values_dim3, &
                edge_vector_x(vec_e_num), edge_vector_y(vec_e_num), edge_vector_z(vec_e_num))
        end if
        if (ivirt == 1        .and. jvirt == 1       ) then
            vec_e_num = 5
            call this%Calculate_on_edge(values_dim1, values_dim2, values_dim3, &
                edge_vector_x(vec_e_num), edge_vector_y(vec_e_num), edge_vector_z(vec_e_num))
        end if
        if (ivirt == virt_nxp .and. jvirt == 1       ) then
            vec_e_num = 6
            call this%Calculate_on_edge(values_dim1, values_dim2, values_dim3, &
                edge_vector_x(vec_e_num), edge_vector_y(vec_e_num), edge_vector_z(vec_e_num))
        end if
        if (ivirt == virt_nxp .and. jvirt == virt_nyp) then
            vec_e_num = 7
            call this%Calculate_on_edge(values_dim1, values_dim2, values_dim3, &
                edge_vector_x(vec_e_num), edge_vector_y(vec_e_num), edge_vector_z(vec_e_num))
        end if
        if (ivirt == 1        .and. jvirt == virt_nyp) then
            vec_e_num = 8
            call this%Calculate_on_edge(values_dim1, values_dim2, values_dim3, &
                edge_vector_x(vec_e_num), edge_vector_y(vec_e_num), edge_vector_z(vec_e_num))
        end if
        if (jvirt == 1        .and. kvirt == virt_nzp) then
            vec_e_num = 9
            call this%Calculate_on_edge(values_dim1, values_dim2, values_dim3, &
                edge_vector_x(vec_e_num), edge_vector_y(vec_e_num), edge_vector_z(vec_e_num))
        end if
        if (ivirt == virt_nxp .and. kvirt == virt_nzp    ) then
            vec_e_num = 10
            call this%Calculate_on_edge(values_dim1, values_dim2, values_dim3, &
                edge_vector_x(vec_e_num), edge_vector_y(vec_e_num), edge_vector_z(vec_e_num))
        end if
        if (jvirt == virt_nyp .and. kvirt == virt_nzp) then
            vec_e_num = 11
            call this%Calculate_on_edge(values_dim1, values_dim2, values_dim3, &
                edge_vector_x(vec_e_num), edge_vector_y(vec_e_num), edge_vector_z(vec_e_num))
        end if
        if (ivirt == 1        .and. kvirt == virt_nzp) then
            vec_e_num = 12
            call this%Calculate_on_edge(values_dim1, values_dim2, values_dim3, &
                edge_vector_x(vec_e_num), edge_vector_y(vec_e_num), edge_vector_z(vec_e_num))
        end if

    end subroutine Calculate_egde

    subroutine Calculate_face (this, v_quantity, angle, edge_num)
        class (slip_vertex_3d_t), intent (in out)     :: this      

        class(quantity_t), intent (in out) :: v_quantity  
        real (8)                           , intent (in)     :: angle    
        integer                            , intent (in)     :: edge_num 
        integer :: nxp, nyp, nzp
        integer :: virt_nxp, virt_nyp, virt_nzp
        integer :: vec_e_num
        integer, dimension(:),allocatable :: ivirt,jvirt,kvirt

        real(8), dimension (:,:,:), pointer :: values_dim1

        real(8), dimension (:,:,:), pointer :: values_dim2

        real(8), dimension (:,:,:), pointer :: values_dim3
        real(8), dimension (:), pointer :: edge_vector_x
        real(8), dimension (:), pointer :: edge_vector_y
        real(8), dimension (:), pointer :: edge_vector_z
        integer, dimension (:), allocatable :: boundary_type
        integer :: i, j, k

        call v_quantity%Point_to_data(values_dim1, values_dim2, values_dim3)
        call v_quantity%boundary_params%Point_to_edge_vectors(edge_vector_x, edge_vector_y, edge_vector_z)

        nxp = v_quantity%d1
        nyp = v_quantity%d2
        nzp = v_quantity%d3

        virt_nxp = v_quantity%boundary_params%parallel_params%virt_nxp
        virt_nyp = v_quantity%boundary_params%parallel_params%virt_nyp
        virt_nzp = v_quantity%boundary_params%parallel_params%virt_nzp
        ivirt = v_quantity%boundary_params%parallel_params%i_virt
        jvirt = v_quantity%boundary_params%parallel_params%j_virt
        kvirt = v_quantity%boundary_params%parallel_params%k_virt

        if (v_quantity%boundary_params%dimension == 2) return

        boundary_type = v_quantity%boundary_params%boundary_type
        select case (edge_num)
            case(1) 
                i = 1

                if (boundary_type(5) == 2) then
                    k = 1
                    i = 1
                    vec_e_num = 4
                    do j = 1, nyp
                        call this%Calculate_on_edge(values_dim1(i,j,k), values_dim2(i,j,k), values_dim3(i,j,k), &
                                edge_vector_x(vec_e_num), edge_vector_y(vec_e_num), edge_vector_z(vec_e_num))
                    end do
                end if

                if (boundary_type(6) == 2) then
                    i = 1
                    k = nzp
                    vec_e_num = 12
                    do j = 1, nyp
                        call this%Calculate_on_edge(values_dim1(i,j,k), values_dim2(i,j,k), values_dim3(i,j,k), &
                                edge_vector_x(vec_e_num), edge_vector_y(vec_e_num), edge_vector_z(vec_e_num))
                    end do
                end if

            case(2) 

                if (boundary_type(5) == 2) then
                    k = 1
                    i = nxp
                    vec_e_num = 2
                    do j = 1, nyp
                        call this%Calculate_on_edge(values_dim1(i,j,k), values_dim2(i,j,k), values_dim3(i,j,k), &
                                edge_vector_x(vec_e_num), edge_vector_y(vec_e_num), edge_vector_z(vec_e_num))
                    end do
                end if
                if (boundary_type(6) == 2) then
                    i = nxp
                    k = nzp
                    vec_e_num = 10
                    do j = 1, nyp
                        call this%Calculate_on_edge(values_dim1(i,j,k), values_dim2(i,j,k), values_dim3(i,j,k), &
                                edge_vector_x(vec_e_num), edge_vector_y(vec_e_num), edge_vector_z(vec_e_num))
                    end do
                end if
            case(3) 
                if (boundary_type(1) == 2) then
                    i = 1
                    j = 1
                    vec_e_num = 5
                    do k = 1, nzp
                        call this%Calculate_on_edge(values_dim1(i,j,k), values_dim2(i,j,k), values_dim3(i,j,k), &
                                edge_vector_x(vec_e_num), edge_vector_y(vec_e_num), edge_vector_z(vec_e_num))
                    end do
                end if

                if (boundary_type(2) == 2) then
                    i = nxp
                    j = 1
                    vec_e_num = 6
                    do k = 1, nzp
                        call this%Calculate_on_edge(values_dim1(i,j,k), values_dim2(i,j,k), values_dim3(i,j,k), &
                                edge_vector_x(vec_e_num), edge_vector_y(vec_e_num), edge_vector_z(vec_e_num))
                    end do
                end if
            case(4) 
                if (boundary_type(2) == 2) then
                    j = nyp
                    i = nxp
                    vec_e_num = 7
                    do k = 1, nzp
                        call this%Calculate_on_edge(values_dim1(i,j,k), values_dim2(i,j,k), values_dim3(i,j,k), &
                                edge_vector_x(vec_e_num), edge_vector_y(vec_e_num), edge_vector_z(vec_e_num))
                    end do
                end if

                if (boundary_type(1) == 2) then
                    j = nyp
                    i = 1
                    vec_e_num = 8
                    do k = 1, nzp
                        call this%Calculate_on_edge(values_dim1(i,j,k), values_dim2(i,j,k), values_dim3(i,j,k), &
                                edge_vector_x(vec_e_num), edge_vector_y(vec_e_num), edge_vector_z(vec_e_num))
                    end do
                end if

            case(5) 
                if (boundary_type(3) == 2) then
                    k = 1
                    j = 1
                    vec_e_num = 1
                    do i = 1, nxp
                        call this%Calculate_on_edge(values_dim1(i,j,k), values_dim2(i,j,k), values_dim3(i,j,k), &
                            edge_vector_x(vec_e_num), edge_vector_y(vec_e_num), edge_vector_z(vec_e_num))
                    end do
                end if

                if (boundary_type(4) == 2) then
                    k = 1
                    j = nyp
                    vec_e_num = 3
                    do i = 1, nxp
                        call this%Calculate_on_edge(values_dim1(i,j,k), values_dim2(i,j,k), values_dim3(i,j,k), &
                            edge_vector_x(vec_e_num), edge_vector_y(vec_e_num), edge_vector_z(vec_e_num))
                    end do
                end if


            case(6) 
                if (boundary_type(3) == 2) then
                    k = nzp
                    j = 1
                    vec_e_num = 9
                    do i = 1, nxp
                        call this%Calculate_on_edge(values_dim1(i,j,k), values_dim2(i,j,k), values_dim3(i,j,k), &
                            edge_vector_x(vec_e_num), edge_vector_y(vec_e_num), edge_vector_z(vec_e_num))
                    end do
                end if

                if (boundary_type(4) == 2) then
                    j = nyp
                    k = nzp
                    vec_e_num = 11
                    do i = 1, nxp
                        call this%Calculate_on_edge(values_dim1(i,j,k), values_dim2(i,j,k), values_dim3(i,j,k), &
                            edge_vector_x(vec_e_num), edge_vector_y(vec_e_num), edge_vector_z(vec_e_num))
                    end do
                end if
        end select
    end subroutine Calculate_face

    subroutine Calculate_corner (this, v_quantity, angle, edge_num)
        use geometry_module , only : Check_corner_vertex
        implicit none
        class (slip_vertex_3d_t), intent (in out)     :: this      

        class(quantity_t), intent (in out) :: v_quantity  
        real (8)                           , intent (in)     :: angle    
        integer                            , intent (in)     :: edge_num 
        integer :: nxp, nyp, nzp
        integer :: virt_nxp, virt_nyp, virt_nzp
        integer, dimension(:),allocatable :: ivirt,jvirt,kvirt

        real(8), dimension (:,:,:), pointer :: values_dim1

        real(8), dimension (:,:,:), pointer :: values_dim2

        real(8), dimension (:,:,:), pointer :: values_dim3
        integer :: i, j, k

        call v_quantity%Point_to_data(values_dim1, values_dim2, values_dim3)

        nxp = v_quantity%d1
        nyp = v_quantity%d2
        nzp = v_quantity%d3

        virt_nxp = v_quantity%boundary_params%parallel_params%virt_nxp
        virt_nyp = v_quantity%boundary_params%parallel_params%virt_nyp
        virt_nzp = v_quantity%boundary_params%parallel_params%virt_nzp
        ivirt = v_quantity%boundary_params%parallel_params%i_virt
        jvirt = v_quantity%boundary_params%parallel_params%j_virt
        kvirt = v_quantity%boundary_params%parallel_params%k_virt

        if (v_quantity%boundary_params%dimension == 2) return
        if (v_quantity%boundary_params%mesh_type /= 2) return
        select case (edge_num)
            case(1) 
                i = 1
                do k = 1, nzp
                    do j = 1, nyp
                        if (Check_corner_vertex(ivirt(i), jvirt(j), kvirt(k), virt_nxp, virt_nyp, virt_nzp) .eqv. .true.) then
                            values_dim1(i, j, k) = 0d0
                            values_dim2(i, j, k) = 0d0
                            values_dim3(i, j, k) = 0d0
                        end if
                    end do
                end do

            case(2) 
                i = nxp
                do k = 1, nzp
                    do j = 1, nyp
                        if (Check_corner_vertex(ivirt(i), jvirt(j), kvirt(k), virt_nxp, virt_nyp, virt_nzp) .eqv. .true.) then
                            values_dim1(i, j, k) = 0d0
                            values_dim2(i, j, k) = 0d0
                            values_dim3(i, j, k) = 0d0
                        end if
                    end do
                end do

            case(3) 
                j = 1
                do k = 1, nzp
                    do i = 1, nxp
                        if (Check_corner_vertex(ivirt(i), jvirt(j), kvirt(k), virt_nxp, virt_nyp, virt_nzp) .eqv. .true.) then
                            values_dim1(i, j, k) = 0d0
                            values_dim2(i, j, k) = 0d0
                            values_dim3(i, j, k) = 0d0
                        end if
                    end do
                end do

            case(4) 
                j = nyp
                do k = 1, nzp
                    do i = 1, nxp
                        if (Check_corner_vertex(ivirt(i), jvirt(j), kvirt(k), virt_nxp, virt_nyp, virt_nzp) .eqv. .true.) then
                            values_dim1(i, j, k) = 0d0
                            values_dim2(i, j, k) = 0d0
                            values_dim3(i, j, k) = 0d0
                        end if
                    end do
                end do

            case(5) 
                k = 1
                do j = 1, nyp
                    do i = 1, nxp
                        if (Check_corner_vertex(ivirt(i), jvirt(j), kvirt(k), virt_nxp, virt_nyp, virt_nzp) .eqv. .true.) then
                            values_dim1(i, j, k) = 0d0
                            values_dim2(i, j, k) = 0d0
                            values_dim3(i, j, k) = 0d0
                        end if
                    end do
                end do
            case(6) 
                k = nzp
                do j = 1, nyp
                    do i = 1, nxp
                        if (Check_corner_vertex(ivirt(i), jvirt(j), kvirt(k), virt_nxp, virt_nyp, virt_nzp) .eqv. .true.) then
                            values_dim1(i, j, k) = 0d0
                            values_dim2(i, j, k) = 0d0
                            values_dim3(i, j, k) = 0d0
                        end if
                    end do
                end do
        end select
    end subroutine Calculate_corner


end module slip_vertex_3d_module
