

module mesh_3d_module
use data_module, only : data_t
use datafile_module, only : datafile_t
    use coordinates_module              , only : coordinates_t
    use parallel_parameters_module      , only : parallel_parameters_t
    use communication_module            , only : communication_t
    use communication_parameters_module , only : communication_parameters_t
    use boundary_parameters_module      , only : boundary_parameters_t
    use vertex_boundary_condition_module, only : vertex_boundary_condition_t, vertex_bc_wrapper_t
use mesh_base_module, only : mesh_base_t
implicit none
    private

    public :: mesh_3d_t
	type, extends(mesh_base_t) :: mesh_3d_t
        private

        type (data_t), dimension(:), pointer, public :: avg_coordinates    



        real(8), dimension(:, :), allocatable        :: theta               
        real(8), dimension(:, :), allocatable        :: phi                 
        integer, dimension(:, :), pointer            :: start_layer_index_p 
        integer, dimension(:)  , allocatable, public :: n_layers_p   

    contains

        procedure, public    :: Calculate_angles => Calculate_angles_3d

        procedure, public :: Set_communication_mesh_3d

        procedure, public :: Set_communication => Set_communication_mesh_3d

        procedure, public    :: Calculate_average_coordinates => Calculate_average_coordinates_3d

        procedure, public    :: Read_and_init_mesh_3d         
        procedure, private   :: Build_angles_pyramid_3d
        procedure, private   :: Build_mesh => Build_mesh_3d
        procedure, private   :: Build_angles => Build_angles_3d           
        procedure, private   :: Build_coordinates => Build_coordinates_3d

        procedure, public :: Write_mesh_3d
        procedure, public    :: Write_mesh_abstract => Write_mesh_3d

        procedure, public :: Read_mesh_3d
        procedure, public    :: Read_mesh_abstract => Read_mesh_3d

    end type mesh_3d_t

    interface mesh_3d_t

        module procedure Constructor_3d

    end interface mesh_3d_t

contains

    type(mesh_3d_t) function Constructor_3d(df,bc_v, bc_params, parallel_params)
        implicit none
        type(datafile_t), intent(in)                :: df
        type(parallel_parameters_t), pointer, intent(inout)  :: parallel_params
        type (vertex_bc_wrapper_t) , dimension(:), pointer, intent(in out) :: bc_v           
        type(boundary_parameters_t), pointer, intent(in) :: bc_params

        type(data_t)        :: xx, yy, zz, xx_avg, yy_avg, zz_avg
        integer :: nxp, nyp, nzp

        nxp = parallel_params%nxp
        nyp = parallel_params%nyp
        nzp = parallel_params%nzp

        call Constructor_3d%Init_mesh_base(df, parallel_params)
        allocate(Constructor_3d%coordinates)
        Constructor_3d%coordinates = coordinates_t(0d0, nxp + 1, nyp + 1, nzp + 1,bc_v, bc_params)



        allocate (data_t :: Constructor_3d%avg_coordinates (3))


        allocate (Constructor_3d%n_layers_p (1))
        Constructor_3d%n_layers_p        = df%number_layers_k


        allocate (Constructor_3d%theta (nxp + 1, nyp + 1))

        allocate (Constructor_3d%phi (nxp + 1, nyp + 1))

        allocate(Constructor_3d%start_layer_index_p,source= df%start_layer_index_p)
        Constructor_3d%avg_coordinates(1) = data_t(nxp + 1, nyp + 1, nzp + 1)
        Constructor_3d%avg_coordinates(2) = data_t(nxp + 1, nyp + 1, nzp + 1)
        Constructor_3d%avg_coordinates(3) = data_t (nxp + 1, nyp + 1, nzp + 1)
        Constructor_3d%dimension      = 3

        call Constructor_3d%Read_and_init_mesh_3d(df)
    end function


    subroutine Calculate_angles_3d(this)
        implicit none
        class (mesh_3d_t)                :: this 
    end subroutine Calculate_angles_3d

    subroutine Read_and_init_mesh_3d(this, df)
        use constants_module, only: PI
        implicit none
        type(datafile_t), intent(in)     :: df

        class (mesh_3d_t), intent(in out)                  :: this            

        integer              ::  layer_num                
        integer              ::  mesh_id                  
        integer              ::  i, j                     
        integer              ::  temp_cells_num           
        integer              ::  total_counter            
        integer              ::  contour_counter          
        integer              ::  number_of_contours       
        integer              ::  contour_segment          
        integer              ::  contour_segments         
        real(8), allocatable ::  contour_params(:,:,:,:)  

        real(8)              ::  y_rad, x_rad, y_power, x_power, y_c, x_c, x1, y1, x2, y2, teta1, &
            teta2, xh1, &
            & xh2, yh1, yh2      


        logical              ::  same_material_in_layer   
        logical              ::  axis_geometry            

        axis_geometry = .false.

        same_material_in_layer = .true.
        contour_segments = 1 

        if (this%mesh_type == 1) then 

        end if

        call this%Build_mesh (df, 1, 1) 

        if (this%mesh_type  == 2) then
            this%sphere_factor = 2d0 * PI
            if (axis_geometry .eqv. .false.) then 
                this%sphere_factor = 1.d0
            end if
        end if
        return
    end subroutine Read_and_init_mesh_3d


    subroutine Build_mesh_3d(this, df, sw_is, sw_js)
        implicit none

        type(datafile_t)  , intent(in)         :: df
        class (mesh_3d_t)    , intent(in out)  :: this
        integer           , intent(in)         :: sw_is, sw_js

        real(8), dimension(:, :), allocatable           :: rr
        integer                                         :: i,j,k
        select case(this%mesh_type)
            case (2) 
                call this%Build_angles(df%theta0, df%phi0,sw_js, sw_is, df%contour_j_type, df%contour_k_type)
            case (1)                  
                write(*,*)"Building pyramid"
                call this%Build_angles_pyramid_3d(df%theta0, df%phi0,sw_js, sw_is, df%contour_j_type, df%contour_k_type)
            case default
                write(*,*) "ERROR unrecognize mesh_type: ", this%mesh_type
                stop
        end select


        call this%Build_coordinates(df%contour_i_type, df%contour_i_line, df%zone_i_type, df%zone_i_d, sw_is)

    end subroutine Build_mesh_3d


    subroutine Build_coordinates_3d(this, contour_i_type, contour_line,zone_type,zone_dr, sw_is)
        use geometry_module, only : Layer_cut_line_3d,Layer_cut_elipse_3d, Calculate_constant,Calculate_geometry_series, Root
        implicit none
        class(mesh_3d_t)              , intent(inout)         :: this
        integer       , dimension(:), intent(in)              :: contour_i_type
        integer                       , intent(in)            :: sw_is
        real(8)       , dimension(:,:), intent(in)            :: contour_line
        real(8)       , dimension(:), intent(in)              :: zone_dr 
        integer       , dimension(:), intent(in)              :: zone_type

        integer                                               :: n_part_max
        real(8)                                               :: r0, y11, x11, y12, x12, rr1, rr2, rr3, xh1, xh2&
            , xh3, yh1, yh2, yh3, teta
        real(8)                                               :: x_layer_up(1, 11), x_layer_down(1, 11), x_layer_up1(1, 11) 
        integer                                               :: i, j, j1, j2, i1, i2, mm, ll, n, j11, j22, k1, k2
        integer                                               :: tmp_start, tmp_end
        integer                                               :: k_start, k_end,k
        integer                                               :: start, end
        integer                                               :: counter,tmp
        real(8)                                               :: sin_teta, dr, qq,ra1,ra2,x1,x2,y1,y2,z1,z2
        real(8), dimension(:), pointer   :: rr
        real(8)       , dimension(:,:),allocatable          :: contour_i_line

        allocate(contour_i_line, source=contour_line)
        call this%Fix_contour_i(contour_i_type, contour_i_line)


        call this%parallel_params%Get_virtual_from_to (tmp_start, tmp_end, 3)
        tmp_start = tmp_start + 1
        tmp_end = tmp_end - 1

        n_part_max = 1
        allocate(rr(0:this%start_layer_index_r(this%n_layers_r(1) + 1, 1)))
        i1 = 1         
        i2 = this%nxp       
        j1 = 1         
        j2 = this%nyp       

        do n = 1, this%n_layers_r(1)
            do mm = 1, 1
                do ll = 1, size(contour_i_line(n, :))
                    x_layer_down(mm, ll) = contour_i_line(n, ll) 
                    x_layer_up(mm, ll)   = contour_i_line(n + 1, ll) 
                    if (this%n_layers_r(1) > n) then
                        x_layer_up1(mm, ll) = contour_i_line(n + 2, ll)
                    end if
                    k1 = this%start_layer_index_r(n, 1) + (sw_is - 1)
                    k2 = this%start_layer_index_r(n + 1, 1) + (sw_is - 1)
                    k_start = k1
                    k_end = k2
                    if (k1 < tmp_start) then 
                        k_start = tmp_start
                    end if
                    if (k2 > tmp_end) then
                        k_end = tmp_end
                    end if
                    tmp = this%parallel_params%get_physical_index(k_start, 3)
                end do 
            end do  
            do j = j1, j2
                do i = i1, i2
                    teta = this%theta(i, j)
                    r0 = 1000.0d0
                    y11 = 0.0d0
                    x11 = 0.0d0
                    y12 = y11 + r0 * cos(teta)
                    x12 = r0 * sin(teta)
                    if (this%mesh_type == 2) then
                        y11 = teta
                        x11 = -r0
                        y12 = teta
                        x12 = r0
                    end if

                    select case(contour_i_type(n))
                        case (0) 
                            call Layer_cut_line_3d (y11, x11, y12, x12, x_layer_down, n_part_max, yh1, xh1)
                            call Layer_cut_line_3d (y11, x11, y12, x12, x_layer_up  , n_part_max, yh2, xh2)
                        case (1)                  
                            call Layer_cut_elipse_3d (y11, x11, y12, x12, x_layer_down, n_part_max, yh1, xh1)
                            call Layer_cut_elipse_3d (y11, x11, y12, x12, x_layer_up  , n_part_max, yh2, xh2)
                        case default
                            write(*,*) "ERROR no contour type number: ", contour_i_type(n)
                    end select

                    if (this%n_layers_r(1) > n) then
                        select case(contour_i_type(n))
                            case (0) 
                                call Layer_cut_line_3d (y11, x11, y12, x12, x_layer_up1, n_part_max, yh3, xh3)
                            case (1)                  
                                call Layer_cut_elipse_3d (y11, x11, y12, x12, x_layer_up1, n_part_max, yh3, xh3)
                            case default
                                write(*,*) "ERROR no contour type number: ", contour_i_type(n)
                        end select
                    end if

                    if (this%mesh_type == 2) then
                        rr1 = xh1
                        rr2 = xh2
                        if (this%n_layers_r(1) > n) rr3 = xh3
                    else
                        rr1 = sqrt(xh1**2 + yh1**2)
                        rr2 = sqrt(xh2**2 + yh2**2)
                        if (this%n_layers_r(1) > n) rr3 = sqrt(xh3**2 + yh3**2)
                    end if

                    select case(zone_type(n))
                        case (0)

                            call Calculate_constant(rr, rr1, rr2, k_start, k2, tmp)
                        case (2)
                            dr = zone_dr(n)
                            if (dr < 0d0) then
                                x1 = this%coordinates%data(1)%values(i, j, k1)
                                x2 = this%coordinates%data(1)%values(i, j, k1 - 1)
                                y1 = this%coordinates%data(2)%values(i, j, k1)
                                y2 = this%coordinates%data(2)%values(i, j, k1 - 1)
                                z1 = this%coordinates%data(3)%values(i, j, k1)
                                z2 = this%coordinates%data(3)%values(i, j, k1 - 1)
                                ra1 = sqrt(x1**2 + y1**2 + z1**2) 
                                ra2 = sqrt(x2**2 + y2**2 + z2**2) 
                                dr = ra1 - ra2
                                if (this%mesh_type == 2) then
                                    dr = abs(x1 - x2)
                                end if
                                dr = dr * abs(zone_dr(n))
                            end if
                            qq = root(1.05d0, dr, (rr2-rr1), k2 - k1)
                            rr(k_start) = rr1

                            call Calculate_geometry_series(rr, dr, qq, k_start + 1, k2, tmp)
                        case default
                            write(*,*) "ERROR no zone type number: ", zone_type(n)
                    end select

                    k = tmp
                    if (this%mesh_type == 2) then
                        do counter = k_start, k_end
                            this%coordinates%data(1)%values(i, j, k) = this%theta(i, j)
                            this%coordinates%data(2)%values(i, j, k) = this%phi(i, j)
                            this%coordinates%data(3)%values(i, j, k) = rr(counter)
!if(this%parallel_params%my_rank == 1) write(*,*) i,j,k,counter,k1
                            k = k + 1
                        end do
                    else
                        do counter = k_start, k_end
                            this%coordinates%data(1)%values(i, j, k) = rr(counter) * sin(this%theta(i,j)) * cos(this%phi(i, j))
                            this%coordinates%data(2)%values(i, j, k) = rr(counter) * sin(this%theta(i,j)) * sin(this%phi(i, j))
                            this%coordinates%data(3)%values(i, j, k) = rr(counter) * cos(this%theta(i,j))
                            k = k + 1
                        end do
                    end if
                end do
            end do
        end do
        deallocate(rr)
        deallocate(contour_i_line)
    end subroutine Build_coordinates_3d

    subroutine Build_angles_3d(this ,theta0_df, phi0_df,sw_is,sw_js, contour_j_type,contour_k_type)
        use constants_module, only : PI
        use geometry_module, only : Calculate_constant

        implicit none
        class(mesh_3d_t)                     , intent(inout)  :: this
        integer, dimension(:) ,allocatable   , intent(in)     :: contour_j_type
        integer, dimension(:) ,allocatable   , intent(in)     :: contour_k_type
        integer                              , intent(in)     :: sw_js,sw_is
        real(8), dimension(:,:), allocatable , intent(in)     :: theta0_df
        real(8), dimension(:,:), allocatable , intent(in)     :: phi0_df

        integer                                        :: i,i1,i2, j1, j2, n, j,nx ,ny
        real(8)                                        :: theta1, theta2 , theta_from, theta_to, theta_all
        real(8)                                        :: phi1, phi2 , phi_from, phi_to, phi_all
        real(8), dimension(:), pointer                 :: virt_theta, virt_phi
        real(8), dimension(:,:), allocatable           :: theta0
        real(8), dimension(:,:), allocatable           :: phi0
        integer :: i_start, i_end, j_start, j_end

        allocate(theta0, source=theta0_df)
        allocate(phi0  , source=phi0_df)

        call this%Fix_contour_angle(theta0, contour_j_type)
        call this%Fix_contour_angle(phi0  , contour_k_type)

        theta_from = theta0(1, 1)
        theta_to   = theta0(this%n_layers_t(1) + 1, 1)
        theta_all = theta_from - theta_to

        phi_from = phi0(1, 1)
        phi_to   = phi0(this%n_layers_p(1) + 1, 1)
        phi_all = phi_from - phi_to

        this%sphere_factor = 1d0
        nx = this%nxp - 1
        ny = this%nyp - 1


        allocate (virt_theta(0 : this%parallel_params%virt_nxp + 1))
        allocate (virt_phi(0 : this%parallel_params%virt_nyp + 1))

        virt_theta = 0d0
        virt_phi = 0d0


        do n = 1, this%n_layers_t(1)
            i1 = this%start_layer_index_t(n, 1) + (sw_is - 1)
            i2 = this%start_layer_index_t(n + 1, 1) + (sw_is - 1)
            theta1 = theta0(n, 1)
            theta2 = theta0(n + 1, 1)
            call Calculate_constant(virt_theta, theta1, theta2, i1, i2)


        end do 




        do n = 1, this%n_layers_p(1)
            j1 = this%start_layer_index_p(n, 1) + (sw_js - 1)
            j2 = this%start_layer_index_p(n + 1, 1) + (sw_js - 1)
            phi1 = phi0(n, 1)
            phi2 = phi0(n + 1, 1)

            call Calculate_constant(virt_phi, phi1, phi2, j1, j2)

        end do



        call this%parallel_params% Get_virtual_from_to(i_start, i_end, 1)
        call this%parallel_params% Get_virtual_from_to(j_start, j_end, 2)

        do j = 1, this%nyp
            do i = 1, this%nxp
                this%theta(i, j) = virt_theta(this%parallel_params%i_virt(i))
                this%phi(i, j)   = virt_phi(this%parallel_params%j_virt(j))
            end do
        end do

        deallocate(virt_theta) 
        deallocate(virt_phi)

        deallocate(theta0)
        deallocate(phi0)



    end subroutine Build_angles_3d

    subroutine Build_angles_pyramid_3d(this ,theta0_df, phi0_df,sw_is,sw_js, contour_j_type,contour_k_type)
        use constants_module, only : PI
        use geometry_module, only : Calculate_constant, Root, Triangle_method

        implicit none
        class(mesh_3d_t)                     , intent(inout)  :: this
        integer, dimension(:) ,allocatable   , intent(in)     :: contour_j_type
        integer, dimension(:) ,allocatable   , intent(in)     :: contour_k_type
        integer                              , intent(in)     :: sw_js,sw_is
        real(8), dimension(:,:), allocatable , intent(in)     :: theta0_df
        real(8), dimension(:,:), allocatable , intent(in)     :: phi0_df

        integer                                        :: i,i1,i2, j1, j2, n, j,nx ,ny
        real(8)                                        :: theta1, theta2 , theta_from, theta_to, theta_all
        real(8)                                        :: phi1, phi2 , phi_from, phi_to, phi_all
        real(8), dimension(:), pointer                 :: virt_theta, virt_phi
        real(8), dimension(:,:), allocatable           :: theta0
        real(8), dimension(:,:), allocatable           :: phi0
        integer :: i_start, i_end, j_start, j_end, virt_nxp,virt_nyp

        real(8) :: theta_top_1_1         
        real(8) :: theta_top_nxp_nyp     
        real(8) :: theta_top_1_nyp       
        real(8) :: theta_top_nxp_1       
        real(8) :: phi_top_1_1           
        real(8) :: phi_top_nxp_nyp       
        real(8) :: phi_top_1_nyp         
        real(8) :: phi_top_nxp_1         
        real(8) :: x_top_1_1             
        real(8) :: x_top_nxp_nyp         
        real(8) :: x_top_1_nyp           
        real(8) :: x_top_nxp_1           
        real(8) :: y_top_1_1             
        real(8) :: y_top_nxp_nyp         
        real(8) :: y_top_1_nyp           
        real(8) :: y_top_nxp_1           
        real(8) :: z_top_1_1             
        real(8) :: z_top_nxp_nyp         
        real(8) :: z_top_1_nyp           
        real(8) :: z_top_nxp_1           
        real(8) :: face_1_1              
        real(8) :: face_nxp_nyp          
        real(8) :: face_1_nyp            
        real(8) :: face_nxp_1            
        real(8) :: face_i,face_j
        real(8) :: x_top, y_top, z_top, rr, asca, bsca, csca, big_a, big_b, big_c, area_pyramid, tmp_len, &
            tmp_mat(3, 3), tmp_mat1(3,3), tmp_mat2(3,3), unitcell_area, face_i_1, face_i_nyp, face_1_j, &
            face_nxp_j, c_area_min, c_area_max, c_area_sum, tur_dtet1, tur_dtet2, qq1, qq2, fact1, fact2, &
            dtet_last, tet, car, car_min, car_max, car_sum, face_v, face_o, c_area
        real(8), dimension(:), allocatable :: face_i_1_virt
        real(8), dimension(:), allocatable :: face_i_nyp_virt
        real(8), dimension(:), allocatable :: face_1_j_virt
        real(8), dimension(:), allocatable :: face_nxp_j_virt

        real(8) :: theta_opening 

        integer :: k, jj, i_face_1,i_face_2, virt_nx, virt_ny, ivirt, jvirt
        allocate(theta0, source=theta0_df)
        allocate(phi0  , source=phi0_df)

        call this%Fix_contour_angle(theta0, contour_j_type)
        call this%Fix_contour_angle(phi0  , contour_k_type)

        theta_from = theta0(1, 1)
        theta_to   = theta0(this%n_layers_t(1) + 1, 1)
        theta_all = theta_from - theta_to

        phi_from = phi0(1, 1)
        phi_to   = phi0(this%n_layers_p(1) + 1, 1)
        phi_all = phi_from - phi_to

        this%sphere_factor = 1d0
        nx = this%nxp - 1
        ny = this%nyp - 1
        virt_nx  = this%parallel_params%virt_nx
        virt_ny  = this%parallel_params%virt_ny
        virt_nxp = this%parallel_params%virt_nxp
        virt_nyp = this%parallel_params%virt_nyp

        allocate (face_i_1_virt(this%parallel_params%virt_nxp))
        allocate (face_1_j_virt(this%parallel_params%virt_nyp))
        allocate (face_i_nyp_virt(this%parallel_params%virt_nxp))
        allocate (face_nxp_j_virt(this%parallel_params%virt_nyp))

        theta_opening     = 1d-3
        theta_top_1_1     = theta_opening
        theta_top_nxp_nyp = theta_opening
        theta_top_1_nyp   = theta_opening
        theta_top_nxp_1   = theta_opening
        phi_top_1_1       = -PI
        phi_top_nxp_nyp   = 0d0
        phi_top_1_nyp     = PI/2d0
        phi_top_nxp_1     = -PI/2d0

        x_top_1_1      = sin(theta_top_1_1) * cos(phi_top_1_1)
        y_top_1_1      = sin(theta_top_1_1) * sin(phi_top_1_1)
        z_top_1_1      = cos(theta_top_1_1)

        x_top_1_nyp    = sin(theta_top_1_nyp) * cos(phi_top_1_nyp)
        y_top_1_nyp    = sin(theta_top_1_nyp) * sin(phi_top_1_nyp)
        z_top_1_nyp    = cos(theta_top_1_nyp)

        x_top_nxp_1    = sin(theta_top_nxp_1) * cos(phi_top_nxp_1)
        y_top_nxp_1    = sin(theta_top_nxp_1) * sin(phi_top_nxp_1)
        z_top_nxp_1    = cos(theta_top_nxp_1)

        x_top_nxp_nyp  = sin(theta_top_nxp_nyp) * cos(phi_top_nxp_nyp)
        y_top_nxp_nyp  = sin(theta_top_nxp_nyp) * sin(phi_top_nxp_nyp)
        z_top_nxp_nyp  = cos(theta_top_nxp_nyp)

        unitcell_area = 4d0 * PI / 60d0 / dble(virt_nx * virt_ny)
        asca = acos(x_top_1_1    *x_top_1_nyp + y_top_1_1*y_top_1_nyp     + z_top_1_1*z_top_1_nyp)
        bsca = acos(x_top_nxp_nyp*x_top_1_nyp + y_top_nxp_nyp*y_top_1_nyp + z_top_nxp_nyp*z_top_1_nyp)
        csca = acos(x_top_nxp_nyp*x_top_1_1   + y_top_nxp_nyp*y_top_1_1   + z_top_nxp_nyp*z_top_1_1)

        big_a = acos( (cos(asca) - cos(bsca)*cos(csca)) / sin(bsca)/sin(csca) )
        big_b = acos( (cos(bsca) - cos(asca)*cos(csca)) / sin(asca)/sin(csca) )
        big_c = acos( (cos(csca) - cos(asca)*cos(bsca)) / sin(asca)/sin(bsca) )

        area_pyramid = big_a + big_b + big_c - PI

        asca = acos(x_top_1_1*x_top_nxp_1     + y_top_1_1*y_top_nxp_1     + z_top_1_1*z_top_nxp_1)
        bsca = acos(x_top_nxp_nyp*x_top_nxp_1 + y_top_nxp_nyp*y_top_nxp_1 + z_top_nxp_nyp*z_top_nxp_1)
        csca = acos(x_top_nxp_nyp*x_top_1_1   + y_top_nxp_nyp*y_top_1_1   + z_top_nxp_nyp*z_top_1_1)

        big_a = acos( (cos(asca) - cos(bsca)*cos(csca)) / sin(bsca)/sin(csca) )
        big_b = acos( (cos(bsca) - cos(asca)*cos(csca)) / sin(asca)/sin(csca) )
        big_c = acos( (cos(csca) - cos(asca)*cos(bsca)) / sin(asca)/sin(bsca) )

        area_pyramid = area_pyramid + big_a + big_b + big_c - PI

        this%sphere_factor = 4d0 * PI / area_pyramid


        fact1 = 1d0
        fact2 = 1d0
        tet = acos(x_top_1_1*x_top_nxp_1 + y_top_1_1*y_top_nxp_1 + z_top_1_1*z_top_nxp_1)
        tur_dtet1 = tet/dble(virt_nx) * fact1
        qq1 = Root(1.05d0, tur_dtet1, tet, virt_nx)

        face_i_1_virt(1) = 0d0
        dtet_last = tur_dtet1 / qq1
        do i = 2, virt_nxp
            face_i_1_virt(i) = face_i_1_virt(i - 1) + dtet_last * qq1
            dtet_last = dtet_last * qq1
        end do
        do i = 1, virt_nxp
            face_i_1_virt(i) = face_i_1_virt(i)/tet
        end do

        face_1_j_virt(1) = 0d0
        dtet_last = tur_dtet1 / qq1


        do i = 2, virt_nyp
            face_1_j_virt(i) = face_1_j_virt(i - 1) + dtet_last * qq1
            dtet_last = dtet_last * qq1
        end do
        do i = 1, virt_nyp
            face_1_j_virt(i) = face_1_j_virt(i)/tet
        end do

        tet = acos(x_top_nxp_nyp*x_top_nxp_1 + y_top_nxp_nyp*y_top_nxp_1 + z_top_nxp_nyp*z_top_nxp_1)
        tur_dtet2 = tet / dble(virt_ny) * fact2
        qq2 = Root(1.05d0, tur_dtet2, tet, virt_ny)
        face_i_nyp_virt(1) = 0d0
        dtet_last = tur_dtet2 / qq2
        do i = 2, virt_nxp
            face_i_nyp_virt(i) = face_i_nyp_virt(i - 1) + dtet_last * qq2
            dtet_last = dtet_last * qq2
        end do
        do i = 1, virt_nxp
            face_i_nyp_virt(i) = face_i_nyp_virt(i)/tet
        end do

        face_nxp_j_virt(1) = 0d0
        dtet_last = tur_dtet1 / qq2
        do i = 2, virt_nyp
            face_nxp_j_virt(i) = face_nxp_j_virt(i - 1) + dtet_last * qq2
            dtet_last = dtet_last * qq2
        end do
        do i = 1, virt_nyp
            face_nxp_j_virt(i) = face_nxp_j_virt(i)/tet
        end do

        do i = 1, this%nxp
            ivirt = this%parallel_params%i_virt(i)
            face_i = dble(ivirt - 1) / dble(virt_nxp - 1)
            do j = 1, this%nyp
                jvirt  = this%parallel_params%j_virt(j)
                face_j = dble(jvirt - 1) / dble(virt_nyp - 1)

                face_i_1    = face_i_1_virt  (ivirt)
                face_i_nyp  = face_i_nyp_virt(ivirt)
                face_1_j    = face_1_j_virt  (jvirt)
                face_nxp_j  = face_nxp_j_virt(jvirt)

                call Triangle_method(x_top_1_1    , y_top_1_1    , z_top_1_1,&
                                  x_top_1_nyp  , y_top_1_nyp  , z_top_1_nyp, &
                                  x_top_nxp_1  , y_top_nxp_1  , z_top_nxp_1, &
                                  x_top_nxp_nyp, y_top_nxp_nyp, z_top_nxp_nyp,&
                                  face_i_1, face_i_nyp, face_1_j, face_nxp_j,     &
                                  x_top,y_top, z_top)
                rr = sqrt(x_top*x_top + y_top*y_top + z_top*z_top)
                x_top = x_top / rr
                y_top = y_top / rr
                z_top = z_top / rr
                this%theta(i, j) = acos(z_top)
                this%phi  (i, j) = atan2(y_top, x_top)
            end do
        end do
        deallocate (face_i_1_virt)
        deallocate (face_1_j_virt)
        deallocate (face_i_nyp_virt)
        deallocate (face_nxp_j_virt)

    end subroutine Build_angles_pyramid_3d


    subroutine Calculate_average_coordinates_3d(this)
        implicit none
        class(mesh_3d_t), intent(inout)   :: this

        integer                           :: i, j, k
        real(8), dimension(:, :, :), pointer ::x, y, z
        call this%Point_to_data (x, y, z)
        do k = 1, this%nzp -1
            do j = 1, this%nyp - 1
                do i = 1, this%nxp - 1
                    this%avg_coordinates(1)%values(i, j, k) = 0.125d0 * ( x(i, j, k) + x (i + 1, j, k) + x(i + 1, j + 1, k)&
                        + x(i, j + 1, k) + x(i, j, k + 1) + x(i + 1, j, k + 1) + x(i + 1, j + 1, k + 1) + x(i, j + 1, k + 1))
                    this%avg_coordinates(2)%values(i, j, k) = 0.125d0 * ( y(i, j, k) + y(i + 1, j, k) + y(i + 1, j + 1, k)&
                        + y(i, j + 1, k) + y(i, j, k + 1) + y(i + 1, j, k + 1) + y(i + 1, j + 1, k + 1) + y(i, j + 1, k + 1))
                    this%avg_coordinates(3)%values(i, j, k) = 0.125d0 * ( z(i, j, k) + z (i + 1, j, k) + z(i + 1, j + 1, k)&
                        + z(i, j + 1, k) + z(i, j, k + 1) + z(i + 1, j, k + 1) + z(i + 1, j + 1, k + 1) + z(i, j + 1, k + 1))
                end do
            end do
        end do
    end subroutine Calculate_average_coordinates_3d
   

    subroutine Set_communication_mesh_3d(this, comm, comm_params)
        class (mesh_3d_t)            :: this 
        type(communication_t), pointer            :: comm
        type(communication_parameters_t), pointer :: comm_params

        call this%Set_communication_mesh_base(comm, comm_params)

        call this%avg_coordinates(1)%Set_communication(comm, comm_params)
        call this%avg_coordinates(2)%Set_communication(comm, comm_params)
        call this%avg_coordinates(3)%Set_communication(comm, comm_params)

    end subroutine Set_communication_mesh_3d

    subroutine Write_mesh_3d(this, unit, iostat, iomsg)
        class (mesh_3d_t), intent(in) :: this  
        integer,      intent(in)     :: unit
        integer,      intent(out)    :: iostat
        character(*), intent(in out)  :: iomsg

#ifdef DEBUG
        write(*,*) '@@@ in Write_mesh_3d @@@'
#endif

        call this%Write_mesh_base(unit, iostat=iostat, iomsg=iomsg)

        write(unit, iostat=iostat, iomsg=iomsg) &
            shape(this%theta), &
            shape(this%phi), &
            shape(this%start_layer_index_p), &
            size(this%n_layers_p)

        write(unit, iostat=iostat, iomsg=iomsg) &
            this%theta, &
            this%phi, &
            this%start_layer_index_p, &
            this%n_layers_p

#ifdef DEBUG
        write(*,*) &
            'shape_theta', &
            shape(this%theta), &
            'shape_phi', &
            shape(this%phi), &
            'shape_start_layer_index_p',&
            shape(this%start_layer_index_p), &
            'size_n_layers_p', &
            size(this%n_layers_p)

        write(*,*) &
            'theta', &
            this%theta, &
            'phi', &
            this%phi, &
            'start_layer_index_p',&
            this%start_layer_index_p, &
            'n_layers_p', &
            this%n_layers_p, &
            '###'

        write(*,*) '@@@ end Write_mesh_3d @@@'
#endif

    end subroutine Write_mesh_3d

    subroutine Read_mesh_3d(this, unit, iostat, iomsg)
        class (mesh_3d_t), intent(in out) :: this  
        integer,      intent(in)     :: unit
        integer,      intent(out)    :: iostat
        character(*), intent(in out)  :: iomsg
        integer, dimension(2) :: shape_theta
        integer, dimension(2) :: shape_phi
        integer, dimension(2) :: shape_start_layer_index_p
        integer :: size_n_layers_p

#ifdef DEBUG
        write(*,*) '@@@ in Read_mesh_3d @@@'
#endif

        call this%Read_mesh_base(unit, iostat=iostat, iomsg=iomsg)

        read(unit, iostat=iostat, iomsg=iomsg) &
            shape_theta, &
            shape_phi, &
            shape_start_layer_index_p, &
            size_n_layers_p

        deallocate(this%theta)
        allocate(this%theta( &
            shape_theta(1), &
            shape_theta(2)))
        deallocate(this%phi)
        allocate(this%phi( &
            shape_phi(1), &
            shape_phi(2)))
        deallocate(this%start_layer_index_p)
        allocate(this%start_layer_index_p( &
            shape_start_layer_index_p(1), &
            shape_start_layer_index_p(2)))
        deallocate(this%n_layers_p)
        allocate(this%n_layers_p(size_n_layers_p))

        read(unit, iostat=iostat, iomsg=iomsg) &
            this%theta, &
            this%phi, &
            this%start_layer_index_p, &
            this%n_layers_p

#ifdef DEBUG
        write(*,*) &
            'shape_theta', &
            shape(this%theta), &
            'shape_phi', &
            shape(this%phi), &
            'shape_start_layer_index_p',&
            shape(this%start_layer_index_p), &
            'size_n_layers_p', &
            size(this%n_layers_p)

        write(*,*) &
            'theta', &
            this%theta, &
            'phi', &
            this%phi, &
            'start_layer_index_p',&
            this%start_layer_index_p, &
            'n_layers_p', &
            this%n_layers_p, &
            '###'

        write(*,*) '@@@ end Read_mesh_3d @@@'
#endif

    end subroutine Read_mesh_3d



end module mesh_3d_module
