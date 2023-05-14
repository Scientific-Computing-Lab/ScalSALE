
module mesh_module
use data_module, only : data_t
use datafile_module, only : datafile_t
use mesh_base_module, only : mesh_base_t
    use coordinates_module              , only : coordinates_t
    use communication_parameters_module , only : communication_parameters_t
    use parallel_parameters_module      , only : parallel_parameters_t
    use boundary_parameters_module      , only : boundary_parameters_t
    use vertex_boundary_condition_module, only : vertex_boundary_condition_t, vertex_bc_wrapper_t
    use communication_module            , only : communication_t
implicit none
    private



    type, extends(mesh_base_t), public :: mesh_t
        private

        real(8), dimension(:)   , allocatable :: theta               

    contains

        procedure, public    :: Calculate_angles

        procedure, public :: Set_communication_mesh

        procedure, public :: Set_communication => Set_communication_mesh

        procedure, private    :: Read_and_init_mesh1           

        procedure, private   :: Build_mesh
        procedure, private   :: Build_coordinates
        procedure, private   :: Build_angles           

        procedure, public :: Write_mesh
        procedure, public    :: Write_mesh_abstract => Write_mesh

        procedure, public :: Read_mesh
        procedure, public    :: Read_mesh_abstract => Read_mesh

    end type mesh_t

    interface mesh_t

        module procedure Constructor_2d
    end interface mesh_t

contains

    type(mesh_t) function Constructor_2d(df, bc_v, bc_params, parallel_params)
        implicit none
        type(datafile_t), intent(in)     :: df
        type(parallel_parameters_t), pointer, intent(inout) :: parallel_params
        type (vertex_bc_wrapper_t) , dimension(:), pointer, intent(in out) :: bc_v           
        type(boundary_parameters_t), pointer, intent(in) :: bc_params

        integer :: nxp, nyp, nzp

        nxp = parallel_params%nxp
        nyp = parallel_params%nyp
        nzp = parallel_params%nzp
        call Constructor_2d%Init_mesh_base(df, parallel_params)
        allocate(Constructor_2d%coordinates)
        Constructor_2d%coordinates = coordinates_t(0d0, nxp + 1, nyp + 1, 1,bc_v, bc_params)



        allocate (Constructor_2d%theta (nyp + 1))
        Constructor_2d%dimension      = 2

        call Constructor_2d%Read_and_init_mesh1(df)
    end function
    
    subroutine Calculate_angles(this)
        implicit none
        class (mesh_t)                :: this 
    end subroutine Calculate_angles

    subroutine Read_and_init_mesh1(this, df)
        use constants_module, only: PI
        implicit none
        type(datafile_t), intent(in)     :: df

        class (mesh_t), intent(in out)                  :: this            

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

        same_material_in_layer = .true.
        contour_segments = 1 


        call this%Build_mesh (df, 1, 1) 

        if (this%mesh_type  == 2) then
            this%sphere_factor = 2d0 * PI
            if (this%cyl == 0) then 
                this%sphere_factor = 1.d0
            end if
        end if
        return
    end subroutine Read_and_init_mesh1


    subroutine Build_mesh(this, df, sw_is, sw_js)
        implicit none

        type(datafile_t)  , intent(in)         :: df
        class (mesh_t)    , intent(in out)     :: this
        integer           , intent(in)         :: sw_is, sw_js

        integer       , dimension(:),allocatable        :: contour_i_type
        real(8)       , dimension(:,:),allocatable            :: contour_i_line
        real(8), dimension(:, :), allocatable           :: rr
        integer                                         :: i, j

        call this%Build_angles(df%theta0, sw_js, df%contour_j_type)
        call this%Build_coordinates (df%contour_i_type, df%contour_i_line, df%zone_i_type, df%zone_i_d, sw_is)
        this%r_factor%values = this%coordinates%data(1)%values * this%cyl + (1 - this%cyl) 
    end subroutine Build_mesh



    subroutine Build_coordinates(this, contour_i_type, contour_line, zone_type, zone_dr, sw_is)
        use geometry_module, only : Layer_cut_line,Layer_cut_elipse, Calculate_constant,Root,Calculate_geometry_series
        implicit none
        class(mesh_t)               , intent(inout)         :: this
        integer                     , intent(in)            :: sw_is
        integer       , dimension(:), intent(in)            :: contour_i_type
        integer       , dimension(:), intent(in)            :: zone_type
        real(8)       , dimension(:), intent(in)            :: zone_dr
        real(8)       , dimension(:,:), intent(in)          :: contour_line

        integer                                               :: n_part_max
        real(8)                                               :: r0, y11, x11, y12, x12, rr1, rr2, rr3, xh1, xh2&
            , xh3, yh1, yh2, yh3, teta
        real(8)                                               :: x_layer_up(1, 11), x_layer_down(1, 11), x_layer_up1(1, 11) 
        integer                                               :: i, j, j1, j2, i1, i2, mm, ll, n, j11, j22
        real(8)                                               :: dr, qq,ra1,ra2,x1,x2,y1,y2
        real(8), dimension(:), pointer :: rr
        real(8)       , dimension(:,:),allocatable          :: contour_i_line
        integer :: d1, d2

        n_part_max = 1

        allocate(contour_i_line, source=contour_line)
        call this%Fix_contour_i(contour_i_type, contour_i_line)

        allocate(rr(0:this%start_layer_index_r(this%n_layers_r(1) + 1, 1)))


        j1 = this%start_layer_index_t(1, 1) + (sw_is - 1)
        j2 = this%start_layer_index_t(this%n_layers_t(1) + 1, 1) + (sw_is - 1)
        do n = 1, this%n_layers_r(1)
            do mm = 1, n_part_max
                do ll = 1, size(contour_i_line(n, :))
                    x_layer_down(mm, ll) = contour_i_line(n, ll) 
                    x_layer_up(mm, ll)   = contour_i_line(n + 1, ll) 

                    if (n < this%n_layers_r(1)) then
                        x_layer_up1(mm, ll) = contour_i_line(n + 2, ll)
                    end if
                    i1 = this%start_layer_index_r(n, 1) + (sw_is - 1)
                    i2 = this%start_layer_index_r(n + 1, 1) + (sw_is - 1)
                end do 
            end do  

            j11 = 1
            j22 = this%nyp
            do j = j11 ,j22
                if (j < j1 .or. j > j2) continue
                teta = this%theta(j)
                r0 = 1000.0d0
                if (this%mesh_type == 2) then
                    y11 = teta
                    x11 = -r0
                    y12 = teta
                    x12 = r0
                else
                    y11 = 0.0
                    x11 = 0.0
                    y12 = r0*cos(teta)
                    x12 = r0*sin(teta)
                end if
                select case(contour_i_type(n))
                    case (0) 
                        call Layer_cut_line (y11, x11, y12, x12, x_layer_down, n_part_max, yh1, xh1)
                        call Layer_cut_line (y11, x11, y12, x12, x_layer_up  , n_part_max, yh2, xh2)

                    case (1)                  
                        call Layer_cut_elipse (y11, x11, y12, x12, x_layer_down, n_part_max, yh1, xh1)
                        call Layer_cut_elipse (y11, x11, y12, x12, x_layer_up  , n_part_max, yh2, xh2)
                    case default
                        write(*,*) "ERROR no contour type number: ", contour_i_type(n)
                end select

                if (n < this%n_layers_r(1)) then
                    select case(contour_i_type(n))
                        case (0) 
                            call Layer_cut_line (y11, x11, y12, x12, x_layer_up1, n_part_max, yh3, xh3)
                        case (1)                  
                            call Layer_cut_elipse (y11, x11, y12, x12, x_layer_up1, n_part_max, yh3, xh3)
                        case default
                            write(*,*) "ERROR no contour type number: ", contour_i_type(n)
                    end select
                end if
                if (this%mesh_type == 2) then
                    rr1 = xh1
                    rr2 = xh2
                    if (this%n_layers_r(1) > n) rr3 = xh3
                else
                    rr1 = sqrt(xh1**2 + (yh1 - 0.0)**2)
                    rr2 = sqrt(xh2**2 + (yh2 - 0.0)**2)
                    if (this%n_layers_r(1) > n) rr3 = sqrt(xh3**2 + (yh3 - 0.0)**2)
                end if

                select case(zone_type(n))
                    case (0)
                        call Calculate_constant(rr, rr1, rr2, i1, i2) 
                    case (2)
                        dr = zone_dr(n)
                        if (dr < 0d0) then
                            x1 = this%coordinates%data(1)%values(i1  , j, 1)
                            x2 = this%coordinates%data(1)%values(i1-1, j, 1)
                            y1 = this%coordinates%data(2)%values(i1  , j, 1)
                            y2 = this%coordinates%data(2)%values(i1-1, j, 1)
                            ra1 = sqrt(x1**2 + y1**2)
                            ra2 = sqrt(x2**2 + y2**2)
                            dr = ra1 - ra2
                            if (this%mesh_type == 2) then
                                dr = abs(x1 - x2)
                            end if
                            dr = dr * abs(zone_dr(n))
                        end if
                        rr(i1) = rr1
                        qq = root(1.05d0, dr, (rr2-rr1), (i2-i1))
                        call Calculate_geometry_series(rr, dr, qq, i1 + 1, i2)
                    case default
                        write(*,*) "ERROR no zone type number: ", zone_type(n)
                end select

                if (this%mesh_type == 2) then
                    do i = i1, i2
                        this%coordinates%data(1)%values(i, j, 1) = rr(i)
                        this%coordinates%data(2)%values(i, j, 1) = this%theta(j)
                    end do
                else
                    do i = i1, i2
                        this%coordinates%data(1)%values(i, j, 1) = rr(i) * sin(teta)
                        this%coordinates%data(2)%values(i, j, 1) = rr(i) * cos(teta)
                    end do
                end if


            end do 
        end do 
        deallocate(rr)
        deallocate(contour_i_line)
    end subroutine Build_coordinates

    subroutine Build_angles(this, theta0_df, sw_js, contour_j_type)
        use constants_module, only : PI
        use geometry_module, only : Calculate_constant

        implicit none
        class(mesh_t)                       , intent(inout)  :: this
        integer, dimension(:)  ,              intent(in)     :: contour_j_type
        real(8), dimension(:,:), allocatable, intent(in)     :: theta0_df
        integer                             , intent(in)     :: sw_js

        integer                                        :: i, j1, j2, n, j,nx ,ny
        real(8)                                        :: theta1, theta2 , theta_from, theta_to, theta_all
        real(8), dimension(:), pointer             :: virt_theta
        real(8), dimension(:,:), allocatable             :: theta0

        allocate(theta0, source=theta0_df)
        call this%Fix_contour_angle(theta0, contour_j_type)

        theta_from = theta0(1, 1)
        theta_to   = theta0(this%n_layers_t(1) + 1, 1)
        this%sphere_factor = 4.d0*PI/(cos(theta_to) - cos(theta_from))

        theta_all = theta_from - theta_to

        nx = this%nxp - 1
        ny = this%nyp - 1

        allocate (virt_theta(0 : this%nyp + 1))
        virt_theta = 0d0

        do n = 1, this%n_layers_t(1)
            j1 = this%start_layer_index_t(n, 1) + (sw_js - 1)
            j2 = this%start_layer_index_t(n + 1, 1) + (sw_js - 1)
            theta1 = theta0(n, 1)
            theta2 = theta0(n + 1, 1)

            call Calculate_constant(virt_theta, theta1, theta2, j1, j2)
        end do 

        j1 = this%start_layer_index_t(1, 1) + (sw_js - 1)
        j2 = this%start_layer_index_t(this%n_layers_t(1) + 1, 1) + (sw_js - 1)

        do j = 1, this%nyp + 1 
            this%theta(j) = virt_theta(j)
        end do
        deallocate(virt_theta) 
        deallocate(theta0)

    end subroutine Build_angles

    subroutine Set_communication_mesh(this, comm, comm_params)
        class (mesh_t)            :: this 
        type(communication_t), pointer            :: comm
        type(communication_parameters_t), pointer :: comm_params


        call this%Set_communication_mesh_base(comm, comm_params)


    end subroutine Set_communication_mesh

    subroutine Write_mesh(this, unit, iostat, iomsg)
        class (mesh_t), intent(in) :: this  
        integer,      intent(in)     :: unit
        integer,      intent(out)    :: iostat
        character(*), intent(in out)  :: iomsg

#ifdef DEBUG
        write(*,*) '@@@ in Write_mesh @@@'
#endif

        call this%Write_mesh_base(unit, iostat=iostat, iomsg=iomsg)

        write(unit, iostat=iostat, iomsg=iomsg) &
            this%theta

#ifdef DEBUG
        write(*,*) &
            'theta', &
            this%theta, &
            '###'

        write(*,*) '@@@ end Write_mesh @@@'
#endif

    end subroutine Write_mesh

    subroutine Read_mesh(this, unit, iostat, iomsg)
        class (mesh_t), intent(in out) :: this  
        integer,      intent(in)     :: unit
        integer,      intent(out)    :: iostat
        character(*), intent(in out)  :: iomsg

#ifdef DEBUG
        write(*,*) '@@@ in Read_mesh @@@'
#endif

        call this%Read_mesh_base(unit, iostat=iostat, iomsg=iomsg)

        read(unit, iostat=iostat, iomsg=iomsg) &
            this%theta

#ifdef DEBUG
        write(*,*) &
            'theta', &
            this%theta, &
            '###'

        write(*,*) '@@@ end Read_mesh @@@'
#endif

    end subroutine Read_mesh

end module mesh_module
