
module velocity_module
    use vertex_quantity_module          , only : vertex_quantity_t
    use data_module                     , only : data_t
    use boundary_parameters_module      , only : boundary_parameters_t
    use omp_lib

    implicit none
    private
    public :: velocity_t

    type, extends(vertex_quantity_t) :: velocity_t
        private

        type (data_t), pointer, public :: dvelocity_x_dx  
        type (data_t), pointer, public :: dvelocity_y_dx  
        type (data_t), pointer, public :: dvelocity_z_dx  
        type (data_t), pointer, public :: dvelocity_x_dy  
        type (data_t), pointer, public :: dvelocity_y_dy  
        type (data_t), pointer, public :: dvelocity_z_dy  
        type (data_t), pointer, public :: dvelocity_x_dz  
        type (data_t), pointer, public :: dvelocity_y_dz  
        type (data_t), pointer, public :: dvelocity_z_dz  


        logical                        :: no_move_flag    
        integer                        :: no_move_index      

        integer, public                :: start_impose_sph_sym      
        integer, public                :: end_impose_sph_sym        
    contains


        procedure, public :: Calculate_derivatives

        procedure, public :: Initialize

        procedure, public :: Clean_velocity

        procedure, public :: Write_velocity
        procedure, public :: Write_quantity_abstract => Write_velocity

        procedure, public :: Read_velocity
        procedure, public :: Read_quantity_abstract => Read_velocity

        procedure, public :: Impose_spherical_symmetry

        procedure, public :: Impose_no_move_3d

    end type


    interface velocity_t

        module procedure Constructor
        module procedure Constructor_no_derivatives
    end interface velocity_t

contains

    type(velocity_t) function Constructor(initial_data, nx, ny, nz, dimension, bc, bc_params)
        use vertex_boundary_condition_module, only : vertex_boundary_condition_t, vertex_bc_wrapper_t
        use general_utils_module            , only : Get_axises_num

        real(8)                                           , intent(in) :: initial_data 
        integer                                           , intent(in) :: nx           
        integer                                           , intent(in) :: ny           
        integer                                           , intent(in) :: nz              
        integer                                           , intent(in) :: dimension
        type (vertex_bc_wrapper_t) , dimension(:), pointer, intent(in) :: bc           
        type(boundary_parameters_t), pointer, intent(in) :: bc_params

        integer                                                        :: axises_num
        Constructor%no_move_index = -1
        Constructor%no_move_flag = .false.
        axises_num = Get_axises_num (nx, ny, nz)
        if (dimension == 3) then
            allocate(Constructor%dvelocity_x_dx)
            allocate(Constructor%dvelocity_y_dx)
            allocate(Constructor%dvelocity_z_dx)
            Constructor%dvelocity_x_dx = data_t(nx, ny, nz)
            Constructor%dvelocity_y_dx = data_t(nx, ny, nz)
            Constructor%dvelocity_z_dx = data_t(nx, ny, nz)

            allocate(Constructor%dvelocity_x_dy)
            allocate(Constructor%dvelocity_y_dy)
            allocate(Constructor%dvelocity_z_dy)
            Constructor%dvelocity_x_dy = data_t(nx, ny, nz)
            Constructor%dvelocity_y_dy = data_t(nx, ny, nz)
            Constructor%dvelocity_z_dy = data_t(nx, ny, nz)

            allocate(Constructor%dvelocity_x_dz)
            allocate(Constructor%dvelocity_y_dz)
            allocate(Constructor%dvelocity_z_dz)
            Constructor%dvelocity_x_dz = data_t(nx, ny, nz)
            Constructor%dvelocity_y_dz = data_t(nx, ny, nz)
            Constructor%dvelocity_z_dz = data_t(nx, ny, nz)
        end if
        call Constructor%Init_vertex_quantity_init_val (initial_data, nx, ny, nz, axises_num, bc, bc_params)
    end function

    type(velocity_t) function Constructor_no_derivatives(initial_data, nx, ny, nz, bc, bc_params)
        use vertex_boundary_condition_module, only : vertex_boundary_condition_t, vertex_bc_wrapper_t
        use general_utils_module            , only : Get_axises_num

        real(8)                                           , intent(in) :: initial_data 
        integer                                           , intent(in) :: nx           
        integer                                           , intent(in) :: ny           
        integer                                           , intent(in) :: nz           
        type (vertex_bc_wrapper_t) , dimension(:), pointer, intent(in) :: bc           
        type(boundary_parameters_t), pointer, intent(in) :: bc_params

        integer                                                        :: axises_num
        Constructor_no_derivatives%no_move_index = -1
        Constructor_no_derivatives%no_move_flag = .false.
        axises_num = Get_axises_num (nx, ny, nz)

        call Constructor_no_derivatives%Init_vertex_quantity_init_val (initial_data, nx, ny, nz, axises_num, bc, bc_params)
    end function


    subroutine Clean_velocity (this)
        class (velocity_t), intent(in out) :: this 

        call this%Clean_vertex_quantity ()
    end subroutine Clean_velocity


    subroutine Calculate_derivatives(this, coordinates, total_vof, total_volume, nx, ny, nz, emf, use_GPU)
        use datafile_module, only : datafile_t
        use geometry_module, only : Vector_grad
        use geometry_module, only : Vector_grad_planes, Vector_grad_vec
        use volume_module                   , only : volume_t
        use vof_module                      , only : vof_t
        use coordinates_module                      , only : coordinates_t
        implicit none       
        class (velocity_t), intent(in out) :: this
        class(coordinates_t), intent(in out) :: coordinates
        type(vof_t)       , intent(in)     :: total_vof
        type(volume_t)    , intent(in)     :: total_volume

        integer, intent(in) :: nx, ny, nz
        real(8), intent(in) :: emf    
        logical, intent(in) :: use_GPU     

        real(8), dimension(:, :, :), pointer :: dvel_x_dx   
        real(8), dimension(:, :, :), pointer :: dvel_x_dy   
        real(8), dimension(:, :, :), pointer :: dvel_x_dz   
        real(8), dimension(:, :, :), pointer :: dvel_y_dx   
        real(8), dimension(:, :, :), pointer :: dvel_y_dy   
        real(8), dimension(:, :, :), pointer :: dvel_y_dz   
        real(8), dimension(:, :, :), pointer :: dvel_z_dx   
        real(8), dimension(:, :, :), pointer :: dvel_z_dy   
        real(8), dimension(:, :, :), pointer :: dvel_z_dz   

        real(8), dimension(:, :, :), pointer :: velocity_x   
        real(8), dimension(:, :, :), pointer :: velocity_y   
        real(8), dimension(:, :, :), pointer :: velocity_z   
        real(8), dimension(:, :, :), pointer :: x            
        real(8), dimension(:, :, :), pointer :: y            
        real(8), dimension(:, :, :), pointer :: z            

        real(8), dimension(:, :, :), pointer :: vol          
        real(8), dimension(:, :, :), pointer :: vof          

        real(8) :: x1, x2, x3, x4, x5, x6, x7, x8  
        real(8) :: y1, y2, y3, y4, y5, y6, y7, y8  
        real(8) :: z1, z2, z3, z4, z5, z6, z7, z8  

        real(8) :: u1, u2, u3, u4, u5, u6, u7, u8  
        real(8) :: v1, v2, v3, v4, v5, v6, v7, v8  
        real(8) :: w1, w2, w3, w4, w5, w6, w7, w8  
        real(8) :: start_time

        integer :: i, j, k  
        integer :: ip,im,jp,jm, km,kp
        integer,save :: cntr = 0


        real(8) :: y4y1z2z1z4z1y2y1, y5y1z4z1z5z1y4y1, y2y1z5z1z2z1y5y1, &
            y1y2z3z2z1z2y3y2, y6y2z1z2z6z2y1y2, y3y2z6z2z3z2y6y2, &
            y2y3z4z3z2z3y4y3, y7y3z2z3z7z3y2y3, y4y3z7z3z4z3y7y3, &
            y3y4z1z4z3z4y1y4, y8y4z3z4z8z4y3y4, y1y4z8z4z1z4y8y4, &
            y6y5z8z5z6z5y8y5, y1y5z6z5z1z5y6y5, y8y5z1z5z8z5y1y5, &
            y7y6z5z6z7z6y5y6, y2y6z7z6z2z6y7y6, y5y6z2z6z5z6y2y6, &
            y8y7z6z7z8z7y6y7, y3y7z8z7z3z7y8y7, y6y7z3z7z6z7y3y7, &
            y5y8z7z8z5z8y7y8, y4y8z5z8z4z8y5y8, y7y8z4z8z7z8y4y8

        real(8) :: z4z1x2x1x4x1z2z1, z5z1x4x1x5x1z4z1, z2z1x5x1x2x1z5z1, &
            z1z2x3x2x1x2z3z2, z6z2x1x2x6x2z1z2, z3z2x6x2x3x2z6z2, &
            z2z3x4x3x2x3z4z3, z7z3x2x3x7x3z2z3, z4z3x7x3x4x3z7z3, &
            z3z4x1x4x3x4z1z4, z8z4x3x4x8x4z3z4, z1z4x8x4x1x4z8z4, &
            z6z5x8x5x6x5z8z5, z1z5x6x5x1x5z6z5, z8z5x1x5x8x5z1z5, &
            z7z6x5x6x7x6z5z6, z2z6x7x6x2x6z7z6, z5z6x2x6x5x6z2z6, &
            z8z7x6x7x8x7z6z7, z3z7x8x7x3x7z8z7, z6z7x3x7x6x7z3z7, &
            z5z8x7x8x5x8z7z8, z4z8x5x8x4x8z5z8, z7z8x4x8x7x8z4z8

        real(8) :: x4x1y2y1y4y1x2x1, x5x1y4y1y5y1x4x1, x2x1y5y1y2y1x5x1, &
            x1x2y3y2y1y2x3x2, x6x2y1y2y6y2x1x2, x3x2y6y2y3y2x6x2, &
            x2x3y4y3y2y3x4x3, x7x3y2y3y7y3x2x3, x4x3y7y3y4y3x7x3, &
            x3x4y1y4y3y4x1x4, x8x4y3y4y8y4x3x4, x1x4y8y4y1y4x8x4, &
            x6x5y8y5y6y5x8x5, x1x5y6y5y1y5x6x5, x8x5y1y5y8y5x1x5, &
            x7x6y5y6y7y6x5x6, x2x6y7y6y2y6x7x6, x5x6y2y6y5y6x2x6, &
            x8x7y6y7y8y7x6x7, x3x7y8y7y3y7x8x7, x6x7y3y7y6y7x3x7, &
            x5x8y7y8y5y8x7x8, x4x8y5y8y4y8x5x8, x7x8y4y8y7y8x4x8
        
        real(8), dimension(24) :: tmp_arr_y
        real(8), dimension(24) :: tmp_arr_z
        real(8), dimension(24) :: tmp_arr_x   
        real(8), dimension(24) :: co_arr
        integer :: ndev
        logical :: do_offload
        logical :: is_accel
        integer :: SIZE

        call this%dvelocity_x_dx%Point_to_data(dvel_x_dx)
        call this%dvelocity_x_dy%Point_to_data(dvel_x_dy)
        call this%dvelocity_x_dz%Point_to_data(dvel_x_dz)
        call this%dvelocity_y_dx%Point_to_data(dvel_y_dx)
        call this%dvelocity_y_dy%Point_to_data(dvel_y_dy)
        call this%dvelocity_y_dz%Point_to_data(dvel_y_dz)
        call this%dvelocity_z_dx%Point_to_data(dvel_z_dx)
        call this%dvelocity_z_dy%Point_to_data(dvel_z_dy)
        call this%dvelocity_z_dz%Point_to_data(dvel_z_dz)

        call this        %Point_to_data(velocity_x, velocity_y, velocity_z)
        call coordinates %Point_to_data(x, y, z)
        call total_vof   %Point_to_data(vof)
        call total_volume%Point_to_data(vol)
        start_time = omp_get_wtime()
      
       !SIZE = 10000 !should be updated from the datafile
       !ndev = omp_get_num_devices()
       !do_offload = (ndev > 0 .and. nz*ny*nx > SIZE)
       do_offload=use_GPU
       is_accel = .TRUE.
       !$omp target if(do_offload) map(is_accel) map(to:x, y, z, velocity_x, velocity_y, velocity_z, vof, vol) map(from:dvel_x_dx, dvel_x_dy, dvel_x_dz, dvel_y_dx, dvel_y_dy, dvel_y_dz, dvel_z_dx, dvel_z_dy, dvel_z_dz)
          IF (omp_is_initial_device()) is_accel = .FALSE.
       !$omp teams distribute parallel do collapse(3) num_threads(4) private(ip,jp,kp,x1,x2,x3,x4,x5,x6,x7,x8,y1,y2,y3,y4,y5,y6,y7,y8,z1,z2,z3,z4,z5,z6,z7,z8,u1,u2,u3,u4,u5,u6,u7,u8) &
       !$omp& private(v1,v2,v3,v4,v5,v6,v7,v8,w1,w2,w3,w4,w5,w6,w7,w8, tmp_arr_x, tmp_arr_y, tmp_arr_z, co_arr) 
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    ip = i + 1
                    jp = j + 1
                    kp = k + 1
                    !if (vof(i,j,k) < emf) cycle 
 
                     x1 = x(i, j, k)
                     x2 = x(ip, j, k)
                     x3 = x(ip, jp, k)
                     x4 = x(i, jp, k)
                     x5 = x(i, j, kp)
                     x6 = x(ip, j, kp)
                     x7 = x(ip, jp, kp)
                     x8 = x(i, jp, kp)
 
                     y1 = y(i, j, k)
                     y2 = y(ip, j, k)
                     y3 = y(ip, jp, k)
                     y4 = y(i, jp, k)
                     y5 = y(i, j, kp)
                     y6 = y(ip, j, kp)
                     y7 = y(ip, jp, kp)
                     y8 = y(i, jp, kp)
 
                     z1 = z(i, j, k)
                     z2 = z(ip, j, k)
                     z3 = z(ip, jp, k)
                     z4 = z(i, jp, k)
                     z5 = z(i, j, kp)
                     z6 = z(ip, j, kp)
                     z7 = z(ip, jp, kp)
                     z8 = z(i, jp, kp)
 
                     u1 = velocity_x(i, j, k)
                     u2 = velocity_x(ip, j, k)
                     u3 = velocity_x(ip, jp, k)
                     u4 = velocity_x(i, jp, k)
                     u5 = velocity_x(i, j, kp)
                     u6 = velocity_x(ip, j, kp)
                     u7 = velocity_x(ip, jp, kp)
                     u8 = velocity_x(i, jp, kp)
 
                     v1 = velocity_y(i, j, k)
                     v2 = velocity_y(ip, j, k)
                     v3 = velocity_y(ip, jp, k)
                     v4 = velocity_y(i, jp, k)
                     v5 = velocity_y(i, j, kp)
                     v6 = velocity_y(ip, j, kp)
                     v7 = velocity_y(ip, jp, kp)
                     v8 = velocity_y(i, jp, kp)
 
                     w1 = velocity_z(i, j, k)
                     w2 = velocity_z(ip, j, k)
                     w3 = velocity_z(ip, jp, k)
                     w4 = velocity_z(i, jp, k)
                     w5 = velocity_z(i, j, kp)
                     w6 = velocity_z(ip, j, kp)
                     w7 = velocity_z(ip, jp, kp)
                     w8 = velocity_z(i, jp, kp)
 
 
                     call Vector_grad_planes(y(i, j, k), y(ip, j, k), y(ip, jp, k), y(i, jp, k), y(i, j, kp), y(ip, j, kp), y(ip, jp, kp), y(i, jp, kp), &
                         z(i, j, k), z(ip, j, k), z(ip, jp, k), z(i, jp, k), z(i, j, kp), z(ip, j, kp), z(ip, jp, kp), z(i, jp, kp), tmp_arr_y)
 
                     call Vector_grad_planes(z(i, j, k), z(ip, j, k), z(ip, jp, k), z(i, jp, k), z(i, j, kp), z(ip, j, kp), z(ip, jp, kp), z(i, jp, kp), &
                          x(i, j, k), x(ip, j, k), x(ip, jp, k), x(i, jp, k), x(i, j, kp), x(ip, j, kp), x(ip, jp, kp), x(i, jp, kp), tmp_arr_z)   
                     
                     call Vector_grad_planes(x(i, j, k), x(ip, j, k), x(ip, jp, k), x(i, jp, k), x(i, j, kp), x(ip, j, kp), x(ip, jp, kp), x(i, jp, kp), &
                          y(i, j, k), y(ip, j, k), y(ip, jp, k), y(i, jp, k), y(i, j, kp), y(ip, j, kp), y(ip, jp, kp), y(i, jp, kp), tmp_arr_x)
 
                     call Vector_grad_vec(u1, u2, u3, u4, u5, u6, u7, u8, co_arr)
 
                     dvel_x_dx(i, j, k) = 1d0 / (12d0 * vol(i, j, k)) * &
                           ((co_arr(1) * tmp_arr_y(1) + co_arr(2) * tmp_arr_y(2) + co_arr(3) * tmp_arr_y(3)) +&
                           (co_arr(4) * tmp_arr_y(4) + co_arr(5) * tmp_arr_y(5) + co_arr(6) * tmp_arr_y(6)) +&
                           (co_arr(7) * tmp_arr_y(7) + co_arr(8) * tmp_arr_y(8) + co_arr(9) * tmp_arr_y(9)) +&
                           (co_arr(10) * tmp_arr_y(10) + co_arr(11) * tmp_arr_y(11) + co_arr(12) * tmp_arr_y(12)) +&
                           (co_arr(13) * tmp_arr_y(13) + co_arr(14) * tmp_arr_y(14) + co_arr(15) * tmp_arr_y(15)) +&
                           (co_arr(16) * tmp_arr_y(16) + co_arr(17) * tmp_arr_y(17) + co_arr(18) * tmp_arr_y(18)) +&
                           (co_arr(19) * tmp_arr_y(19) + co_arr(20) * tmp_arr_y(20) + co_arr(21) * tmp_arr_y(21)) +&
                           (co_arr(22) * tmp_arr_y(22) + co_arr(23) * tmp_arr_y(23) + co_arr(24) * tmp_arr_y(24)))
 
 
                     dvel_x_dy(i, j, k) = 1d0 / (12d0 * vol(i, j, k)) * &
                          ((co_arr(1) * tmp_arr_z(1) + co_arr(2) * tmp_arr_z(2) + co_arr(3) * tmp_arr_z(3)) +&
                          (co_arr(4) * tmp_arr_z(4) + co_arr(5) * tmp_arr_z(5) + co_arr(6) * tmp_arr_z(6)) +&
                          (co_arr(7) * tmp_arr_z(7) + co_arr(8) * tmp_arr_z(8) + co_arr(9) * tmp_arr_z(9)) +&
                          (co_arr(10) * tmp_arr_z(10) + co_arr(11) * tmp_arr_z(11) + co_arr(12) * tmp_arr_z(12)) +&
                          (co_arr(13) * tmp_arr_z(13) + co_arr(14) * tmp_arr_z(14) + co_arr(15) * tmp_arr_z(15)) +&
                          (co_arr(16) * tmp_arr_z(16) + co_arr(17) * tmp_arr_z(17) + co_arr(18) * tmp_arr_z(18)) +&
                          (co_arr(19) * tmp_arr_z(19) + co_arr(20) * tmp_arr_z(20) + co_arr(21) * tmp_arr_z(21)) +&
                          (co_arr(22) * tmp_arr_z(22) + co_arr(23) * tmp_arr_z(23) + co_arr(24) * tmp_arr_z(24)))
 
                     dvel_x_dz(i, j, k) = 1d0 / (12d0 * vol(i, j, k)) *  &
                          ((co_arr(1) * tmp_arr_x(1) + co_arr(2) * tmp_arr_x(2) + co_arr(3) * tmp_arr_x(3)) +&
                          (co_arr(4) * tmp_arr_x(4) + co_arr(5) * tmp_arr_x(5) + co_arr(6) * tmp_arr_x(6)) +&
                          (co_arr(7) * tmp_arr_x(7) + co_arr(8) * tmp_arr_x(8) + co_arr(9) * tmp_arr_x(9)) +&
                          (co_arr(10) * tmp_arr_x(10) + co_arr(11) * tmp_arr_x(11) + co_arr(12) * tmp_arr_x(12)) +&
                          (co_arr(13) * tmp_arr_x(13) + co_arr(14) * tmp_arr_x(14) + co_arr(15) * tmp_arr_x(15)) +&
                          (co_arr(16) * tmp_arr_x(16) + co_arr(17) * tmp_arr_x(17) + co_arr(18) * tmp_arr_x(18)) +&
                          (co_arr(19) * tmp_arr_x(19) + co_arr(20) * tmp_arr_x(20) + co_arr(21) * tmp_arr_x(21)) +&
                          (co_arr(22) * tmp_arr_x(22) + co_arr(23) * tmp_arr_x(23) + co_arr(24) * tmp_arr_x(24)))
  
 
                     call Vector_grad_vec(v1, v2, v3, v4, v5, v6, v7, v8, co_arr)
 
 
                     dvel_y_dx(i, j, k) = 1d0 / (12d0 * vol(i, j, k)) *  &
                         ((co_arr(1) * tmp_arr_y(1) + co_arr(2) * tmp_arr_y(2) + co_arr(3) * tmp_arr_y(3)) +&
                         (co_arr(4) * tmp_arr_y(4) + co_arr(5) * tmp_arr_y(5) + co_arr(6) * tmp_arr_y(6)) +&
                         (co_arr(7) * tmp_arr_y(7) + co_arr(8) * tmp_arr_y(8) + co_arr(9) * tmp_arr_y(9)) +&
                         (co_arr(10) * tmp_arr_y(10) + co_arr(11) * tmp_arr_y(11) + co_arr(12) * tmp_arr_y(12)) +&
                         (co_arr(13) * tmp_arr_y(13) + co_arr(14) * tmp_arr_y(14) + co_arr(15) * tmp_arr_y(15)) +&
                         (co_arr(16) * tmp_arr_y(16) + co_arr(17) * tmp_arr_y(17) + co_arr(18) * tmp_arr_y(18)) +&
                         (co_arr(19) * tmp_arr_y(19) + co_arr(20) * tmp_arr_y(20) + co_arr(21) * tmp_arr_y(21)) +&
                         (co_arr(22) * tmp_arr_y(22) + co_arr(23) * tmp_arr_y(23) + co_arr(24) * tmp_arr_y(24)))
 
 
                     dvel_y_dy(i, j, k) = 1d0 / (12d0 * vol(i, j, k)) *  &
                         ((co_arr(1) * tmp_arr_z(1) + co_arr(2) * tmp_arr_z(2) + co_arr(3) * tmp_arr_z(3)) +&
                         (co_arr(4) * tmp_arr_z(4) + co_arr(5) * tmp_arr_z(5) + co_arr(6) * tmp_arr_z(6)) +&
                         (co_arr(7) * tmp_arr_z(7) + co_arr(8) * tmp_arr_z(8) + co_arr(9) * tmp_arr_z(9)) +&
                         (co_arr(10) * tmp_arr_z(10) + co_arr(11) * tmp_arr_z(11) + co_arr(12) * tmp_arr_z(12)) +&
                         (co_arr(13) * tmp_arr_z(13) + co_arr(14) * tmp_arr_z(14) + co_arr(15) * tmp_arr_z(15)) +&
                         (co_arr(16) * tmp_arr_z(16) + co_arr(17) * tmp_arr_z(17) + co_arr(18) * tmp_arr_z(18)) +&
                         (co_arr(19) * tmp_arr_z(19) + co_arr(20) * tmp_arr_z(20) + co_arr(21) * tmp_arr_z(21)) +&
                         (co_arr(22) * tmp_arr_z(22) + co_arr(23) * tmp_arr_z(23) + co_arr(24) * tmp_arr_z(24)))
 
 
                     dvel_y_dz(i, j, k) = 1d0 / (12d0 * vol(i, j, k)) *  &
                         ((co_arr(1) * tmp_arr_x(1) + co_arr(2) * tmp_arr_x(2) + co_arr(3) * tmp_arr_x(3)) +&
                         (co_arr(4) * tmp_arr_x(4) + co_arr(5) * tmp_arr_x(5) + co_arr(6) * tmp_arr_x(6)) +&
                         (co_arr(7) * tmp_arr_x(7) + co_arr(8) * tmp_arr_x(8) + co_arr(9) * tmp_arr_x(9)) +&
                         (co_arr(10) * tmp_arr_x(10) + co_arr(11) * tmp_arr_x(11) + co_arr(12) * tmp_arr_x(12)) +&
                         (co_arr(13) * tmp_arr_x(13) + co_arr(14) * tmp_arr_x(14) + co_arr(15) * tmp_arr_x(15)) +&
                         (co_arr(16) * tmp_arr_x(16) + co_arr(17) * tmp_arr_x(17) + co_arr(18) * tmp_arr_x(18)) +&
                         (co_arr(19) * tmp_arr_x(19) + co_arr(20) * tmp_arr_x(20) + co_arr(21) * tmp_arr_x(21)) +&
                         (co_arr(22) * tmp_arr_x(22) + co_arr(23) * tmp_arr_x(23) + co_arr(24) * tmp_arr_x(24)))
 
 
 
                     call Vector_grad_vec(w1, w2, w3, w4, w5, w6, w7, w8, co_arr)
 
                     dvel_z_dx(i, j, k) = 1d0 / (12d0 * vol(i, j, k)) *  &
                           ((co_arr(1) * tmp_arr_y(1) + co_arr(2) * tmp_arr_y(2) + co_arr(3) * tmp_arr_y(3)) +&
                           (co_arr(4) * tmp_arr_y(4) + co_arr(5) * tmp_arr_y(5) + co_arr(6) * tmp_arr_y(6)) +&
                           (co_arr(7) * tmp_arr_y(7) + co_arr(8) * tmp_arr_y(8) + co_arr(9) * tmp_arr_y(9)) +&
                           (co_arr(10) * tmp_arr_y(10) + co_arr(11) * tmp_arr_y(11) + co_arr(12) * tmp_arr_y(12)) +&
                           (co_arr(13) * tmp_arr_y(13) + co_arr(14) * tmp_arr_y(14) + co_arr(15) * tmp_arr_y(15)) +&
                           (co_arr(16) * tmp_arr_y(16) + co_arr(17) * tmp_arr_y(17) + co_arr(18) * tmp_arr_y(18)) +&
                           (co_arr(19) * tmp_arr_y(19) + co_arr(20) * tmp_arr_y(20) + co_arr(21) * tmp_arr_y(21)) +&
                           (co_arr(22) * tmp_arr_y(22) + co_arr(23) * tmp_arr_y(23) + co_arr(24) * tmp_arr_y(24)))
 
 
                     dvel_z_dy(i, j, k) = 1d0 / (12d0 * vol(i, j, k)) *  &
                           ((co_arr(1) * tmp_arr_z(1) + co_arr(2) * tmp_arr_z(2) + co_arr(3) * tmp_arr_z(3)) +&
                           (co_arr(4) * tmp_arr_z(4) + co_arr(5) * tmp_arr_z(5) + co_arr(6) * tmp_arr_z(6)) +&
                           (co_arr(7) * tmp_arr_z(7) + co_arr(8) * tmp_arr_z(8) + co_arr(9) * tmp_arr_z(9)) +&
                           (co_arr(10) * tmp_arr_z(10) + co_arr(11) * tmp_arr_z(11) + co_arr(12) * tmp_arr_z(12)) +&
                           (co_arr(13) * tmp_arr_z(13) + co_arr(14) * tmp_arr_z(14) + co_arr(15) * tmp_arr_z(15)) +&
                           (co_arr(16) * tmp_arr_z(16) + co_arr(17) * tmp_arr_z(17) + co_arr(18) * tmp_arr_z(18)) +&
                           (co_arr(19) * tmp_arr_z(19) + co_arr(20) * tmp_arr_z(20) + co_arr(21) * tmp_arr_z(21)) +&
                           (co_arr(22) * tmp_arr_z(22) + co_arr(23) * tmp_arr_z(23) + co_arr(24) * tmp_arr_z(24)))
 
                     dvel_z_dz(i, j, k) = 1d0 / (12d0 * vol(i, j, k)) *  &

                           ((co_arr(1) * tmp_arr_x(1) + co_arr(2) * tmp_arr_x(2) + co_arr(3) * tmp_arr_x(3)) +&
                           (co_arr(4) * tmp_arr_x(4) + co_arr(5) * tmp_arr_x(5) + co_arr(6) * tmp_arr_x(6)) +&
                           (co_arr(7) * tmp_arr_x(7) + co_arr(8) * tmp_arr_x(8) + co_arr(9) * tmp_arr_x(9)) +&
                           (co_arr(10) * tmp_arr_x(10) + co_arr(11) * tmp_arr_x(11) + co_arr(12) * tmp_arr_x(12)) +&
                           (co_arr(13) * tmp_arr_x(13) + co_arr(14) * tmp_arr_x(14) + co_arr(15) * tmp_arr_x(15)) +&
                           (co_arr(16) * tmp_arr_x(16) + co_arr(17) * tmp_arr_x(17) + co_arr(18) * tmp_arr_x(18)) +&
                           (co_arr(19) * tmp_arr_x(19) + co_arr(20) * tmp_arr_x(20) + co_arr(21) * tmp_arr_x(21)) +&
                           (co_arr(22) * tmp_arr_x(22) + co_arr(23) * tmp_arr_x(23) + co_arr(24) * tmp_arr_x(24)))

                end do
            end do
        end do      
        !$omp end teams distribute parallel do    
        !$omp end target
        write(*,*) "Calculate_derivatives executed on accelerator:", is_accel
        write(*,*) "Calculate_derivatives time: ", omp_get_wtime() - start_time
        !$omp target update if(do_offload) from(dvel_x_dx, dvel_x_dy, dvel_x_dz, dvel_y_dx, dvel_y_dy, dvel_y_dz, dvel_z_dx, dvel_z_dy, dvel_z_dz)
        cntr = cntr  + 1
    end subroutine Calculate_derivatives

    subroutine Impose_spherical_symmetry(this, coords)                                                                                     
        use coordinates_module                      , only : coordinates_t
        implicit none
        class (velocity_t)   , intent(inout)  :: this 
        class (coordinates_t), intent(inout)  :: coords

        real(8), dimension(:, :, :), pointer :: velocity_x   
        real(8), dimension(:, :, :), pointer :: velocity_y   
        real(8), dimension(:, :, :), pointer :: velocity_z   
        real(8), dimension(:, :, :), pointer :: x            
        real(8), dimension(:, :, :), pointer :: y            
        real(8), dimension(:, :, :), pointer :: z            

        real(8) :: vel_x_tmp 
        real(8) :: vel_y_tmp 
        real(8) :: vel_z_tmp 

        real(8), dimension(:), allocatable :: vel_radial_vec            
        integer, dimension(:), allocatable :: index_vertex_vec          

        integer :: nxp, nyp, nzp, virt_k_start, tmp_i, k1,i,j,k
        integer :: virt_nxp, virt_nyp, virt_nzp
        integer :: nzp_rad 
        integer :: from_index, to_index

        virt_nzp = this%parallel_params%virt_nzp 
        call this%parallel_params%get_virtual_from_to(virt_k_start, tmp_i, 3)
        nxp = this%d1
        nyp = this%d2
        nzp = this%d3
        from_index = this%start_impose_sph_sym
        to_index   = this%end_impose_sph_sym

        if ( (from_index <= 0) .and. (to_index > nzp)) return 

        allocate(vel_radial_vec  (1:virt_nzp), source=0d0)
        allocate(index_vertex_vec(1:virt_nzp), source=0)
        call this  %Point_to_data(velocity_x, velocity_y, velocity_z)
        call coords%Point_to_data(x, y, z)
        do k = 1, nzp
            k1 = virt_k_start + k
            do j = 1, nyp
                do i = 1, nxp
                            vel_x_tmp = velocity_x(i, j, k)
                            vel_y_tmp = velocity_y(i, j, k)
                            vel_z_tmp = velocity_z(i, j, k)
                        vel_radial_vec(k1) = vel_radial_vec(k1) + (vel_x_tmp * x(i, j, k)  &
                                                                +  vel_y_tmp * y(i, j, k)  &
                                                                +  vel_z_tmp * z(i, j, k)) &
                                                                / (x(i,j,k) * x(i,j,k) + y(i,j,k)*y(i,j,k) + z(i,j,k) * z(i,j,k))

                        index_vertex_vec(k1) = index_vertex_vec(k1) + 1

                end do
            end do
        end do

        do k = 1, nzp
            k1 = virt_k_start + k
            if (k <= from_index .or. k >= to_index) then
                vel_radial_vec(k1) = vel_radial_vec(k1) / dble(index_vertex_vec(k1))
                do j = 1, nyp
                    do i = 1, nxp
                        velocity_x(i, j, k) = vel_radial_vec(k1) * x(i, j, k)
                        velocity_y(i, j, k) = vel_radial_vec(k1) * y(i, j, k)
                        velocity_z(i, j, k) = vel_radial_vec(k1) * z(i, j, k)
                    end do
                end do
            end if
        end do
        deallocate(vel_radial_vec)
        deallocate(index_vertex_vec)
    end subroutine Impose_spherical_symmetry

    subroutine Impose_no_move_3d(this, coords)
        use coordinates_module                      , only : coordinates_t
        implicit none
        class (velocity_t)   , intent(inout)  :: this 
        class (coordinates_t), intent(inout)  :: coords

        integer :: i, j, k, nxp, nyp, nzp
        real(8) :: radial_vel 
        real(8) :: vel_x_tmp 
        real(8) :: vel_y_tmp 
        real(8) :: vel_z_tmp 
        real(8) :: vel
        real(8), dimension(:, :, :), pointer :: velocity_x   
        real(8), dimension(:, :, :), pointer :: velocity_y   
        real(8), dimension(:, :, :), pointer :: velocity_z   
        real(8), dimension(:, :, :), pointer :: x            
        real(8), dimension(:, :, :), pointer :: y            
        real(8), dimension(:, :, :), pointer :: z            

        if (this%no_move_flag .eqv. .false.) return

        call this  %Point_to_data(velocity_x, velocity_y, velocity_z)
        call coords%Point_to_data(x, y, z)
        radial_vel = 1d0
        nxp = this%d1
        nyp = this%d2
        nzp = this%d3
        k = this%no_move_index
        do j = 1, nyp
            do i = 1, nxp
                          vel_x_tmp = velocity_x(i, j, k)
                          vel_y_tmp = velocity_y(i, j, k)
                          vel_z_tmp = velocity_z(i, j, k)
                vel = (vel_x_tmp*x(i,j,k) + vel_y_tmp*y(i,j,k) +vel_z_tmp*z(i,j,k)) &
                      / (x(i,j,k)*x(i,j,k) + y(i,j,k)*y(i,j,k) + z(i,j,k)*z(i,j,k))
                radial_vel = min(radial_vel, vel)
            end do
        end do

        if (radial_vel < -1d-10) then
            this%no_move_flag = .false.
        end if
        if (this%no_move_flag .eqv. .false.) return

        do j = 1, nyp
            do i = 1, nxp
                velocity_x(i, j, k) = 0d0
                velocity_y(i, j, k) = 0d0
                velocity_z(i, j, k) = 0d0
            end do
        end do
    end subroutine Impose_no_move_3d

    subroutine Initialize(this,  i_sphere, i_sphere_up, no_move_layer, start_index)
        use parallel_parameters_module      , only : parallel_parameters_t
        implicit none
        class (velocity_t)   , intent(inout)             :: this 
        integer              , intent(in)                :: no_move_layer 
        integer                             , intent(in) :: i_sphere             
        integer                             , intent(in) :: i_sphere_up          
        integer, dimension(:,:), allocatable, intent(in) :: start_index

        integer :: ind, no_move_index
        integer :: k_start, k_end, inside, i_sphere2

        if (no_move_layer < 0) then
            this%no_move_index = -1
            this%no_move_flag = .false.
        else
            no_move_index = start_index(no_move_layer, 1)

            this%no_move_index = no_move_index
            inside = this%parallel_params%Inside_physical_world(no_move_index, axis=3)
            select case (inside)
                case (1)
                    this%no_move_flag = .false.
                case (0)
                    this%no_move_flag = .true.
                case(-1)
                    this%no_move_flag = .false.
            end select
            this%no_move_index = this%parallel_params%get_physical_index(no_move_index, 3)
        end if

        inside = this%parallel_params%Inside_physical_world(i_sphere, axis=3)
        select case (inside)
            case (1)
                this%end_impose_sph_sym = this%d3
            case (0)
                this%end_impose_sph_sym = this%parallel_params%get_physical_index(i_sphere, 3)
            case (-1)
                this%end_impose_sph_sym = huge(1)
        end select

        i_sphere2 = this%parallel_params%virt_nzp + 1 - i_sphere_up
        inside = this%parallel_params%Inside_physical_world(i_sphere2, axis=3)
        select case (inside)
            case (1)
                this%start_impose_sph_sym = -1
            case (0)
                this%start_impose_sph_sym = this%parallel_params%get_physical_index(i_sphere2, 3)
            case (-1)
                this%start_impose_sph_sym = 1
        end select
    end subroutine Initialize

    subroutine Write_velocity(this, unit, iostat, iomsg)
        class (velocity_t), intent(in) :: this  
        integer,      intent(in)     :: unit
        integer,      intent(out)    :: iostat
        character(*), intent(in out)  :: iomsg

#ifdef DEBUG
        write(*,*) '@@@ in Write_velocity @@@'
#endif

        call this%Write_vertex_quantity(unit, iostat=iostat, iomsg=iomsg)

#ifdef DEBUG
        write(*,*) '@@@ end Write_velocity @@@'
#endif

    end subroutine Write_velocity

    subroutine Read_velocity(this, unit, iostat, iomsg)
        class (velocity_t), intent(in out) :: this  
        integer,      intent(in)     :: unit
        integer,      intent(out)    :: iostat
        character(*), intent(in out)  :: iomsg

#ifdef DEBUG
        write(*,*) '@@@ in Read_velocity @@@'
#endif

        call this%Read_vertex_quantity(unit, iostat=iostat, iomsg=iomsg)

#ifdef DEBUG
        write(*,*) '@@@ end Read_velocity @@@'
#endif

    end subroutine Read_velocity

end module velocity_module
