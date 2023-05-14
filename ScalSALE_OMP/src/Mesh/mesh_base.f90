

module mesh_base_module
    use data_module                     , only : data_t
    use datafile_module                 , only : datafile_t
    use coordinates_module              , only : coordinates_t
    use parallel_parameters_module      , only : parallel_parameters_t
    use communication_module            , only : communication_t
    use communication_parameters_module , only : communication_parameters_t
    use boundary_parameters_module      , only : boundary_parameters_t
    use vertex_boundary_condition_module, only : vertex_boundary_condition_t, vertex_bc_wrapper_t

    implicit none
    private

    public :: mesh_base_t

    type, abstract :: mesh_base_t
        private
        integer, public :: nxp, nyp, nzp 
        integer, public :: nx, ny, nz 
        integer, public :: dimension     

        type (coordinates_t)        , pointer, public :: coordinates    
        type (data_t), public                                 :: r_factor       

        type (communication_t)      ,pointer,public   :: mpi_comm       
        type (parallel_parameters_t), pointer, public :: parallel_params

        real(8), public :: cyl   
        real(8), public :: omcyl 

        integer, dimension(:, :), allocatable, public     :: start_layer_index_r 

        integer, dimension(:, :), allocatable, public     :: start_layer_index_t  

        integer, dimension(:)   , allocatable, public :: n_layers_r   

        integer, dimension(:)   , allocatable, public :: n_layers_t   

        real(8), public                       :: sphere_factor     

        integer, public                       :: mesh_type  

    contains

        procedure, public    :: Init_mesh_base

        procedure, public    :: Fix_contour_i

        procedure, public    :: Fix_contour_angle

        procedure, public    :: Set_communication_mesh_base

        procedure, public    :: Set_communication

        procedure, public    :: Point_to_r_factor

        procedure, public    :: Exchange_virtual_space_blocking

        procedure, public    :: Clean_mesh

        procedure, private   :: Ptr_coordinates_1d

        procedure, private   :: Ptr_coordinates_2d

        procedure, private   :: Ptr_coordinates_3d
        generic,   public    :: Point_to_data  =>     &
            Ptr_coordinates_1d,   &
            Ptr_coordinates_2d,   &
            Ptr_coordinates_3d

        procedure, public    :: Calculate_average_coordinates

        procedure, public :: Write_mesh_abstract
        procedure, public :: Write_mesh_base
        generic :: write(unformatted) => Write_mesh_base

        procedure, public :: Read_mesh_abstract
        procedure, public :: Read_mesh_base
        generic :: read(unformatted) => Read_mesh_base

    end type mesh_base_t


contains


    subroutine Init_mesh_base(this, df, parallel_params)
        implicit none
        class(mesh_base_t), intent(in out) :: this
        type(datafile_t), intent(in)       :: df
        type(parallel_parameters_t), pointer, intent(inout)     :: parallel_params

        integer :: nxp, nyp, nzp
        this%dimension = df%dimension
        allocate (this%coordinates)
        nxp = parallel_params%nxp
        nyp = parallel_params%nyp
        nzp = parallel_params%nzp
        this%parallel_params => parallel_params


        this%mesh_type = df%mesh_type

        allocate (this%n_layers_r(1))
        allocate (this%n_layers_t(1))
        this%n_layers_r        = df%number_layers_i
        this%n_layers_t        = df%number_layers_j

        allocate(this%start_layer_index_r, source=df%start_layer_index_r)
        allocate(this%start_layer_index_t, source=df%start_layer_index_t)


        this%nxp           = nxp
        this%nyp           = nyp
        this%nzp           = nzp
        this%nx            = parallel_params%nx
        this%ny            = parallel_params%ny
        this%nz            = parallel_params%nz

        this%r_factor       = data_t (nxp + 1, nyp + 1, nzp)
        this%cyl            = df%cyl
        this%omcyl          = 1 - this%cyl


    end subroutine

    subroutine Set_communication_mesh_base(this, comm, comm_params)
        class (mesh_base_t),intent(in out)            :: this 
        type(communication_t), pointer            :: comm
        type(communication_parameters_t), pointer :: comm_params

        this%mpi_comm => comm
        call this%coordinates%Set_communication(comm, comm_params)
        call this%r_factor%Set_communication(comm, comm_params)

    end subroutine Set_communication_mesh_base

    subroutine Fix_contour_i(this, contour_i_type, contour_i_line)
        use geometry_module , only : Elipse_cut_line
        use constants_module, only : PI
        implicit none
        class (mesh_base_t),intent(in out)            :: this 
        integer       , dimension(:), intent(in)            :: contour_i_type
        real(8)       , dimension(:,:), intent(inout)       :: contour_i_line

        real(8) :: a,b,al,be,yc,xc,teta1,teta2,y1,x1,y2,x2,xh1,xh2,yh1,yh2
        integer :: i,n, ll
        n = size(contour_i_type)

        do i = 1, n
            if (contour_i_type(i) == 1) then 
                a  = contour_i_line(i, 5)
                b  = contour_i_line(i, 6)
                al = contour_i_line(i, 7)
                be = contour_i_line(i, 8)
                yc = contour_i_line(i, 9)
                xc = contour_i_line(i, 10)
                teta1 = contour_i_line(i, 1)
                teta2 = contour_i_line(i, 2)
                if (contour_i_line(i, 4) == 1d0) then 
                    teta1 = teta1 * acos(-1d0)
                    teta2 = teta2 * acos(-1d0)
                elseif (contour_i_line(i, 4) == 2d0) then 
                    teta1 = teta1 / 180d0 * PI
                    teta2 = teta2 / 180d0 * PI
                end if
                y1 = yc
                x1 = xc
                y2 = yc + 100.0d0 * cos(teta1)
                x2 = xc + 100.0d0 * sin(teta1)
                call Elipse_cut_line(a, b, al, be, yc, xc, y1, x1, y2, x2, yh1, xh1)
                y1 = yc
                x1 = xc
                y2 = yc + 100.0d0 * cos(teta2)
                x2 = xc + 100.0d0 * sin(teta2)
                call Elipse_cut_line(a, b, al, be, yc, xc, y1, x1, y2, x2, yh2, xh2)
                contour_i_line(i, 1) = yh1
                contour_i_line(i, 2) = xh1
                contour_i_line(i, 3) = yh2
                contour_i_line(i, 4) = xh2
            end if
        end do

    end subroutine Fix_contour_i

    subroutine Fix_contour_angle(this, virt, contour_type)
        use geometry_module , only : Elipse_cut_line
        use constants_module, only : PI
        implicit none
        class (mesh_base_t),intent(in out)            :: this 
        integer       , dimension(:), intent(in)            :: contour_type
        real(8)       , dimension(:,:), intent(inout)       :: virt

        integer :: i,n, ll
        n = size(contour_type)
        do i = 1, n
            if (contour_type(i) == 2) then      
                virt(i, 1) = virt(i, 1) * PI
            else if (contour_type(1) == 1) then
                virt(i, 1) = virt(i, 1) * PI/180.0d0
            else if (contour_type(1) == 0) then
                virt(i, 1) = virt(i, 1)
            end if
        end do 
    end subroutine Fix_contour_angle

    subroutine Set_communication(this, comm, comm_params)
        class (mesh_base_t)            :: this 
        type(communication_t), pointer            :: comm
        type(communication_parameters_t), pointer :: comm_params

    end subroutine Set_communication

    subroutine Ptr_coordinates_1d (this, dim_num, ptr)
        implicit none

        class (mesh_base_t)               , intent(in out) :: this 
        integer                           , intent(in)     :: dim_num

        real(8), dimension(:,:,:), pointer, intent(out)    :: ptr 

        call this%coordinates%Point_to_data (dim_num, ptr)
    end subroutine Ptr_coordinates_1d

    subroutine Ptr_coordinates_2d (this, ptr_x, ptr_y)
        implicit none
        class (mesh_base_t)               , intent(in out) :: this    
        real(8), dimension(:,:,:), pointer, intent(out)    :: ptr_x   
        real(8), dimension(:,:,:), pointer, intent(out)    :: ptr_y   

        call this%coordinates%Point_to_data (ptr_x, ptr_y)
    end subroutine Ptr_coordinates_2d

    subroutine Ptr_coordinates_3d (this, ptr_x, ptr_y, ptr_z)
        implicit none

        class (mesh_base_t)               , intent(in out) :: this    
        real(8), dimension(:,:,:), pointer, intent(out)    :: ptr_x   
        real(8), dimension(:,:,:), pointer, intent(out)    :: ptr_y   
        real(8), dimension(:,:,:), pointer, intent(out)    :: ptr_z   

        call this%coordinates%Point_to_data (ptr_x, ptr_y, ptr_z)
    end subroutine Ptr_coordinates_3d

    subroutine Point_to_r_factor(this, ptr)
        implicit none
        class (mesh_base_t)               , intent(in)  :: this  
        real(8), dimension(:,:,:), pointer, intent(out) :: ptr   

        call this%r_factor%Point_to_data (ptr)
    end subroutine Point_to_r_factor

    subroutine Clean_mesh (this)
        implicit none
        class (mesh_base_t), intent(in out) :: this 
        integer                             :: i

        deallocate (this% coordinates)
        deallocate (this% start_layer_index_r)
        deallocate (this% start_layer_index_t)
    end subroutine Clean_mesh


    subroutine Calculate_average_coordinates(this)
        implicit none
        class (mesh_base_t), intent(in out) :: this 
    end subroutine Calculate_average_coordinates


    subroutine Exchange_virtual_space_blocking(this, ghost_width)
        class (mesh_base_t), intent(in out) :: this  
        integer, optional :: ghost_width
        integer :: ghost_width_1
        integer :: i
        if (.not. present(ghost_width)) then
            ghost_width_1 = 1   
        else
            ghost_width_1 = ghost_width
        end if

        call this%coordinates%Exchange_virtual_space_blocking(ghost_width_1)

    end subroutine Exchange_virtual_space_blocking
   
    subroutine Write_mesh_abstract(this, unit, iostat, iomsg)
        class (mesh_base_t), intent(in) :: this  
        integer,      intent(in)     :: unit
        integer,      intent(out)    :: iostat
        character(*), intent(in out)  :: iomsg
    end subroutine Write_mesh_abstract

    subroutine Write_mesh_base(this, unit, iostat, iomsg)
        class (mesh_base_t), intent(in) :: this  
        integer,      intent(in)     :: unit
        integer,      intent(out)    :: iostat
        character(*), intent(in out) :: iomsg

#ifdef DEBUG
        write(*,*) "@@@ in Write_mesh_base @@@"
#endif

        call this%coordinates%Write_quantity_abstract(unit,iostat=iostat, iomsg=iomsg)

        write(unit, iostat=iostat, iomsg=iomsg) &
            shape(this%start_layer_index_r),&
            shape(this%start_layer_index_t), &
            size(this%n_layers_r), &
            size(this%n_layers_t)

        write(unit, iostat=iostat, iomsg=iomsg) &
            this%nxp ,&
            this%nyp, &
            this%nzp, &
            this%dimension, &
            this%r_factor, &
            this%cyl, &
            this%omcyl, &
            this%start_layer_index_r, &
            this%start_layer_index_t, &
            this%n_layers_r, &
            this%n_layers_t, &
            this%sphere_factor, &
            this%mesh_type

#ifdef DEBUG
        write(*,*) &
            'shape_start_layer_index_r', &
            shape(this%start_layer_index_r), &
            'shape_start_layer_index_t', &
            shape(this%start_layer_index_t), &
            'size_n_layers_r', &
            size(this%n_layers_r), &
            'size_n_layers_t', &
            size(this%n_layers_t)

        write(*,*) &
            'nxp', &
            this%nxp, &
            'nyp', &
            this%nyp, &
            'nzp', &
            this%nzp, &
            'dimension', &
            this%dimension, &
            'cyl', &
            this%cyl, &
            'omcyl', &
            this%omcyl, &
            'start_layer_index_r', &
            this%start_layer_index_r, &
            'start_layer_index_t', &
            this%start_layer_index_t, &
            'n_layers_r', &
            this%n_layers_r, &
            'n_layers_t', &
            this%n_layers_t, &
            'sphere_factor', &
            this%sphere_factor, &
            'mesh_type', &
            this%mesh_type, &
            '###'

        write(*,*) "@@@ end Write_mesh_base @@@"
#endif

    end subroutine Write_mesh_base

    subroutine Read_mesh_abstract(this, unit, iostat, iomsg)
        class (mesh_base_t), intent(in out) :: this  
        integer,      intent(in)     :: unit
        integer,      intent(out)    :: iostat
        character(*), intent(in out)  :: iomsg
    end subroutine Read_mesh_abstract

    subroutine Read_mesh_base(this, unit, iostat, iomsg)
        class (mesh_base_t), intent(in out) :: this  
        integer,      intent(in)     :: unit
        integer,      intent(out)    :: iostat
        character(*), intent(in out)  :: iomsg
        integer, dimension(2) :: shape_start_layer_index_r
        integer, dimension(2) :: shape_start_layer_index_t
        integer :: size_n_layers_r
        integer :: size_n_layers_t


#ifdef DEBUG
        write(*,*) "@@@ in Read_mesh_base @@@"
#endif
        call this%coordinates%Read_quantity_abstract(unit,iostat=iostat, iomsg=iomsg)

        read(unit, iostat=iostat, iomsg=iomsg) &
            shape_start_layer_index_r, &
            shape_start_layer_index_t, &
            size_n_layers_r, &
            size_n_layers_t


        deallocate(this%start_layer_index_r)
        allocate(this%start_layer_index_r( &
            shape_start_layer_index_r(1), &
            shape_start_layer_index_r(2)))
        deallocate(this%start_layer_index_t)
        allocate(this%start_layer_index_t( &
            shape_start_layer_index_t(1), &
            shape_start_layer_index_t(2)))
        deallocate(this%n_layers_r)
        allocate(this%n_layers_r(size_n_layers_r))
        deallocate(this%n_layers_t)
        allocate(this%n_layers_t(size_n_layers_t))

        read(unit, iostat=iostat, iomsg=iomsg) &
            this%nxp, &
            this%nyp, &
            this%nzp, &
            this%dimension, &
            this%r_factor, &
            this%cyl, &
            this%omcyl, &
            this%start_layer_index_r, &
            this%start_layer_index_t, &
            this%n_layers_r, &
            this%n_layers_t, &
            this%sphere_factor, &
            this%mesh_type

#ifdef DEBUG
        write(*,*) &
            'shape_start_layer_index_r', &
            shape(this%start_layer_index_r), &
            'shape_start_layer_index_t', &
            shape(this%start_layer_index_t), &
            'size_n_layers_r', &
            size(this%n_layers_r), &
            'size_n_layers_t', &
            size(this%n_layers_t)

        write(*,*) &
            'nxp', &
            this%nxp, &
            'nyp', &
            this%nyp, &
            'nzp', &
            this%nzp, &
            'dimension', &
            this%dimension, &
            'cyl', &
            this%cyl, &
            'omcyl', &
            this%omcyl, &
            'start_layer_index_r', &
            this%start_layer_index_r, &
            'start_layer_index_t', &
            this%start_layer_index_t, &
            'n_layers_r', &
            this%n_layers_r, &
            'n_layers_t', &
            this%n_layers_t, &
            'sphere_factor', &
            this%sphere_factor, &
            'mesh_type', &
            this%mesh_type, &
            '###'

        write(*,*) "@@@ end Read_mesh_base @@@"
#endif

    end subroutine Read_mesh_base
   
end module mesh_base_module
