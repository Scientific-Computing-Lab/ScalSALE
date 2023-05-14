
module quantity_module
    use data_module, only : data_t, Get_copy
    use communication_module, only : communication_t
    use communication_parameters_module  , only : communication_parameters_t
    use boundary_parameters_module       , only : boundary_parameters_t
    use parallel_parameters_module       , only : parallel_parameters_t
    use data_4d_module, only : data_4d_t
    implicit none
    private
    public :: quantity_t

    type, abstract :: quantity_t
        private

        type (data_t) , dimension(:), pointer, public :: data
        type (data_4d_t), pointer, public :: data_4d
        type(boundary_parameters_t), pointer, public :: boundary_params
        type(parallel_parameters_t), pointer, public :: parallel_params

        integer                     , public              :: d1    
        integer                     , public              :: d2    
        integer                     , public              :: d3    
        integer                     , public              :: number_of_axises  
    contains


        procedure, public  :: Init_quantity_init_arr
        procedure, public  :: Init_quantity_no_bc

        procedure, public  :: Init_quantity_init_val

        procedure, public  :: Init_quantity_no_init



        procedure, private  :: Ptr_coordinates_1d

        procedure, private  :: Ptr_coordinates_1d_first

        procedure, private  :: Ptr_coordinates_2d

        procedure, private  :: Ptr_coordinates_3d

        procedure, private  :: Ptr_coordinates_4d


        generic,   public    :: Point_to_data  =>         &
            Ptr_coordinates_1d,       &
            Ptr_coordinates_1d_first, &
            Ptr_coordinates_2d,       &
            Ptr_coordinates_3d,       &
            Ptr_coordinates_4d

        procedure, public :: Clean_quantity

        procedure, public :: Set_communication

        procedure, public :: debug_check_nan

        procedure, public :: Exchange_virtual_space_blocking

        procedure, public :: Exchange_virtual_space_nonblocking

        procedure, public :: Exchange_end
        procedure, public :: Init_quantity_init_val_4d

        procedure, public :: debug_print
        procedure, public :: Write_quantity_abstract
        procedure, public :: Write_quantity
        generic :: write(unformatted) => Write_quantity

        procedure, public :: Read_quantity_abstract
        procedure, public :: Read_quantity
        generic :: read(unformatted) => Read_quantity


    end type

    public    :: Get_data_copy


contains


    subroutine Init_quantity_init_arr(this, initial_data, d1, d2, d3, axises_num, bc_params)
        implicit none
        class(quantity_t)        , intent(in out) :: this      
        real(8), dimension(:,:,:), intent(in)     :: initial_data  
        integer                  , intent(in)     :: d1            
        integer                  , intent(in)     :: d2            
        integer                  , intent(in)     :: d3            
        type(boundary_parameters_t), pointer, intent(in) :: bc_params

        integer                  , intent(in)     :: axises_num
        integer                                   :: i
        nullify(this%data_4d)
        nullify(this%data)
        allocate (data_t :: this%data (axises_num))
        do i=1, axises_num
            this%data(i) = data_t (initial_data, d1, d2, d3)
        end do
        this%d1 = d1 - 1
        this%d2 = d2 - 1
        this%d3 = d3 - 1
        this%number_of_axises = axises_num
        this%boundary_params => bc_params
        nullify(this%data_4d)

    end subroutine

    subroutine Init_quantity_init_val(this, initial_val, d1, d2, d3, axises_num, bc_params)
        implicit none
        class(quantity_t) , intent(in out) :: this      
        real(8)           , intent(in)     :: initial_val  
        integer           , intent(in)     :: d1            
        integer           , intent(in)     :: d2            
        integer           , intent(in)     :: d3            
        type(boundary_parameters_t), pointer, intent(in) :: bc_params

        integer           , intent(in)     :: axises_num
        integer                            :: i
        nullify(this%data_4d)
        nullify(this%data)
        allocate (data_t :: this%data (axises_num))
        do i=1, axises_num
            this%data(i) = data_t (initial_val, d1, d2, d3)
        end do
        this%d1 = d1 - 1
        this%d2 = d2 - 1
        this%d3 = d3 - 1
        this%number_of_axises = axises_num
        this%boundary_params => bc_params
        nullify(this%data_4d)
    end subroutine

    subroutine Init_quantity_init_val_4d(this, initial_val, d1, d2, d3, d4, axises_num, bc_params)
        implicit none
        class(quantity_t) , intent(in out) :: this
        real(8)           , intent(in)     :: initial_val
        integer           , intent(in)     :: d1
        integer           , intent(in)     :: d2
        integer           , intent(in)     :: d3
        integer           , intent(in)     :: d4
        type(boundary_parameters_t), pointer, intent(in) :: bc_params

        integer           , intent(in)     :: axises_num
        integer                            :: i
        nullify(this%data_4d)
        nullify(this%data)
        allocate (data_4d_t :: this%data_4d)
        this%data_4d = data_4d_t(initial_val, d1,d2,d3,d4)
        this%d1 = d1 - 1
        this%d2 = d2 - 1
        this%d3 = d3 - 1
        this%number_of_axises = axises_num
        this%boundary_params => bc_params
        nullify(this%data)
    end subroutine

    subroutine Init_quantity_no_bc(this, initial_val, d1, d2, d3, d4, axises_num)
        implicit none
        class(quantity_t) , intent(in out) :: this
        real(8)           , intent(in)     :: initial_val
        integer           , intent(in)     :: d1
        integer           , intent(in)     :: d2
        integer           , intent(in)     :: d3
        integer           , intent(in)     :: d4

        integer           , intent(in)     :: axises_num
        integer                            :: i
        nullify(this%data_4d)
        nullify(this%data)
        allocate (data_4d_t :: this%data_4d)
        this%data_4d = data_4d_t(initial_val, d1,d2,d3,d4)
        this%d1 = d1 - 1
        this%d2 = d2 - 1
        this%d3 = d3 - 1
        this%number_of_axises = axises_num
        nullify(this%data)
    end subroutine

    subroutine Init_quantity_no_init(this, d1, d2, d3, axises_num, bc_params)
        implicit none
        class(quantity_t) , intent(in out) :: this      
        integer           , intent(in)     :: d1            
        integer           , intent(in)     :: d2            
        integer           , intent(in)     :: d3            
        type(boundary_parameters_t), pointer, intent(in) :: bc_params

        integer           , intent(in)     :: axises_num
        integer                            :: i
        nullify(this%data_4d)
        nullify(this%data)
        allocate (data_t :: this%data (axises_num))
        do i=1, axises_num
            this%data(i) = data_t (d1, d2, d3)
        end do
        this%d1 = d1 - 1
        this%d2 = d2 - 1
        this%d3 = d3 - 1
        this%number_of_axises = axises_num

        this%boundary_params => bc_params

    end subroutine




    subroutine Clean_quantity (this)
        class (quantity_t), intent(in out) :: this   
        integer                            :: i

        do i=1, this%number_of_axises
            call this%data(i)%Clean_data ()
        end do
        deallocate (this%data)
    end subroutine Clean_quantity

    subroutine Ptr_coordinates_1d (this, axis, ptr)
        class (quantity_t)                , intent(in) :: this  
        real(8), dimension(:,:,:), pointer, intent(out)    :: ptr   
        integer                           , intent(in)     :: axis  

        call this%data(axis)%Point_to_data (ptr)
    end subroutine Ptr_coordinates_1d

    subroutine Ptr_coordinates_1d_first (this, ptr)
        class (quantity_t)                , intent(in) :: this  
        real(8), dimension(:,:,:), pointer, intent(out)    :: ptr   
        call this%data(1)%Point_to_data (ptr)
    end subroutine Ptr_coordinates_1d_first

    function Get_data_copy (this, i)
        class (quantity_t), intent(in)     :: this 
        integer           , intent(in)     :: i
        real(8), dimension(:,:,:), pointer :: Get_data_copy  

        Get_data_copy = Get_copy(this%data(i))
    end function Get_data_copy


    subroutine Ptr_coordinates_2d (this, ptr_x, ptr_y)
        class (quantity_t)                , intent(in out) :: this         
        real(8), dimension(:,:,:), pointer, intent(out)    :: ptr_x, ptr_y 

        call this%data(1)%Point_to_data (ptr_x)
        call this%data(2)%Point_to_data (ptr_y)
    end subroutine Ptr_coordinates_2d

    subroutine Ptr_coordinates_3d (this, ptr_x, ptr_y, ptr_z)
        class (quantity_t)         , intent(in out) :: this                       
        real(8), dimension(:,:,:), pointer, intent(out)    :: ptr_x, ptr_y, ptr_z 

        call this%data(1)%Point_to_data (ptr_x)
        call this%data(2)%Point_to_data (ptr_y)
        call this%data(3)%Point_to_data (ptr_z)
    end subroutine Ptr_coordinates_3d


    subroutine Ptr_coordinates_4d (this, ptr_x)
        class (quantity_t)         , intent(in out) :: this
        real(8), dimension(:,:,:,:), pointer, intent(out)    :: ptr_x

        call this%data_4d%Point_to_data (ptr_x)
    end subroutine Ptr_coordinates_4d


    subroutine Set_communication (this, comm, comm_params)
        implicit none
        class (quantity_t), intent(in out) :: this  
        type(communication_t), pointer            :: comm
        type(communication_parameters_t), pointer :: comm_params
        integer :: i
        if (associated(this%data)) then
            do i=1, this%number_of_axises
                call this%data(i)%Set_communication (comm, comm_params)
            end do
        else
            call this%data_4d%Set_communication(comm, comm_params)
        end if
        this%parallel_params => comm%parallel_params
    end subroutine Set_communication


    subroutine Exchange_virtual_space_blocking (this, ghost_width)
        implicit none
        class (quantity_t), intent(in out) :: this  
        integer, optional , intent(in)     :: ghost_width
        integer                            :: tmp_ghost_width
        integer :: i

        if (.not. present(ghost_width)) then
            tmp_ghost_width = 1
        else
            tmp_ghost_width = ghost_width
        end if
        if (associated(this%data)) then
            do i=1, this%number_of_axises
                call this%data(i)%Exchange_virtual_space_blocking (tmp_ghost_width)
            end do
        else

            call this%data_4d%Exchange_virtual_space_blocking (tmp_ghost_width)
        end if


    end subroutine Exchange_virtual_space_blocking

    subroutine Exchange_virtual_space_nonblocking (this, ghost_width)
        class (quantity_t), intent(in out) :: this  
        integer           , optional       :: ghost_width
        integer                            :: tmp_ghost_width
        integer :: i

        if (.not. present(ghost_width)) then
            tmp_ghost_width = 1
        else
            tmp_ghost_width = ghost_width
        end if


        if (associated(this%data)) then
            do i=1, this%number_of_axises
                call this%data(i)%Exchange_virtual_space_nonblocking (tmp_ghost_width)
            end do
        else
            call this%data_4d%Exchange_virtual_space_nonblocking (tmp_ghost_width)
        end if
    end subroutine Exchange_virtual_space_nonblocking

    subroutine Exchange_end (this)
        class (quantity_t), intent(in out) :: this  
        integer :: i

!        do i=1, this%number_of_axises
!            call this%data(i)%Exchange_end ()
!        end do

        if (associated(this%data)) then
            do i=1, this%number_of_axises
                call this%data(i)%Exchange_end ()
            end do
        else
            call this%data_4d%Exchange_end ()
        end if

    end subroutine Exchange_end


    subroutine debug_check_nan(this, caller)
        implicit none
        class (quantity_t), intent(in out) :: this  
        CHARACTER(*), intent(in):: caller
        integer :: size_d, i
        size_d = size(this%data)
        do i = 1,size_d
            call this%data(i)%debug_check_nan(caller)
        end do
    end subroutine debug_check_nan

    subroutine debug_print(this, caller, axis,flag)
        implicit none
        class (quantity_t), intent(in out) :: this  
        integer, optional :: flag
        integer, optional :: axis
        CHARACTER(*) caller
        integer :: width
        integer :: size_d, i, axis_type
        character(1) :: axis_num

        if (.not. present(flag)) then
            width = 0   
        else
            width = flag   
        end if

        if (.not. present(axis)) then
            size_d = size(this%data)
            do i = 1,size_d
                write(axis_num,'(i1)') i
                call this%data(i)%debug_print(caller // " Axis:" // axis_num, width)
            end do
        else
            write(axis_num,'(i1)')axis
            call this%data(i)%debug_print(caller // " Axis:" // axis_num, width)
        end if
    end subroutine debug_print

    subroutine Write_quantity_abstract(this, unit, iostat, iomsg)
        class (quantity_t), intent(in) :: this  
        integer,      intent(in)     :: unit
        integer,      intent(out)    :: iostat
        character(*), intent(in out)  :: iomsg
    end subroutine Write_quantity_abstract

    subroutine Write_quantity(this, unit, iostat, iomsg)
        class (quantity_t), intent(in) :: this  
        integer,      intent(in)     :: unit
        integer,      intent(out)    :: iostat
        character(*), intent(in out) :: iomsg

#ifdef DEBUG
        write(*,*) "@@@ in Write_quantity @@@"
#endif

        write(unit, iostat=iostat, iomsg=iomsg) &
            this%data, &
            this%d1, &
            this%d2, &
            this%d3, &
            this%number_of_axises

#ifdef DEBUG
        write(*,*) &
            'd1', &
            this%d1, &
            'd2', &
            this%d2, &
            'd3', &
            this%d3, &
            'number_of_axises', &
            this%number_of_axises, &
            '###'

        write(*,*) "@@@ end Write_quantity @@@"
#endif

    end subroutine Write_quantity

    subroutine Read_quantity_abstract(this, unit, iostat, iomsg)
        class (quantity_t), intent(in out) :: this  
        integer,      intent(in)     :: unit
        integer,      intent(out)    :: iostat
        character(*), intent(in out)  :: iomsg
    end subroutine Read_quantity_abstract

    subroutine Read_quantity(this, unit, iostat, iomsg)
        class (quantity_t), intent(in out) :: this  
        integer,      intent(in)     :: unit
        integer,      intent(out)    :: iostat
        character(*), intent(in out)  :: iomsg

#ifdef DEBUG
        write(*,*) "@@@ in Read_quantity @@@"
#endif

        read(unit, iostat=iostat, iomsg=iomsg) &
            this%data, &
            this%d1, &
            this%d2, &
            this%d3, &
            this%number_of_axises

#ifdef DEBUG
        write(*,*) &
            'd1', &
            this%d1, &
            'd2', &
            this%d2, &
            'd3', &
            this%d3, &
            'number_of_axises', &
            this%number_of_axises ,&
            '###'

        write(*,*) "@@@ end Read_quantity @@@"
#endif

    end subroutine Read_quantity
end module quantity_module
