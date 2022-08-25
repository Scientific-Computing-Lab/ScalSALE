
module vertex_quantity_module
    use data_module, only : data_t
    use vertex_boundary_condition_module, only : vertex_boundary_condition_t, vertex_bc_wrapper_t
    use quantity_module, only : quantity_t
    use boundary_parameters_module      , only : boundary_parameters_t

    implicit none
    private
    public :: vertex_quantity_t

    type, extends(quantity_t) :: vertex_quantity_t
        private


        type (vertex_bc_wrapper_t) , dimension(:), pointer, public :: boundary_conditions

    contains


        procedure, public   :: Init_vertex_quantity_init_arr

        procedure, public   :: Init_vertex_quantity_init_val

        procedure, public   :: Init_vertex_quantity_no_init

        procedure, public   :: Clean_vertex_quantity

        procedure, public   :: Apply_boundary

        procedure, public :: Write_vertex_quantity

        procedure, public :: Read_vertex_quantity

    end type


contains

    subroutine Init_vertex_quantity_init_arr(this, initial_data, nx, ny, nz, axises_num, bc, bc_params)
        class (vertex_quantity_t)                , intent(in out) :: this         
        real(8), dimension(:,:,:)                , intent(in)     :: initial_data 
        integer                                  , intent(in)     :: nx           
        integer                                  , intent(in)     :: ny           
        integer                                  , intent(in)     :: nz           
        type(boundary_parameters_t), pointer, intent(in) :: bc_params

        integer                                  , intent(in)     :: axises_num
        type (vertex_bc_wrapper_t) , dimension(:), pointer, intent(in)     :: bc           
        integer                                                   :: i

        this%boundary_conditions => bc


        call this%Init_quantity_init_arr (initial_data, nx, ny, nz, axises_num, bc_params)



    end subroutine

    subroutine Init_vertex_quantity_init_val(this, initial_val, nx, ny, nz, axises_num, bc, bc_params)
        class (vertex_quantity_t)                , intent(in out) :: this         
        real(8)                                  , intent(in)     :: initial_val  
        integer                                  , intent(in)     :: nx           
        integer                                  , intent(in)     :: ny           
        integer                                  , intent(in)     :: nz           
        type(boundary_parameters_t), pointer, intent(in) :: bc_params

        integer                                  , intent(in)     :: axises_num

        type (vertex_bc_wrapper_t) , dimension(:), pointer, intent(in)     :: bc           

        integer                                                   :: i
        this%boundary_conditions => bc

        call this%Init_quantity_init_val (initial_val, nx, ny, nz, axises_num, bc_params)

    end subroutine

    subroutine Init_vertex_quantity_no_init(this, nx, ny, nz, axises_num, bc, bc_params)
        class (vertex_quantity_t)                , intent(in out) :: this         
        integer                                  , intent(in)     :: nx           
        integer                                  , intent(in)     :: ny           
        integer                                  , intent(in)     :: nz           
        type(boundary_parameters_t), pointer, intent(in) :: bc_params

        integer                                  , intent(in)     :: axises_num

        type (vertex_bc_wrapper_t) , dimension(:), pointer, intent(in)     :: bc           

        integer                                                   :: i

        this%boundary_conditions => bc


        call this%Init_quantity_no_init (nx, ny, nz, axises_num, bc_params)
    end subroutine


    subroutine Clean_vertex_quantity (this)
        class (vertex_quantity_t), intent(in out) :: this 

        call this%Clean_quantity ()
    end subroutine Clean_vertex_quantity


    subroutine Apply_boundary(this, coordinates_values, is_blocking)
        implicit none
        class (vertex_quantity_t)        , intent(in out) :: this
        type(data_t), dimension(:), pointer,intent(inout)        :: coordinates_values
        logical, optional :: is_blocking 
        integer, dimension(:)            , pointer        :: ed_num
        integer, dimension(:)            , allocatable        :: ed_num_test
        integer                                           :: i,  edge
        logical :: is_blocking_local

        if (.not. present(is_blocking)) then
            is_blocking_local = .True.  
        else
            is_blocking_local = is_blocking   
        end if

        call this%boundary_params%Point_to_edges (ed_num)
        do i = 1, size(ed_num)
            edge = ed_num(i)

            call this%boundary_conditions(edge)%bc% Calculate (this, coordinates_values, 0d0, edge)
        end do

        do i = 1, size(ed_num)
            edge = ed_num(i)
            call this%boundary_conditions(edge)%bc% Calculate_face (this, 0d0, edge)
        end do

        do i = 1, size(ed_num)
            edge = ed_num(i)
            call this%boundary_conditions(edge)%bc% Calculate_corner (this, 0d0, edge)
        end do

        if (is_blocking_local) then
            call this%Exchange_virtual_space_blocking()
        else
            call this%Exchange_virtual_space_nonblocking()
        end if
    end subroutine Apply_boundary

    subroutine Write_vertex_quantity(this, unit, iostat, iomsg)
        class (vertex_quantity_t), intent(in) :: this  
        integer,      intent(in)     :: unit
        integer,      intent(out)    :: iostat
        character(*), intent(in out)  :: iomsg

#ifdef DEBUG
        write(*,*) '@@@ in Write_vertex_quantity @@@'
#endif

        call this%Write_quantity(unit, iostat=iostat, iomsg=iomsg)


#ifdef DEBUG

        write(*,*) '@@@ end Write_vertex_quantity @@@'
#endif

    end subroutine Write_vertex_quantity

    subroutine Read_vertex_quantity(this, unit, iostat, iomsg)
        class (vertex_quantity_t), intent(in out) :: this  
        integer,      intent(in)     :: unit
        integer,      intent(out)    :: iostat
        character(*), intent(in out)  :: iomsg

#ifdef DEBUG
        write(*,*) '@@@ in Read_vertex_quantity @@@'
#endif

        call this%Read_quantity(unit, iostat=iostat, iomsg=iomsg)


#ifdef DEBUG

        write(*,*) '@@@ end Read_vertex_quantity @@@'
#endif

    end subroutine Read_vertex_quantity

end module vertex_quantity_module
