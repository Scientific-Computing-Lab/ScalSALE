
module vertex_boundary_condition_module
    use data_module, only : data_t
    use quantity_module, only : quantity_t

    implicit none
    private
    public :: vertex_boundary_condition_t, vertex_bc_wrapper_t

    type, abstract :: vertex_boundary_condition_t

        private
    contains




        procedure(Calculate_interface)       , deferred :: Calculate

        procedure(Calculate_face_interface)  , deferred :: Calculate_face

        procedure(Calculate_corner_interface), deferred :: Calculate_corner

    end type

    abstract interface

        subroutine Calculate_interface (this, v_quantity, coordinates, angle, edge_num)
            import :: vertex_boundary_condition_t, data_t, quantity_t

            class (vertex_boundary_condition_t), intent (in out)     :: this      

            class(quantity_t), intent (in out) :: v_quantity  
            type(data_t), dimension(:), pointer,intent(inout)        :: coordinates
            real (8)                           , intent (in)     :: angle    
            integer                            , intent (in)     :: edge_num 

        end subroutine

        subroutine Calculate_face_interface (this, v_quantity, angle, edge_num)
            import :: vertex_boundary_condition_t, data_t, quantity_t

            class (vertex_boundary_condition_t), intent (in out)     :: this      

            class(quantity_t), intent (in out) :: v_quantity  
            real (8)                           , intent (in)     :: angle    
            integer                            , intent (in)     :: edge_num 

        end subroutine

        subroutine Calculate_corner_interface (this, v_quantity, angle, edge_num)
            import :: vertex_boundary_condition_t, data_t, quantity_t

            class (vertex_boundary_condition_t), intent (in out)     :: this      

            class(quantity_t), intent (in out) :: v_quantity  
            real (8)                           , intent (in)     :: angle    
            integer                            , intent (in)     :: edge_num 

        end subroutine

    end interface

    type vertex_bc_wrapper_t
        class(vertex_boundary_condition_t), pointer :: bc
    end type vertex_bc_wrapper_t


contains


end module vertex_boundary_condition_module
