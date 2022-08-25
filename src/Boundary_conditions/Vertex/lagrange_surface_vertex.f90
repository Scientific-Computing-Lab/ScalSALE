
module lagrange_surface_vertex_module
    use vertex_boundary_condition_module, only : vertex_boundary_condition_t
    use data_module, only : data_t
    use quantity_module, only : quantity_t
    implicit none
    private

    type, extends(vertex_boundary_condition_t), public :: lagrange_surface_vertex_t
        private
    contains

        procedure, public :: Calculate => lagrange_surface_vertex_calculate
        procedure, public :: Calculate_face
        procedure, public :: Calculate_corner
    end type lagrange_surface_vertex_t



contains
    subroutine lagrange_surface_vertex_calculate (this, v_quantity, coordinates, angle, edge_num)
        class (lagrange_surface_vertex_t)           , intent (in out) :: this      
        class(quantity_t), intent (in out) :: v_quantity  
        type(data_t), dimension(:), pointer,intent(inout)        :: coordinates
        real (8)                           , intent (in)     :: angle    
        integer                            , intent (in)     :: edge_num 


    end subroutine lagrange_surface_vertex_calculate

    subroutine Calculate_face (this, v_quantity, angle, edge_num)
        class (lagrange_surface_vertex_t), intent (in out)     :: this      

        class(quantity_t), intent (in out) :: v_quantity  
        real (8)                           , intent (in)     :: angle    
        integer                            , intent (in)     :: edge_num 

    end subroutine

    subroutine Calculate_corner (this, v_quantity, angle, edge_num)
        implicit none
        class (lagrange_surface_vertex_t), intent (in out)     :: this      

        class(quantity_t), intent (in out) :: v_quantity  
        real (8)                           , intent (in)     :: angle    
        integer                            , intent (in)     :: edge_num 

    end subroutine
end module lagrange_surface_vertex_module
