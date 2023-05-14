
module lagrange_surface_cell_module
    use cell_boundary_condition_module, only : cell_boundary_condition_t
    use data_module    , only : data_t
    use quantity_module, only : quantity_t
    implicit none
    private

    type, extends(cell_boundary_condition_t), public :: lagrange_surface_cell_t
        private
    contains

        procedure, public :: Calculate => lagrange_surface_cell_calculate
    end type lagrange_surface_cell_t

contains
    subroutine lagrange_surface_cell_calculate (this, c_quantity, edge_num)
        class (lagrange_surface_cell_t) , intent (in out)     :: this      
        class(quantity_t)   , intent (in out)     :: c_quantity  
        integer             , intent (in)         :: edge_num  
    end subroutine

end module lagrange_surface_cell_module
