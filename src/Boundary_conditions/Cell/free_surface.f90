
module free_surface_module
    use cell_boundary_condition_module, only : cell_boundary_condition_t
    use data_module    , only : data_t
    use quantity_module, only : quantity_t

    implicit none
    private

    type, extends(cell_boundary_condition_t), public :: free_surface_t
        private
    contains

        procedure, public :: Calculate => Free_surface_calculate
    end type free_surface_t

contains
    subroutine Free_surface_calculate (this, c_quantity, edge_num)
        class (free_surface_t) , intent (in out)     :: this      
        class(quantity_t)   , intent (in out)     :: c_quantity  
        integer             , intent (in)         :: edge_num  
    end subroutine

end module free_surface_module
