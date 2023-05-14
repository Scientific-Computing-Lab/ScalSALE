
module cell_boundary_condition_module
    use data_module    , only : data_t
    use quantity_module, only : quantity_t
    implicit none
    private
    public :: cell_boundary_condition_t, cell_bc_wrapper_t
    type, abstract :: cell_boundary_condition_t

    contains


        procedure(Calculate_interface), deferred :: Calculate
    end type

    abstract interface
        subroutine Calculate_interface (this, c_quantity, edge_num)
             import :: cell_boundary_condition_t, data_t, quantity_t

            class (cell_boundary_condition_t) , intent (in out)     :: this      
            class(quantity_t)   , intent (in out)     :: c_quantity  
            integer             , intent (in)         :: edge_num  

        end subroutine

    end interface

    type cell_bc_wrapper_t
        class(cell_boundary_condition_t), pointer :: bc
    end type cell_bc_wrapper_t

contains

end module cell_boundary_condition_module
