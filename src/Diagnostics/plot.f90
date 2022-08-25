
module plot_module
   use diagnostic_module, only : diagnostic_t
   implicit none
   private
   public :: plot

   type, extends(diagnostic_t) :: plot
      private

   contains


      procedure :: apply => apply_plot


   end type



contains

   subroutine apply_plot(this)
      class (plot), intent(in out) :: this  


   end subroutine



end module plot_module
