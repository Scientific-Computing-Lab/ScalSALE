
module diagnostics_group_module
   use hydro_step_module, only : hydro_step_t
   use diagnostic_module, only : diagnostic_t, diagnostic_wrapper_t
   implicit none
   private
   public :: diagnostics_group_t

   type :: diagnostics_group_t
      private

      type (diagnostic_wrapper_t) , dimension(:), allocatable :: diagnostics

      class (hydro_step_t), pointer :: hydro_step              
      integer                       :: diagnostics_group_file  
      integer                       :: group_size              
      integer                       :: start_time              
      real(8)                       :: time_step_size          
      integer                       :: diagnostics_counter = 0 

   contains


      procedure :: Apply

      procedure :: Clean_diagnostics_group

      procedure :: Add_diagnostic


   end type

   interface diagnostics_group_t
      module procedure Constructor
   end interface diagnostics_group_t


contains

   type(diagnostics_group_t) function Constructor (hydro, diag_group_file, group_s, start_t, time_step_s)
      class(hydro_step_t), pointer, intent(in)     :: hydro           
      integer                     , intent(in)     :: diag_group_file 
      integer                     , intent(in)     :: group_s         
      integer                     , intent(in)     :: start_t         
      real(8)                     , intent(in)     :: time_step_s     
      integer                                      :: i

      allocate (diagnostic_wrapper_t :: Constructor%diagnostics (group_s))

      Constructor%hydro_step => hydro
      Constructor%diagnostics_group_file = diag_group_file
      Constructor%group_size = group_s
      Constructor%start_time = start_t
      Constructor%time_step_size = time_step_s
   end function

   subroutine Apply(this)
      class(diagnostics_group_t) :: this 

      integer :: i

      do i=1, this%group_size
         call this%diagnostics(i)%diag%Apply()
      end do
   end subroutine

   subroutine Clean_diagnostics_group (this)
      class (diagnostics_group_t), intent(in out) :: this 

      deallocate (this%diagnostics)
   end subroutine Clean_diagnostics_group


   subroutine Add_diagnostic (this, new_diagnostic)
      class (diagnostics_group_t) , intent(in out) :: this           
      class (diagnostic_t), target, intent(in)     :: new_diagnostic 

      type (diagnostic_wrapper_t)                  :: new_wrapper

      new_wrapper%diag => new_diagnostic
      this%diagnostics(this%diagnostics_counter + 1) = new_wrapper
      this%diagnostics_counter = this%diagnostics_counter + 1
   end subroutine Add_diagnostic



end module diagnostics_group_module

