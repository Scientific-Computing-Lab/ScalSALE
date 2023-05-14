module cr_module
   use hydro_step_module, only : hydro_step_t
   use boundary_parameters_module, only : boundary_parameters_t
   use time_module, only : time_t
   use, intrinsic :: iso_c_binding

   implicit none
   include 'scrf.h'
   private

   type, public :: cr_t
      type (hydro_step_t), pointer :: hydro  
      type (time_t), pointer       :: time   
      type(boundary_parameters_t) , pointer :: bc_params
      character(len=1), dimension(:),allocatable            :: run_name
      integer                                               :: with_cr
      real(8)                      :: finish_time   

      contains
      procedure, public :: Get_ckpt_name
      procedure, public  :: Checkpoint
      procedure, public  :: Restart
      procedure, public  :: Clean_cr
      procedure, private :: Write_ckpt_to_file
      procedure, private :: Read_ckpt_from_file
      procedure, private :: Overwrite_persistent_quantities
   end type cr_t


   interface cr_t
      module procedure Constructor
   end interface cr_t

contains

   type(cr_t) function Constructor(hydro, time, bc_params, run_name, with_cr)
      type (hydro_step_t), pointer, intent(inout)                :: hydro  
      type (time_t), pointer, intent(inout)                      :: time   
      type (boundary_parameters_t), pointer, intent(inout)       :: bc_params   
      character(len=1), dimension(:),allocatable, intent(inout)  :: run_name
      integer, intent(in)                                        :: with_cr
#ifndef TEST
      integer :: ierr, i
      character(len=size(run_name)) :: tmp


      Constructor%with_cr = with_cr
      if (with_cr == 1) then
         Constructor%hydro => hydro
         Constructor%time => time
         Constructor%bc_params => bc_params
         Constructor%run_name = run_name
         Constructor%finish_time = time%finish_time
         call scr_init(ierr)
      end if

#endif
   end function Constructor


    subroutine Get_ckpt_name(this, name, ckpt_name)
      class (cr_t),                                   intent(in out) :: this
      character(len=1), dimension(:),allocatable,     intent(in)     :: name
      character(len=size(this%run_name)+len(name)+1), intent(in out) :: ckpt_name
      integer                                                        :: i, j

      if (this%with_cr == 1) then

         do i=1,size(this%run_name)
            ckpt_name(i:i) = this%run_name(i)
         end do
         ckpt_name(i:i) = '/'
         do j=1,size(name)
            ckpt_name(i+j:i+j) = name(j)
         end do
      end if
   end subroutine Get_ckpt_name


   subroutine Write_ckpt_to_file(this, fname, valid)
      class (cr_t),     intent(in out) :: this
      character(len=*), intent(in)     :: fname
      integer,          intent(out)    :: valid
      integer                          :: fd

#ifdef DEBUG
      write(*,*) 'Writing checkpoint: ', trim(fname)
#endif
      if (this%with_cr == 1) then
         open(newunit=fd, file=fname, status='new', form='unformatted')
         write(fd) &
            this%time, &
            this%hydro
         call flush(fd)
         close(fd)
         valid = 1
      end if
   end subroutine Write_ckpt_to_file


   subroutine Checkpoint(this, name)
      class (cr_t),     intent(in out) :: this
      character(len=*), intent(in)     :: name
#ifndef TEST
      integer                          :: ierr, need_checkpoint, valid, should_exit
      character(len=SCR_MAX_FILENAME)  :: scr_name
      integer :: fd

      if (this%with_cr == 1) then
         call SCR_Need_checkpoint(need_checkpoint, ierr)
         if (need_checkpoint == 1) then
            call SCR_Start_checkpoint(ierr)
            call SCR_Route_file(name , scr_name, ierr)
            write (*,*) '@@ Writing checkpoint to ', trim(scr_name)
            call this%Write_ckpt_to_file(scr_name, valid)
            call SCR_Complete_checkpoint(valid, ierr)
         end if

         call SCR_Should_exit(should_exit, ierr)
         if (should_exit == 1) then
         end if
      end if
#endif
   end subroutine Checkpoint


   subroutine Read_ckpt_from_file(this, fname, valid)
      class (cr_t),     intent(in out) :: this
      character(len=*), intent(in)     :: fname
      integer,          intent(out)    :: valid
      integer                          :: fd
      character(len=25)                :: content

#ifdef DEBUG
      write(*,*) 'Reading checkpoint: ', trim(fname)
#endif
write(*,*) "Reading checkpoint:", trim(fname)
      open(unit=fd, file=fname, status='old', form='unformatted')
      read(fd) &
         this%time, &
         this%hydro
      close(fd)

      valid = 1

   end subroutine Read_ckpt_from_file


   subroutine Overwrite_persistent_quantities(this)
      class (cr_t), intent(in out)                    :: this

      write (*,*) 'Overwriting persistent parameters'
      write (*,*) 'finish_time: ', this%time%finish_time,' -> ',this%finish_time
      this%time%finish_time = this%finish_time
   end subroutine Overwrite_persistent_quantities


   subroutine Restart(this, name)
      class (cr_t), intent(in out)                    :: this
      character(len=*), intent(in)                    :: name
#ifndef TEST
      integer                                         :: ierr, have_restart, valid
      character(len=SCR_MAX_FILENAME)                 :: ds_name
      character(len=SCR_MAX_FILENAME)                 :: scr_name

      if (this%with_cr == 1) then

         call SCR_Have_restart(have_restart, ds_name, ierr)
         if (have_restart == 1) then
            call SCR_Start_restart(ds_name, ierr)
            call SCR_Route_file(name, scr_name, ierr)
            write (*,*) '@@ Restarting from ', trim(scr_name)
            call this%Read_ckpt_from_file(scr_name, valid)
            call this%Overwrite_persistent_quantities()
            call SCR_Complete_restart(valid, ierr)
         end if
      end if
#endif
   end subroutine Restart


   subroutine Clean_cr(this)
     class (cr_t), intent(in out) :: this   
#ifndef TEST
     integer                      :: ierr

     if (this%with_cr == 1) call scr_finalize(ierr)
#endif
   end subroutine Clean_cr

end module
