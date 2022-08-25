
module material_base_module
   use energy_module                 , only : energy_t
   use cell_mass_module              , only : cell_mass_t
   use vof_module                    , only : vof_t
   use cell_boundary_condition_module, only : cell_bc_wrapper_t
   use data_module                   , only : data_t
   use communication_module, only : communication_t
   use communication_parameters_module, only : communication_parameters_t
   use material_quantity_module      , only : material_quantity_t
   implicit none
   private

   type, abstract, public :: material_base_t
      private

      type (material_quantity_t), public          , pointer :: vof
      type (material_quantity_t),public        , pointer :: sie
      type (material_quantity_t),pointer, public      :: cell_mass
      type (material_quantity_t)                , pointer  :: initial_layers_of_mats

      integer,dimension(:), pointer                        :: material_ids
      integer,dimension(:), pointer,public                 :: nrg_calc
      integer, public :: nmats

   contains

      procedure, public :: Init_material_base

      procedure, public :: Set_communication_material_base

      procedure, public :: Clean_material_base

      procedure, public :: Point_to_initial_layers

      procedure, public :: Write_material_abstract
      procedure, public :: Write_material_base
      generic :: write(unformatted) => Write_material_base

      procedure, public :: Read_material_abstract
      procedure, public :: Read_material_base
      generic :: read(unformatted) => Read_material_base

   end type material_base_t


contains

   subroutine Init_material_base(this, nxp, nyp, nzp,nmats, mat_ids, bc_cell, bc_params)
      use boundary_parameters_module, only : boundary_parameters_t
      implicit none
      class(material_base_t)                , intent(in out)       :: this
      integer,dimension(:), allocatable         , intent(in)           :: mat_ids
      integer                               , intent(in)           :: nxp           
      integer                               , intent(in)           :: nyp           
      integer                               , intent(in)           :: nzp           
      integer                               , intent(in)           :: nmats
!      real(8)                               , intent(in)           :: sie_0
      type(cell_bc_wrapper_t), dimension(:), pointer,  intent(in) :: bc_cell   
      type(boundary_parameters_t), pointer, intent(in) :: bc_params
integer :: i
      allocate(this%cell_mass)
      allocate(this%initial_layers_of_mats)
      allocate(this%sie)
      allocate(this%vof)
      allocate(this%material_ids(nmats))
      this%nmats = nmats

do i = 1, nmats
      this%material_ids(i) = mat_ids(i)
end do

      this%vof = material_quantity_t(0d0, nxp, nyp, nzp,nmats, bc_cell, bc_params)
      this%sie = material_quantity_t (0d0, nxp, nyp, nzp, nmats,bc_cell, bc_params)
      this%cell_mass = material_quantity_t (0d0, nxp, nyp, nzp,nmats, bc_cell, bc_params)
      this%initial_layers_of_mats = material_quantity_t (0d0, nxp, nyp, nzp, nmats)

!      if (sie_0(1) == 0) then
!         this%nrg_calc = 1
!      else
!         this%nrg_calc = 0
!      end if



   end subroutine



   subroutine Point_to_initial_layers(this, ptr)
      implicit none
      class (material_base_t)                  , intent(in out) :: this  
      real(8)       , dimension(:,:,:,:), pointer, intent(out)    :: ptr

      call this%initial_layers_of_mats%Point_to_data (ptr)

   end subroutine Point_to_initial_layers



   subroutine Set_communication_material_base(this, comm, comm_params)
      class (material_base_t)            :: this 
      type(communication_t), pointer            :: comm
      type(communication_parameters_t), pointer :: comm_params
      call this%vof%Set_communication(comm, comm_params)

      call this% sie %Set_communication(comm, comm_params)
      call this%cell_mass%Set_communication(comm, comm_params)
      call this%initial_layers_of_mats%Set_communication(comm, comm_params)

   end subroutine Set_communication_material_base




   subroutine Clean_material_base(this)
      class (material_base_t), intent(in out) :: this 

!      call this%vof            %Clean_vof
!      call this%sie            %Clean_energy
!      call this%cell_mass      %Clean_cell_mass
   end subroutine Clean_material_base

   subroutine Write_material_abstract(this, unit, iostat, iomsg)
      class (material_base_t), intent(in) :: this  
      integer,      intent(in)     :: unit
      integer,      intent(out)    :: iostat
      character(*), intent(in out)  :: iomsg
   end subroutine Write_material_abstract

   subroutine Write_material_base(this, unit, iostat, iomsg)
      class (material_base_t), intent(in) :: this  
      integer,      intent(in)     :: unit
      integer,      intent(out)    :: iostat
      character(*), intent(in out) :: iomsg

!#ifdef DEBUG
!      write(*,*) "@@@ in Write_material_base @@@"
!#endif
!
!!      call this%vof%Write_quantity_abstract(unit, iostat=iostat, iomsg=iomsg)
!!      call this%sie%Write_quantity_abstract(unit, iostat=iostat, iomsg=iomsg)
!!      call this%cell_mass%Write_quantity_abstract(unit, iostat=iostat, iomsg=iomsg)
!
!      write(unit, iostat=iostat, iomsg=iomsg) &
!         this%initial_layers_of_mats, &
!         this%material_id
!
!#ifdef DEBUG
!      write(*,*) &
!         'material_id', &
!         this%material_id, &
!         '###'
!
!      write(*,*) "@@@ end Write_material_base @@@"
!#endif

   end subroutine Write_material_base

   subroutine Read_material_abstract(this, unit, iostat, iomsg)
      class (material_base_t), intent(in out) :: this  
      integer,      intent(in)     :: unit
      integer,      intent(out)    :: iostat
      character(*), intent(in out)  :: iomsg
   end subroutine Read_material_abstract

   subroutine Read_material_base(this, unit, iostat, iomsg)
      class (material_base_t), intent(in out) :: this  
      integer,      intent(in)     :: unit
      integer,      intent(out)    :: iostat
      character(*), intent(in out)  :: iomsg
!
!#ifdef DEBUG
!      write(*,*) "@@@ in Read_material_base @@@"
!#endif
!
!      call this%vof%Read_quantity_abstract(unit, iostat=iostat, iomsg=iomsg)
!      call this%sie%Read_quantity_abstract(unit, iostat=iostat, iomsg=iomsg)
!      call this%cell_mass%Read_quantity_abstract(unit, iostat=iostat, iomsg=iomsg)
!
!      read(unit, iostat=iostat, iomsg=iomsg) &
!         this%initial_layers_of_mats, &
!         this%material_id
!
!#ifdef DEBUG
!      write(*,*) &
!         'material_id', &
!         this%material_id, &
!         '###'
!
!      write(*,*) "@@@ end Read_material_base @@@"
!#endif

   end subroutine Read_material_base

end module material_base_module

