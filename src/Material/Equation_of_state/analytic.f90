
module analytic_module
   use equation_of_state_module, only : equation_of_state_t
   use data_module, only : data_t
   implicit none
   private

   type, extends(equation_of_state_t), public :: analytic_t
      private
   contains

      procedure, public :: Calculate => Calculate_analytic
   end type analytic_t



contains

   subroutine Calculate_analytic (this, pressure, sound_vel, density, sie,&
                                      dp_de, dp_drho, dt_de, dt_drho,gamma_gas, atomic_mass, temperature, mat_num, nrg_or_tmp,&
                                      nx, ny, nz, vof, emf)

         class (analytic_t)        , intent (inout) :: this        
         real(8), dimension (:,:,:), pointer, intent (inout) :: pressure    
         real(8), dimension (:,:,:), pointer, intent (inout) :: sound_vel   
         real(8), dimension (:,:,:), pointer, intent (inout) :: density     
         real(8), dimension (:,:,:), pointer, intent (inout) :: sie         
         real(8), dimension (:,:,:), pointer, intent (inout) :: temperature 
         real(8), dimension (:,:,:), pointer, intent (inout) :: dp_de       
         real(8), dimension (:,:,:), pointer, intent (inout) :: dp_drho     
         real(8), dimension (:,:,:), pointer, intent (inout) :: dt_de       
         real(8), dimension (:,:,:), pointer, intent (inout) :: dt_drho     
         real(8)                            , intent (in)    :: atomic_mass 
         real(8)                            , intent (in)    :: gamma_gas   
         integer                            , intent (in   ) :: mat_num     

         integer                            , intent (in   ) :: nrg_or_tmp
         integer                            , intent (in   ) :: nx
         integer                            , intent (in   ) :: ny
         integer                            , intent (in   ) :: nz
         real(8), dimension (:,:,:), pointer, intent (in   ) :: vof
         real(8)                            , intent (in   ) :: emf



      end subroutine Calculate_analytic

end module analytic_module
