
module constants_module
   implicit none

   real (8) :: PI, AVOGADRO, K_BOLTZMAN, SIGMA, SIGMA_BOLTZMAN, &
               C_LIGHT, A_RADIATION, MC2, TEV, RADIUS_ELECTRON, &
               THIRD, TWIFTH
   integer :: N_SEGMENT_MAX


parameter( PI=acos(-1d0) )
   parameter( AVOGADRO = 6.0221367D23 )                        
   parameter( K_BOLTZMAN = 1.3807D-16 )                        
   parameter( SIGMA_BOLTZMAN = 5.6703D-5 )                     
   parameter( C_LIGHT = 2.99792458D10 )                        
   parameter( A_RADIATION = 4D0 * SIGMA_BOLTZMAN / C_LIGHT )   
   parameter( MC2 = 5.11D5 )                                   
   parameter( TEV = 11605D0 )                                  
   parameter( RADIUS_ELECTRON = 2.818D-13 )                    
   parameter( THIRD = 1.D0/3.D0 )                              
   parameter( TWIFTH = .25D0 * THIRD)                          

   parameter( N_SEGMENT_MAX = 12)                               

end module constants_module
