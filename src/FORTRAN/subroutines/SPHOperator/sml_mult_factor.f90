!****************************************************************************
!
!  SUBROUTINE: sml_mutl_factor 
!
!  PURPOSE:  This subroutine determines the multiplaction factor to the 
!           smoothing length that determines the compactness of the kernel 
!            the factors value depends onthe smoothing kernel used for simulation
!
!   CREATED:        9/13/2017       by  GR Lab
!   Last Modified:  04/18/2022       by  PARIKSHIT BOREGOWDA 
!****************************************************************************

subroutine sml_mult_factor(scale_k)

use config_parameter, only: skf
implicit none

real(8) scale_k

    if (skf.eq.1) then 
        scale_k = 2.0D0           
    else if (skf.eq.2) then 
        scale_k = 2.0D0
    else if (skf.eq.3) then 
        scale_k = 3.0D0 
    else if (skf.eq.4) then 
        scale_k = 2.0D0  
    else if (skf .eq. 5) then
        scale_k = 2.0D0
    endif 

end