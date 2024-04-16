!****************************************************************************
!
!  SUBROUTINE: fncnApproxOperator
!
!  PURPOSE: This defines GAMMA matrix, which is the correction factor to be used
!           for calculating first oder reviatives
!           ⟨〈f(x)〉⟩≅1/γ ∑_b〖f(x_b) W(‖x-x_b ‖,h) m_b/ρ_b 〗 
!
!   CREATED:        04/27/2022       by  PARIKSHIT BOREGOWDA
!   Last Modified:  12/04/2022        by  PARIKSHIT BOREGOWDA 
!****************************************************************************
    
subroutine fncnApproxOperator(fncn_a,fncn_b,mass_b,rho_b,weight)
    implicit none


    real(8), intent(in):: fncn_b, mass_b,rho_b,weight
    real(8):: fncn_a


    fncn_a = fncn_a + fncn_b*weight*mass_b/rho_b
end



