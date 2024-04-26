!****************************************************************************
!
!  SUBROUTINE: EulerSchemeSPHOperations
!
!  PURPOSE:  Subroutine to test different SPH operations for Euler Integration Scheme
!
!   CREATED:        04/13/2022       by  PARIKSHIT BOREGOWDA
!   Last Modified:  08/05/2022        by  PARIKSHIT BOREGOWDA 
!****************************************************************************

subroutine TimeIndependentSPHOperations

use config_parameter, only: NumericalSimCase             
use particle_data, only:itimestep, w_aa, w, dwdx, &
        & gamma_discrt, gamma_cont, del_gamma_as, del_gamma, &
        & xi1_mat, beta_mat,gamma_mat,xi_cont_mat, &
        & gamma_mat_inv,xi1_mat_inv,xi_cont_mat_inv, &
        & pBC_edges, pBC_epair_a, pBC_epair_s, &
        & epair_a,epair_s,  pair_i, pair_j

implicit none


! The configaration of the input geometry and features are verified by
! printing the input data to be graphically displayed
    call output_initial_config

! Subroutine to recaclulate kernels, kernel derivative and 
! all correction factors that dpend on position of particle
    call PositionDependentFactors
    
! Account for periodic boundary, by settign up periodic bc and copying particles   
    if (Allocated(pBC_edges)) then
        call output_initial_config
    endif
    
!Run suitable time independent numerical Simulation    
    if(NumericalSimCase .eq. 0) then
        call testForSPHApproxLapFormulations
    endif
  
    
!deallocate all correction factors 
    deallocate(gamma_cont, epair_a, epair_s, pair_i, pair_j, w,dwdx)
    
    deallocate(w_aa, gamma_discrt,del_gamma_as, del_gamma, xi1_mat, beta_mat,gamma_mat,xi_cont_mat, &
                   & gamma_mat_inv,xi1_mat_inv,xi_cont_mat_inv )


    if (Allocated(pBC_edges)) then
        deallocate(pBC_epair_a, pBC_epair_s)
    endif


end