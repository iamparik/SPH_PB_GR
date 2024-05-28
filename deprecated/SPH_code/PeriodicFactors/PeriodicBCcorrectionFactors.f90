!****************************************************************************
!
!  SUBROUTINE: PeriodicBCcorrectionFactors
!
!  PURPOSE: For simulations having periodic boundary condition
!           the duplicate particles from periodicity are updated with
!           their original particle's correction factor
!
!   CREATED:        08/12/2022       by  PARIKSHIT BOREGOWDA
!   Last Modified:  08/12/2022       by  PARIKSHIT BOREGOWDA 
!****************************************************************************
subroutine PeriodicBCcorrectionFactors
    use particle_data, only: nperiodic,pBC_duplicate_pair,&
        & gamma_discrt, gamma_cont, gamma_cont_prev, del_gamma, &
        & xi_cont_mat, xi_cont_mat_inv, xi1_mat, xi1_mat_inv, &
        & beta_mat, gamma_mat, gamma_mat_inv, w_aa, b_MLS
    implicit none

    integer a


    do a = 1, nperiodic
        gamma_discrt(pBC_duplicate_pair(2,a))=gamma_discrt(pBC_duplicate_pair(1,a))
        gamma_cont(pBC_duplicate_pair(2,a))= gamma_cont(pBC_duplicate_pair(1,a))
        if(allocated(gamma_cont_prev)) gamma_cont_prev(pBC_duplicate_pair(2,a))=gamma_cont_prev(pBC_duplicate_pair(1,a))
        del_gamma(:,pBC_duplicate_pair(2,a))=del_gamma(:,pBC_duplicate_pair(1,a))
        xi_cont_mat(:,:,pBC_duplicate_pair(2,a))=xi_cont_mat(:,:,pBC_duplicate_pair(1,a))
        xi_cont_mat_inv(:,:,pBC_duplicate_pair(2,a))=xi_cont_mat_inv(:,:,pBC_duplicate_pair(1,a))
        xi1_mat(:,:,pBC_duplicate_pair(2,a))= xi1_mat(:,:,pBC_duplicate_pair(1,a))
        xi1_mat_inv(:,:,pBC_duplicate_pair(2,a))= xi1_mat_inv(:,:,pBC_duplicate_pair(1,a))
        beta_mat(:,:,pBC_duplicate_pair(2,a))=beta_mat(:,:,pBC_duplicate_pair(1,a))
        gamma_mat(:,:,pBC_duplicate_pair(2,a))=gamma_mat(:,:,pBC_duplicate_pair(1,a))
        gamma_mat_inv(:,:,pBC_duplicate_pair(2,a))=gamma_mat_inv(:,:,pBC_duplicate_pair(1,a))
        w_aa(pBC_duplicate_pair(2,a))=w_aa(pBC_duplicate_pair(1,a))
        if(allocated(b_MLS)) b_MLS(:,pBC_duplicate_pair(2,a))=b_MLS(:,pBC_duplicate_pair(1,a))
    enddo

    
end    
