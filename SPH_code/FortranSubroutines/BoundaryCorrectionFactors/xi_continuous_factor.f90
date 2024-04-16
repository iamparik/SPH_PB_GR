!****************************************************************************
!
!  SUBROUTINE: xi_continuous_factor
!
!  PURPOSE: This defines Xi matrix, which is the correction factor to be used
!           for calculating first oder reviatives
!           〖ξ_a〗=〖I〗γ -〖β_a〗
!
!   CREATED:        06/04/2022       by  PARIKSHIT BOREGOWDA
!   Last Modified:  09/20/2022        by  PARIKSHIT BOREGOWDA 
!****************************************************************************
    
subroutine xi_continuous_factor
   use config_parameter, only: SPH_dim
   use particle_data, only: beta_mat, gamma_cont, xi_cont_mat, ntotal

   implicit none   
   
   integer(4) a,d
   real(8), dimension(:,:,:),allocatable:: gamma_cont_mat
   
   allocate(gamma_cont_mat(SPH_dim,SPH_dim,ntotal))
   
   gamma_cont_mat=0.D0
   
   do a=1,ntotal
        do d=1,SPH_dim
                gamma_cont_mat(d,d,a)= gamma_cont(a)
        enddo
    enddo
   

    if( .NOT. allocated(xi_cont_mat)) allocate(xi_cont_mat(SPH_dim,SPH_dim,ntotal)) 

   xi_cont_mat=0.D0
   
   xi_cont_mat=gamma_cont_mat - beta_mat
   
   deallocate(gamma_cont_mat)

end

