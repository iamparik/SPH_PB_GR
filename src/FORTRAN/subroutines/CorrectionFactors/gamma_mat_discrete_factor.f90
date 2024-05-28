!****************************************************************************
!
!  SUBROUTINE: gamma_mat_discrete_factor
!
!  PURPOSE: This defines GAMMA matrix, which is the correction factor to be used
!           for calculating first oder reviatives
!           〖¯Γ_a〗_ij=〖β ̃_a〗_ij+〖¯ξ_a〗_ij   
!
!   CREATED:        04/27/2022       by  PARIKSHIT BOREGOWDA
!   Last Modified:  04/27/2022        by  PARIKSHIT BOREGOWDA 
!****************************************************************************
    
subroutine gamma_mat_dicrete_factor
   use config_parameter, only: SPH_dim
   use particle_data, only: beta_mat,xi1_mat, gamma_mat, ntotal

   implicit none        

    if( .NOT. allocated(gamma_mat)) allocate(gamma_mat(SPH_dim,SPH_dim,ntotal)) 

   gamma_mat=0.D0
   
   gamma_mat=beta_mat+xi1_mat


   
   
end

