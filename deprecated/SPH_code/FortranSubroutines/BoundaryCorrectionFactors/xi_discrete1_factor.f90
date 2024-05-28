!****************************************************************************
!
!  SUBROUTINE: xi_discrete1_factor
!
!  PURPOSE: This defines the normalization 
!           Xi〖 ξ_a〗_ij=∑_b〖(〖x_b〗_i-〖x_a〗_i )  m_b/ρ_b  ∂_j_a W_ab 〗
!
!   CREATED:        04/27/2022       by  PARIKSHIT BOREGOWDA
!   Last Modified:  04/27/2022        by  PARIKSHIT BOREGOWDA 
!****************************************************************************
subroutine xi_discrete1_factor
    use config_parameter, only: SPH_dim
    use particle_data, only: xi1_mat, niac, pair_i, pair_j, ntotal, &
                & dwdx,mass, rho, x 
    implicit none

    real(8) k,i,j, a, b

    if( .NOT. allocated(xi1_mat)) allocate(xi1_mat(SPH_dim,SPH_dim,ntotal)) 

   xi1_mat=0.D0

    do k=1,niac
        a=pair_i(k)
        b=pair_j(k)
        
        do i=1,SPH_dim 
            do j=1,SPH_dim
                xi1_mat(i,j,a) = xi1_mat(i,j,a) + (x(i,b)-x(i,a))*mass(b)*dwdx(j,k)/rho(b) 
                xi1_mat(i,j,b) = xi1_mat(i,j,b) + (x(i,a)-x(i,b))*mass(a)*(-dwdx(j,k))/rho(a)
            enddo          
        
        enddo
        if(isNAN(xi1_mat(1,1,a))) pause

    enddo

end
    