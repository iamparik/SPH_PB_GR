!****************************************************************************
!
!  SUBROUTINE: beta_discrete_factor
!
!  PURPOSE: This defines the factor beta, with (x'_s_i-x_a_i) taken out of the boundary 
!           integral elements 
!           〖β_a〗_ij=∑_s (x_s_i-x_a_i) (-∂γ_as)  
!
!   CREATED:        04/27/2022       by  PARIKSHIT BOREGOWDA
!   Last Modified:  04/27/2022        by  PARIKSHIT BOREGOWDA 
!****************************************************************************

subroutine beta_discrete_factor
    use config_geometry, only: SPH_dim
    use particle_data, only: beta_mat, eniac, epair_a, epair_s, ntotal, &
                & del_gamma_as,edge,nedge_rel_edge,x 
    implicit none
    
    real(8) k,i,j, a, s
    
    if( .NOT. allocated(beta_mat)) allocate(beta_mat(SPH_dim,SPH_dim,ntotal)) 

    beta_mat=0.D0
    
    do k=1,eniac
        a=epair_a(k)
        s=epair_s(k)
        
        do i=1,SPH_dim 
            do j=1,SPH_dim
                !the below definition is only true for one particle per edge, ingeneral this needs to be changed to include some
                ! numerical itnegration schemes
                beta_mat(i,j,a) = beta_mat(i,j,a) + (x(i,nedge_rel_edge(s))-x(i,a))*(-del_gamma_as(j,k)) 
                
            enddo
        enddo

    enddo

end