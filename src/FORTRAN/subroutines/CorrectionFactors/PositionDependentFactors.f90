!****************************************************************************
!
!  SUBROUTINE: PositionDependentFactors
!
!  PURPOSE:  Subroutine to recaclulate kernels, kernel derivative and 
!               all correction factors that dpend on position of particle
!
!   CREATED:        06/29/2023       by  PARIKSHIT BOREGOWDA
!   Last Modified:  06/29/2023         by  PARIKSHIT BOREGOWDA 
!****************************************************************************

subroutine PositionDependentFactors
    
    use particle_data, only: pBC_edges,nperiodic,pBC_duplicate_pair, w_aa, &
                    & gamma_discrt, del_gamma, gamma_cont, xi1_mat, beta_mat, gamma_mat
    use config_parameter, only: SPH_dim
    
    implicit none
    
    integer(4):: i,a 
    
    ! Nearest Neighbor Particle search algorithm is called to find particle pairs
    ! These Particle pairs interact within the volume integral in SPH formulation
    ! In a discrete sense this requires particles with mass for summation
    call nnps_algorithm

    ! Nearest Neighbor Edge Search algorithm is called to find particles whose Kernels
    ! are truncated at the edge boundaries. This is needed to evaluate boundary integral 
    ! in SPH formulation.
    call nnes_algorithm
        
    if (Allocated(pBC_edges)) then
            call PeriodicBCsetup2D
    endif

    ! Call the function to determine kernel weights  at 0 radius, as 
    ! interpolation of 0th order functions require value of weigths at center of the kernel
    call smoothing_w0
 
    ! update w_aa for periodic particles
    if (Allocated(pBC_edges)) then
        call PeriodicParameter(w_aa)
    endif
        
        
    !! Call relevant correction factors

    ! Call the discrete version of gamma summation(W_ab*m_b/rho_b)
    call gamma_discrete1_factor
    
    ! update gamma_discrt for periodic particles
    if (Allocated(pBC_edges)) then
        call PeriodicParameter(gamma_discrt)
    endif

    ! Call the analytical del_gamm evaluation Boundary_integral(W*n_s)
    ! However, this analytical determination is only for the discretized boundary with linear boundary elements.
    ! So to find the del_gamma for the real boundary, the boundary discretization
    ! needs to refined. So, a Circular boundary's  del_gamma will only be reproduced when the n polygon representing
    ! tends to infinity vertices polygon. However, for a given polygon boudnary approxiamtion, analytical del_gamm calculated below,
    ! gives the accurate value of del_gamma for that n-polygon.
    call dgamma_analytical
    
    ! update gamma_cont for periodic particles
    if (Allocated(pBC_edges)) then  
        do  a= 1,nperiodic
            del_gamma(:,pBC_duplicate_pair(2,a))=del_gamma(:,pBC_duplicate_pair(1,a))
        enddo
    endif


    ! Value of continuos gamma can be evaluated either from the del_gamma like in Ferrand's work,
    ! or by defining a divergence oeprator and evaluating gamma with an analyical formula
    ! call gamma_continuous_Ferrand(dt)
    call gamma_continuous_leroy

    ! update gamma_cont for periodic particles
    if (Allocated(pBC_edges)) then
        call PeriodicParameter(gamma_cont)
    endif
 
    !Xi is the factor that is used to calculate Xi〖 ξ_a〗_ij=∑_b〖(〖x_b〗_i-〖x_a〗_i )  m_b/ρ_b  ∂_j_a W_ab 〗
    call xi_discrete1_factor
    
    !Update Xi for Periodic particles
    if (Allocated(pBC_edges)) then
        do  a= 1,nperiodic
            xi1_mat(:,:,pBC_duplicate_pair(2,a))=xi1_mat(:,:,pBC_duplicate_pair(1,a))
        enddo
    endif

    !Beta is the factor that is used to calculate   〖β ̃_a〗_ij=∑_s∫_(∂(Ω ∩ Ω_w )_s)〖(〖x_s^'〗_i-〖x_a〗_i ) W_(as^' ) n_(s_j )  dS^' 〗 
    call beta_discrete_factor
    
    !Update beta for Periodic particles
    if (Allocated(pBC_edges)) then
        do  a= 1,nperiodic
            beta_mat(:,:,pBC_duplicate_pair(2,a))=beta_mat(:,:,pBC_duplicate_pair(1,a))
        enddo   
    endif

    
    !We call the factor Gamma_mtartix,  〖¯Γ_a〗_ij=〖β ̃_a〗_ij+〖¯ξ_a〗_ij  
    call gamma_mat_dicrete_factor
    
    !Update gamma_mat for Periodic particles
    if (Allocated(pBC_edges)) then
        do  a= 1,nperiodic
            gamma_mat(:,:,pBC_duplicate_pair(2,a))=gamma_mat(:,:,pBC_duplicate_pair(1,a))
        enddo      
    endif
    
    !Find inverse of factors that can be used as correction factors.
    call CorrectionFactorMatrixInversion

            
    ! For periodic BC's copy all factors    
    !if (Allocated(pBC_edges)) then
    !    call PeriodicBCcorrectionFactors
    !endif
        
end