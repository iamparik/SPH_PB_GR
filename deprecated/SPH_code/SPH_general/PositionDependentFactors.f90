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
    
    use particle_data, only: pBC_edges
    
    implicit none
    

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
 

        
        
    !! Call relevant correction factors

    ! Call the discrete version of gamma summation(W_ab*m_b/rho_b)
    call gamma_discrete1_factor



    ! Call the analytical del_gamm evaluation Boundary_integral(W*n_s)
    ! However, this analytical determination is only for the discretized boundary with linear boundary elements.
    ! So to find the del_gamma for the real boundary, the boundary discretization
    ! needs to refined. So, a Circular boundary's  del_gamma will only be reproduced when the n polygon representing
    ! tends to infinity vertices polygon. However, for a given polygon boudnary approxiamtion, analytical del_gamm calculated below,
    ! gives the accurate value of del_gamma for that n-polygon.
    call dgamma_analytical
    


    ! Value of continuos gamma can be evaluated either from the del_gamma like in Ferrand's work,
    ! or by defining a divergence oeprator and evaluating gamma with an analyical formula
    ! call gamma_continuous_Ferrand(dt)
    call gamma_continuous_leroy

 
    !Xi is the factor that is used to calculate Xi〖 ξ_a〗_ij=∑_b〖(〖x_b〗_i-〖x_a〗_i )  m_b/ρ_b  ∂_j_a W_ab 〗
    call xi_discrete1_factor
    

    !Beta is the factor that is used to calculate   〖β ̃_a〗_ij=∑_s∫_(∂(Ω ∩ Ω_w )_s)〖(〖x_s^'〗_i-〖x_a〗_i ) W_(as^' ) n_(s_j )  dS^' 〗 
    call beta_discrete_factor
  

    
    !We call the factor Gamma_mtartix,  〖¯Γ_a〗_ij=〖β ̃_a〗_ij+〖¯ξ_a〗_ij  
    call gamma_mat_dicrete_factor
    
    
    !Find inverse of factors that can be used as correction factors.
    call CorrectionFactorMatrixInversion

            
    ! For periodic BC's copy all factors    
    if (Allocated(pBC_edges)) then
        call PeriodicBCcorrectionFactors
    endif
        
end