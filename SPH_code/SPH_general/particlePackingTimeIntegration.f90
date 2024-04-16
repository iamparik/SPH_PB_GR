!****************************************************************************
!
!  SUBROUTINE: particlePackingTimeIntegration
!
!  PURPOSE:  Subroutine to initialize SPH particle locations in an optimal fashion
!
!   CREATED:        04/13/2022       by  PARIKSHIT BOREGOWDA
!   Last Modified:  07/11/2023        by  PARIKSHIT BOREGOWDA 
!****************************************************************************

subroutine particlePackingTimeIntegration(packagingIterations)

use config_parameter, only: SPH_dim, pi
use config_geometry, only: print_step, save_step, hsml_const, dx_r    
use particle_data, only: nreal, w_aa, w, dwdx, &
        & gamma_discrt, gamma_cont, del_gamma_as, del_gamma, &
        & xi1_mat, beta_mat,gamma_mat,xi_cont_mat, &
        & gamma_mat_inv,xi1_mat_inv,xi_cont_mat_inv, &
        & epair_a,epair_s,  pair_i, pair_j, &
        & b_MLS, ntotal, packableParticle, &
        & FirstCutoffStep,SecondCutoffStep, &
        & xstart, vx, x

implicit none

integer(4) i, a
logical logicalcounter, packingInProgress
real(8) dt, xEdgeTemp(2,2), xrefPoint(2), xEdge_surfNorm(2),gamma_cutoff
integer(4), intent(in):: packagingIterations

logicalcounter= .true. 
packingInProgress =.true.




!intialize the start of the loop
i=0


do while (packingInProgress)
    i=i+1
    if (i .eq. 1) then 
        write(*,*)'______________________________________________'
        write(*,*)'  Starting the Packaging Alorithm'
        write(*,*)'______________________________________________'
    endif    
    
    ! Nearest Neighbor Particle search algorithm is called to find particle pairs
    ! These Particle pairs interact within the volume integral in SPH formulation
    ! In a discrete sense this requires particles with mass for summation
    call nnps_algorithm

    ! Nearest Neighbor Edge Search algorithm is called to find particles whose Kernels
    ! are truncated at the edge boundaries. This is needed to evaluate boundary integral 
    ! in SPH formulation.
    call nnes_algorithm
    
    ! Call the analytical del_gamma evaluation Boundary_integral(W*n_s)
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
    
    ! Call the function to determine kernel weights  at 0 radius, as 
    ! interpolation of 0th order functions require value of weigths at center of the kernel
    call smoothing_w0
    
    ! Call the discrete version of gamma summation(W_ab*m_b/rho_b)
    call gamma_discrete1_factor
    
    ! We define particles that can be packed initially, following Algorithm of Boregowda et. al 2024
    if (i .eq. 1) then
        ALLOCATE(packableParticle(ntotal))
        packableParticle(:) = .true.
        ! We only want to pack particles next to boundary initially
        do a = 1, ntotal
            if(gamma_cont(a) .eq. 1.D0) packableParticle(a) = .false.                
        enddo
        
    endif
    
    if( FirstCutoffStep .and. logicalCounter) then
        
        xEdgeTemp(:,:)= 0.D0
        xrefPoint(:)=0.D0
        xEdge_surfNorm(:)=0.D0
        xEdgeTemp(1,1)= -2.D0*hsml_const !2 only for wendland kernel
        xEdgeTemp(1,2)= 2.D0*hsml_const !2 only for wendland kernel
        xrefPoint(2) = -3.D0*dx_r/5.D0
        xEdge_surfNorm(2)= -1.D0
        
        call gamma_analytical_leroy_2D(xrefPoint,xEdgeTemp,xEdge_surfNorm,gamma_cutoff,hsml_const,pi)
        gamma_cutoff=1-gamma_cutoff
        
        write(*,*) "gamma_cutoff for the particle packing scheme is : " , gamma_cutoff
        
        do a = 1, ntotal
            if(gamma_cont(a) .lt. gamma_cutoff) then 
                packableParticle(a) = .false.  
            else
                packableParticle(a) = .true. 
            endif            
        enddo
        
        logicalCounter =.false.
    endif
    
    
!Call the packaging algorithm
    call ParticlePackingScheme(i)
    
    
    if( SecondCutoffStep) then
        write(*,*) "packing ends at i = ", i
        packingInProgress = .false.
    elseif(i .eq. packagingIterations) then
        write(*,*) "packing ends at i = ", i
        packingInProgress = .false.
    endif
        
    deallocate(gamma_cont,del_gamma_as, del_gamma, epair_a, epair_s, pair_i, pair_j, w,dwdx, w_aa, gamma_discrt)
  
enddo


if (allocated(xstart)) deallocate(xstart)

if (allocated(vx)) deallocate(vx)


end