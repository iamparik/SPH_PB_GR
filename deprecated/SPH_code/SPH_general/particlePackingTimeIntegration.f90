!****************************************************************************
!
!  SUBROUTINE: particlePackingTimeIntegration
!
!  PURPOSE:  Subroutine to initialize SPH particle locations in an optimal fashion
!
!   CREATED:        04/13/2022       by  PARIKSHIT BOREGOWDA
!   Last Modified:  07/11/2023        by  PARIKSHIT BOREGOWDA 
!****************************************************************************

subroutine particlePackingTimeIntegration(quick_converge_step2C)

use config_parameter, only: SPH_dim, pi, DataConfigPath, &
    & print_step, save_step, hsml_const, dx_r, itype_real_min, itype_real_max    
use particle_data, only: nreal, w_aa, w, dwdx, &
        & gamma_discrt, gamma_cont, del_gamma_as, del_gamma, &
        & xi1_mat, beta_mat,gamma_mat, &
        & gamma_mat_inv,xi1_mat_inv, &
        & epair_a,epair_s, eniac, pair_i, pair_j, niac, nedge_rel_edge,  &
        & ntotal, xstart, x, mass, rho, itype, delC, xstart, &
        & edge,etype, hsml, vol, surf_norm

implicit none


logical, intent(in) :: quick_converge_step2C
real(8) dt, xEdgeTemp(2,2), xrefPoint(2), xEdge_surfNorm(2),gamma_cutoff, scale_k, &
    &   dCF, TPD, delC_avg
real(8) PP_Variable_prev, PP_Variable, grad_b_term, maxShift, w_dxr, delr(SPH_dim), &
    & extra_vec(SPH_dim), dstress(SPH_dim), dx_as, w_dxas, ps_pa, PSTShift
integer(4) iterstep, a, b, k,s,d, cutoff_step
logical packing_in_progress
logical, DIMENSION(:),ALLOCATABLE :: packableParticle
real(8),DIMENSION(:,:),ALLOCATABLE :: bdry_push
character (40) :: x2name


call sml_mult_factor(scale_k)

Allocate(xStart(SPH_dim,ntotal))
xStart=x   

!intialize the start of the loop
iterstep=0

!intialize cutoff_step as zero
cutoff_step=0

!initialize packing_in_progress to true 
packing_in_progress = .true.

! start the particle packing loop
do while (packing_in_progress)
    iterstep= iterstep+1
    if (iterstep .eq. 1) then 
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
    if (iterstep .eq. 1) then
        ALLOCATE(packableParticle(ntotal))
        packableParticle(:) = .true.
        
        ! We only want to pack particles next to boundary initially
        do a = 1, ntotal
            if(gamma_cont(a) .eq. 1.D0) packableParticle(a) = .false.                
        enddo
        
        !Initialize the previous particle packing variable to be zero
        PP_Variable_prev=0.D0

    endif
    
    if( cutoff_step .eq. 1 ) then
        
        xEdgeTemp(:,:)= 0.D0
        xrefPoint(:)=0.D0
        xEdge_surfNorm(:)=0.D0
        xEdgeTemp(1,1)= -scale_k*hsml_const 
        xEdgeTemp(1,2)= scale_k*hsml_const 
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
        
        cutoff_step = 2
    endif
    
    
!Call the packaging algorithm
    ! Define the coeffecient of particle shifting technique used for packing particles
    grad_b_term= 0.5D0*hsml_const**2.D0!
    
    maxShift=0.5D0*dx_r !0.2D0* hsml_const !0.25D0*dx_r!0.25*dx_r

    Allocate( bdry_push(SPH_dim, ntotal),delC(SPH_dim,ntotal))
    bdry_push=0.D0
    delC=0.D0

    call kernel(dx_r,delr,hsml_const,w_dxr,extra_vec)

    do k = 1,niac
        a=pair_i(k)
        b=pair_j(k)

        dCF= 1.D0
    
        call CorrectedScaGradPtoP(delC(:,a),delC(:,b),dCF,dCF,dwdx(:,k), mass(a), mass(b), rho(a), rho(b), &
                    & gamma_cont(a), gamma_discrt(a), gamma_mat(:,:,a), gamma_mat_inv(:,:,a), xi1_mat_inv(:,:,a), &
                    & gamma_cont(b), gamma_discrt(b), gamma_mat(:,:,b), gamma_mat_inv(:,:,b), xi1_mat_inv(:,:,b), &
                    & SPH_dim, 1, 1) ! SPH_dim, correctionFactorID, grad_type
    enddo


    do k= 1, eniac
        a=epair_a(k)
        s=epair_s(k)
        b = nedge_rel_edge(s)
    
        dx_as=norm2(x(:,a)-x(:,b))
        call kernel(dx_as,delr,hsml_const,w_dxas,extra_vec)  
        
        dCF= 1.D0
    
        call CorrectedScaGradPtoB(delC(:,a),1.D0,dCF,del_gamma_as(:,k),  &
                            & gamma_cont(a), gamma_discrt(a), gamma_mat(:,:,a), gamma_mat_inv(:,:,a), xi1_mat_inv(:,:,a), &
                            & SPH_dim, 1, 1)
        
        ! Add boundary force function below
        ! A step function twice the force is used here for particles
        ! too close to the boundary
        if ((dx_as .le. dx_r/2.D0) ) then 
            ! use compressive bdry force only if the bdry is less than the particle radius 
            ! but since bdry particles
            !if (dot_product(x(:,a)-x(:,b),surf_norm(:,s)) .le. dx_r/2.D0) ps_pa = 2.D0
            ps_pa = 2.D0
        else
            ps_pa = 1.D0
        endif
        
        dCF=0.5D0*(ps_pa -1.D0)
    
        call CorrectedScaGradPtoB(bdry_push(:,a),1.D0,dCF,del_gamma_as(:,k),  &
                            & gamma_cont(a), gamma_discrt(a), gamma_mat(:,:,a), gamma_mat_inv(:,:,a), xi1_mat_inv(:,:,a), &
                            & SPH_dim, 1, 1)
    enddo  

    !------------------------Particle shift calcualtion from force terms -----------------------------!
    do a=1,ntotal    
        delr=0.D0
        if((itype(a) .le. itype_real_max) .and. (itype(a) .gt. itype_real_min) .and.  packableParticle(a)) then 
            
            dstress(:) = -grad_b_term*(delC(:,a)+bdry_push(:,a))
            PSTShift = min(norm2(dstress(:)), maxShift)
            if(PSTShift .gt. 1D-10*grad_b_term) then
                delr=PSTshift*(dstress(:)/norm2(dstress(:)))
            endif
            x(:,a) = x(:,a)+ delr
        endif  
    enddo

    deallocate(bdry_push)
    
    
!------------------------Calculate parameters useful to simulation -----------------------------!  
 
    !Initialize total particle displacement to zero and average concgradient to zero
    TPD=0.D0
    delC_avg = 0.D0
    do a=1,ntotal        
        if((itype(a) .le. itype_real_max) .and. (itype(a) .gt. itype_real_min)) then
            TPD = TPD + norm2(xStart(:,a)-x(:,a))/nreal
            delC_avg= delC_avg + norm2(delC(:,a))/nreal
        endif
    enddo
    
    call outputPacking(iterstep,100,TPD,delC_avg) !input: iterStep, saveStep, TPD
  
    ! Now we check if for every 1000 steps if the particle packing variable
    ! is unchanged compared to its previous value
    if(mod(iterstep,1000) .eq. 0) then
        ! Select the particle packing variable (like TPD)
        PP_variable =TPD
        
        if((abs(PP_Variable - PP_Variable_prev)) .lt. 1.D-2*PP_Variable) then
            cutoff_step = cutoff_step + 1
        elseif((cutoff_step .eq. 2) .and. quick_converge_step2C) then
            cutoff_step = cutoff_step + 1
        endif
        
        PP_variable_prev = PP_variable
    endif
    
    deALLOCATE(delC)
    

        
    deallocate(gamma_cont,del_gamma_as, del_gamma, epair_a, epair_s, pair_i, pair_j, w,dwdx, w_aa, gamma_discrt)
  
    if( cutoff_step .eq. 3) then
        write(*,*) "packing ends at iterstep = ", iterstep
        packing_in_progress = .False.
    endif
        
enddo

deallocate(xstart)
    
    write(x2name, '(A,A)') DataConfigPath,'/input_PP.dat'
    open (1, file = x2name)
    

    ! write particle data to the file input_PP.dat
    do a=1,nreal
        write(1,'(I10,4(e22.10,2x))') itype(a), (x(d, a), d=1,SPH_dim), vol(a)
    enddo
    
    close(1)

end