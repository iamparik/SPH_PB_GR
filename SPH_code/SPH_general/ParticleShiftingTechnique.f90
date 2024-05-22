!****************************************************************************
!
!  SUBROUTINE: ParticleShiftingTechnique
!
!  PURPOSE:  Subroutine to shift particles and update concentration gradient
!****************************************************************************

subroutine ParticleShiftingTechnique(PSTtype,PSTcoeff)

use config_parameter, only:SPH_dim, itype_real_max, itype_real_min, &
        & hsml_const, dx_r, FScutoff
use particle_data, only: niac,pair_i, pair_j,eniac,epair_a, epair_s, &
    & w, dwdx, delC, surf_norm, del_gamma_as, nedge_rel_edge, ntotal, etotal, &
    & mass, rho, x, vx, itype,FreeSurfaceVar, bdryVal_seg, &
    & gamma_cont, gamma_discrt, gamma_mat, gamma_mat_inv,xi1_mat_inv


implicit none
integer(4), intent(in) :: PSTtype
real(8), intent(in) :: PSTcoeff
real(8) maxShift, grad_b_term, dstress(SPH_dim)
real(8) PSTshift, delr(SPH_dim), extra_vec(SPH_dim), dCF, w_dxr, dx_as, w_dxas, F_b(SPH_dim)
integer(4) k,d,a,b,s
real(8), DIMENSION(:,:,:), allocatable :: grad_vel

! Find the coeffecient term used to multiple the delC term
grad_b_term= PSTCoeff*hsml_const**2.D0

!maximum Particle shift admissible
maxShift=0.2D0*dx_r

extra_vec=0.D0

Allocate( grad_vel(SPH_dim, SPH_dim, ntotal),delC(SPH_dim,ntotal))
grad_vel=0.D0
delC=0.D0

call kernel(dx_r,delr,hsml_const,w_dxr,extra_vec)

do k = 1,niac
    a=pair_i(k)
    b=pair_j(k)

    ! For different PSTtype we change dCF as 
    if (PSTtype .eq. 1) then
        dCF= 1.D0
    elseif (PSTtype .eq. 2) then
        dCF=1.D0 + 0.2D0*(w(k)/w_dxr)**4.D0
    else
        dCF= 1.D0
    endif
    
    call CorrectedScaGradPtoP(delC(:,a),delC(:,b),dCF,dCF,dwdx(:,k), mass(a), mass(b), rho(a), rho(b), &
                & gamma_cont(a), gamma_discrt(a), gamma_mat(:,:,a), gamma_mat_inv(:,:,a), xi1_mat_inv(:,:,a), &
                & gamma_cont(b), gamma_discrt(b), gamma_mat(:,:,b), gamma_mat_inv(:,:,b), xi1_mat_inv(:,:,b), &
                & SPH_dim, 1, 1) ! SPH_dim, correctionFactorID, grad_type
    
    do d=1,SPH_dim 
        ! add an if condition with PCBI for inside and WCBI near bdry if needed
        call CorrectedScaGradPtoP(grad_vel(d,:,a),grad_vel(d,:,b),vx(d,a),vx(d,b),dwdx(:,k), mass(a), mass(b), rho(a), rho(b), &
                & gamma_cont(a), gamma_discrt(a), gamma_mat(:,:,a), gamma_mat_inv(:,:,a), xi1_mat_inv(:,:,a), &
                & gamma_cont(b), gamma_discrt(b), gamma_mat(:,:,b), gamma_mat_inv(:,:,b), xi1_mat_inv(:,:,b), &
                & SPH_dim, 3, 2) ! SPH_dim, correctionFactorID, grad_type
    enddo
enddo


do k= 1, eniac
    a=epair_a(k)
    s=epair_s(k)
    b = nedge_rel_edge(s)
    
    dx_as=norm2(x(:,a)-x(:,b))
    call kernel(dx_as,delr,hsml_const,w_dxas,extra_vec)
    
    ! For different PSTtype we change dCF as follows 
    if (PSTtype .eq. 1) then
        dCF= 1.D0
    elseif (PSTtype .eq. 2) then
        dCF= 1.D0 + 0.2D0*(w_dxas/w_dxr)**4.D0
    else
        dCF= 1.D0
    endif
    
    call CorrectedScaGradPtoB(delC(:,a),1.D0,dCF,del_gamma_as(:,k),  &
                        & gamma_cont(a), gamma_discrt(a), gamma_mat(:,:,a), gamma_mat_inv(:,:,a), xi1_mat_inv(:,:,a), &
                        & SPH_dim, 1, 1)
    
    do d=1,SPH_dim
        
        F_b(:) = bdryVal_seg(1:SPH_dim,s)
        call CorrectedScaGradPtoB(grad_vel(d,:,a),vx(d,a),F_b(d),del_gamma_as(:,k),  &
                        & gamma_cont(a), gamma_discrt(a), gamma_mat(:,:,a), gamma_mat_inv(:,:,a), xi1_mat_inv(:,:,a), &
                        & SPH_dim, 3, 2) ! SPH_dim, correctionFactorID, grad_type
    enddo    
enddo  

! Now perform particle shifting
if(PSTtype .ge. 1) then
    do a=1,ntotal    
        delr=0.D0
        if((itype(a) .le. itype_real_max) .and. (itype(a) .gt. itype_real_min) &
            & .and. ( FreeSurfaceVar(a) .gt. FScutoff) ) then 
            
            dstress(:) = -grad_b_term*delC(:,a)
            PSTShift = min(norm2(dstress(:)), maxShift)
            if(PSTShift .gt. 1D-10*grad_b_term) then
                delr=PSTshift*(dstress(:)/norm2(dstress(:)))
            endif
            x(:,a) = x(:,a)+ delr
            do d=1,SPH_dim
                vx(d,a)= vx(d,a) + dot_product(grad_vel(d,:,a),delr)
            enddo  
        endif  

    enddo
endif

deallocate(grad_vel)

end
