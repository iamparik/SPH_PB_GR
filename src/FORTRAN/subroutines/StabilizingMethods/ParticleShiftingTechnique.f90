!****************************************************************************
!
!  SUBROUTINE: ParticleShiftingTechnique
!
!  PURPOSE:  Subroutine to shift particles and update concentration gradient
!****************************************************************************

subroutine ParticleShiftingTechnique(PSTtype,PSTcoeff, PST_Step_, itimestep_)

use config_parameter, only:SPH_dim, hsml_const, dx_r, FScutoff
use particle_data, only: niac,pair_i, pair_j,eniac,epair_a, epair_s, &
    & w, dwdx, delC, surf_norm, del_gamma_as, ntotal, nreal, etotal, &
    & mass, rho, x, vx, itype,free_surf_particle, bdryVal_seg, &
    & gamma_cont, gamma_discrt, gamma_mat, gamma_mat_inv,xi1_mat_inv, x_prev


implicit none
integer(4), intent(in) :: PSTtype, PST_Step_, itimestep_
real(8), intent(in) :: PSTcoeff
real(8) maxShift, grad_b_term, dstress(SPH_dim)
real(8) PSTshift, delr(SPH_dim), extra_vec(SPH_dim), dCF, w_dxr, dx_as, w_dxas, F_b(SPH_dim)
integer(4) k,d,a,b,s
real(8), DIMENSION(:,:,:), allocatable :: grad_vel
logical :: perform_shifting

! Find the coeffecient term used to multiple the delC term
grad_b_term= PSTCoeff*hsml_const**2.D0

!maximum Particle shift admissible
maxShift=0.2D0*dx_r

extra_vec=0.D0

Allocate( delC(SPH_dim,ntotal))
delC=0.D0

if ((mod(itimestep_,PST_step_).eq.0) .and. PSTtype .ge. 1) then
    perform_shifting = .true.
else
    perform_shifting = .false.
endif



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
    
        call CorrectedScaGradPtoP(delC(:,a),delC(:,b),a,b,dCF,dCF,dwdx(:,k), mass(a), mass(b), rho(a), rho(b), &
                    & gamma_cont(a), gamma_discrt(a), gamma_mat(:,:,a), gamma_mat_inv(:,:,a), xi1_mat_inv(:,:,a), &
                    & gamma_cont(b), gamma_discrt(b), gamma_mat(:,:,b), gamma_mat_inv(:,:,b), xi1_mat_inv(:,:,b), &
                    & free_surf_particle(a),free_surf_particle(b),SPH_dim, 1, 1) ! SPH_dim, correctionFactorID, grad_type
    
    enddo

    do k= 1, eniac
        a=epair_a(k)
        s=epair_s(k)

        ! For different PSTtype we change dCF as follows 
        if (PSTtype .eq. 1) then
            dCF= 1.D0
        elseif (PSTtype .eq. 2) then
            dx_as=dot_product(x(:,a), surf_norm(:,s)) 
            call kernel(dx_as,delr,hsml_const,w_dxas,extra_vec)
            dCF= 1.D0 + 0.2D0*(w_dxas/w_dxr)**4.D0
        else
            dCF= 1.D0
        endif
    
        call CorrectedScaGradPtoB(delC(:,a),a,s,1.D0,dCF,del_gamma_as(:,k),  &
                            & gamma_cont(a), gamma_discrt(a), gamma_mat(:,:,a), gamma_mat_inv(:,:,a), xi1_mat_inv(:,:,a), &
                            & free_surf_particle(a),SPH_dim, 1, 1)
    

    enddo  

    ! Now perform particle shifting
if(perform_shifting) then 
    Allocate( x_prev(SPH_dim,ntotal))
    do a=1,nreal    
        delr=0.D0
        x_prev(:,a)=x(:,a)

        dstress(:) = -grad_b_term*delC(:,a)
        PSTShift = min(norm2(dstress(:)), maxShift)
        if(PSTShift .gt. 1D-10*grad_b_term) then
            delr=PSTshift*(dstress(:)/norm2(dstress(:)))
        endif
        
        ! account for free surface
        delr=dble(1-free_surf_particle(a))*delr
        
        x(:,a) = x_prev(:,a)+ delr     

    enddo
endif


end
