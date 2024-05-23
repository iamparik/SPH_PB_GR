!****************************************************************************
!
!  SUBROUTINE: EulerIntegration_fluid
!
!  PURPOSE:  Subroutine to march in time for fluid flow problem
!
!   CREATED:        08/20/2022       by  PARIKSHIT BOREGOWDA
!   Last Modified:  07/11/2023       by  PARIKSHIT BOREGOWDA 
!****************************************************************************

subroutine ParticlePackingSchemeNew(iterStep)

    use config_parameter, only: SPH_dim,dataPackingPath, itype_real_max, itype_real_min, &
        & c_sound, hsml_const, packagingIterations,dx_r, &
        & save_step, PSTCoeff, delC_Cap
    use particle_data, only: avgVol, ntotal, nreal, etotal, itype,nedge_rel_edge,&
        & dstress, x, vx, del_gamma, gamma_cont, PE,KE,TE, mass, max_vel, gamma_cont,gamma_discrt, &
        & niac, pair_i, pair_j, dwdx, eniac, epair_a,epair_s, del_gamma_as, mass, rho, surf_norm, &
        & delC, delCAvg, delCMax, delCL2, xstart, ga_Max,ga_Avg,ga_L2, g_a_min, g_a_max, packableParticle, &
        & FirstCutoffStep,SecondCutoffStep, PP_variable, PP_Variable_prev
    
    
    implicit none
    integer(4), intent(in) :: iterStep
    integer(4) :: i, j, k,a,b,s,d
    real(8)::grad_b_term, maxShift, PSTShift, ps_pa, energyStepT, TPD , dtsress(SPH_dim)
    character (40) :: xname, x2name
    
    ! Define the coeffecient of particle shifting technique used for packing particles
    grad_b_term= PSTCoeff*hsml_const**2.D0!
    
    maxShift=0.5D0*dx_r !0.2D0* hsml_const !0.25D0*dx_r!0.25*dx_r
    
    extra_vec=0.D0
    
    if( iterstep .eq. 1) then
        Allocate(xStart(SPH_dim,ntotal))
        xStart=x        
    endif

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
        if((itype(a) .le. itype_real_max) .and. (itype(a) .gt. itype_real_min)) then 
            
            dstress(:) = -grad_b_term*(delC(:,a)+bdry_push(:,a))
            PSTShift = min(norm2(dstress(:)), maxShift)
            if(PSTShift .gt. 1D-10*grad_b_term) then
                delr=PSTshift*(dstress(:)/norm2(dstress(:)))
            endif
            x(:,a) = x(:,a)+ delr
        endif  
    enddo

    deallocate(bdry_push)
    

end
     
