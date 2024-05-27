!****************************************************************************
!
!  SUBROUTINE: EulerIntegration_fluid
!
!  PURPOSE:  Subroutine to march in time for fluid flow problem
!
!   CREATED:        08/20/2022       by  PARIKSHIT BOREGOWDA
!   Last Modified:  07/11/2023       by  PARIKSHIT BOREGOWDA 
!****************************************************************************

subroutine ParticlePackingScheme(iterStep)

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
    real(8)::grad_b_term, maxShift, PSTShift, r_c, Pm1, Pm2, energyStepT, TPD
    character (40) :: xname, x2name
    
    ! Define the coeffecient of particle shifting technique used for packing particles
    grad_b_term= PSTCoeff*hsml_const**2.D0!
    
    maxShift=0.5D0*dx_r !0.2D0* hsml_const !0.25D0*dx_r!0.25*dx_r
    
    if( iterstep .eq. 1) then
        Allocate(xStart(SPH_dim,ntotal))
        xStart=x        
    endif
    
    
    ALLOCATE( dstress(SPH_dim,ntotal)) !delC_diff(SPH_dim,ntotal)
    
    if(iterStep .eq. 1)  then
        
        ALLOCATE(vx(SPH_dim,ntotal))
        vx=0.D0

    endif
    
    
    !delC_diff=0.D0
    !delC=0.D0
    dstress=0.D0

    !------------------------Force term calcualtion -----------------------------!
    call ConcGradient
    
    do k= 1, eniac
        a=epair_a(k)
        s=epair_s(k)
        b = nedge_rel_edge(s)       
        
        Pm1=2.D0
        Pm2=2.D0
        if (norm2(x(:,a)-x(:,b)) .le. dx_r/2.D0) then
            !r_c=0.5D0-min(0.5D0,abs(dot_product(x(:,a)-x(:,b), surf_norm(:,s)))/dx_r)
            r_c=min(0.5D0,abs(dot_product(x(:,a)-x(:,b), surf_norm(:,s)))/dx_r)
            r_c= 0.5D0*(Pm2 - r_c*(Pm2-Pm1)/0.5D0 -1.D0)
        else
            r_c=0.D0
        endif
            
                       
        do d=1,SPH_dim
            call fncnApproxOperator(dstress(d,a),-r_c,1.D0,1.D0,del_gamma_as(d,k)/gamma_cont(a))
        enddo  

    enddo
    
    !------------------------Particle shift calcualtion from force terms -----------------------------!
    do a=1,ntotal
        !dstress(:,a)= -grad_b_term*dstress(:,a) - lin_damp_term*vx(:,a) + (1/gamma_cont(a))*grad_b_term*del_gamma(:,a)/gamma_cont(a)
        !delC(:,a)= -delC(:,a) !- lin_damp_term*vx(:,a) !+grad_b_term*del_gamma(:,a)
        if((itype(a) .le. itype_real_max) .and. (itype(a) .gt. itype_real_min)) then
            vx(:,a)=0.D0
            if ( packableParticle(a) ) then    
                dstress(:,a)= -grad_b_term*(delC(:,a) +dstress(:,a))                
                PSTShift=min(norm2(dstress(:,a)), maxShift)                
                if(PSTShift .le. 1D-10*grad_b_term) then
                !if( gamma_cont(a) .le. 0.62D0) then
                    vx(:,a)=0.D0
                else
                    vx(:,a)=PSTshift*(dstress(:,a)/norm2(dstress(:,a)))
                endif   
                x(:,a) = x(:,a) + vx(:,a)             
            endif 
        endif        
    enddo

!------------------------Calculate parameters useful to simulation -----------------------------!  
 
    
    TPD=0.D0
    do a=1,ntotal        
        if((itype(a) .le. itype_real_max) .and. (itype(a) .gt. itype_real_min)) then
            TPD = TPD + norm2(xStart(:,a)-x(:,a))/nreal
        endif
    enddo
    
    call outputPacking(iterstep,100,TPD)
  
    ! Now we check if for every 1000 steps if the TPD is unchanged compared to its previous value
    if( .not. FirstCutoffStep) then
        if(mod(iterstep,1000) .eq. 0) then
            PP_variable = TPD
            if((abs(PP_Variable - PP_Variable_prev)) .lt. 1.D-2*PP_Variable) then
                FirstCutoffStep= .true.
                write(*,*) "pp_variable = ", PP_Variable
            endif
            PP_Variable_prev=PP_Variable
        endif
        if(FirstCutoffStep) PP_Variable_prev=0.D0
    else
        
        if(mod(iterstep,1000) .eq. 0) then
        !We run packing for 1000 simulations further and stop packing scheme
            SecondCutoffStep= .true.
            
        !If below condition is used then the above condition needs to be commented
        ! As below condition requires TPD to be unchanged for every 1000 steps to
        ! to stop packing    
            !if((abs(PP_Variable - PP_Variable_prev)) .lt. 1.D-2*PP_Variable) SecondCutoffStep= .true.
            
            PP_Variable_prev=PP_Variable
        endif
    endif
    
    
    deALLOCATE(delC, dstress)
    !if(iterStep .eq. packagingIterations) deallocate(vx,xstart)


end
     
