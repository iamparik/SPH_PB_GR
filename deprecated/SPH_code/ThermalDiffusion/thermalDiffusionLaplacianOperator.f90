!
!  SUBROUTINE: thermalDiffusionLaplacianOperator
!
!  PURPOSE:  Subroutine that can calls different SPH operators 
!           to find interpolations of different functions
!
!   CREATED:        05/02/2022       by  PARIKSHIT BOREGOWDA
!   Last Modified:  01/28/2023       by  PARIKSHIT BOREGOWDA 
!****************************************************************************
    
subroutine thermalDiffusionLaplacianOperator(fncn,lap_fncn,oprtrTyype) 
    use config_parameter, only: SPH_dim, dataOutputPath, etype_thermal_dirichlet, &
        & etype_thermal_neumann, etype_periodic, dx_r
    use particle_data ,   only: rho, mass, itype, x, nreal,nedge,nghost,ntotal,&
            & niac, pair_i, pair_j, w, w_aa, dwdx,     &
            & eniac, epair_a, epair_s,nedge_rel_edge, edge, surf_norm,  &
            & gamma_cont, gamma_mat, gamma_discrt, del_gamma_as, &
            & gamma_mat_inv, xi1_mat_inv, xi_cont_mat_inv, etotal,etype, &
            & nperiodic,pBC_duplicate_pair, pBC_edges !, bdryVal_temp
    implicit none
    
    real(8) :: fncn(ntotal), lap_fncn(ntotal)
    integer(4) :: oprtrTyype
    integer(4) a, k, d, s
    real(8) kk(ntotal), kk_s(etotal), G(ntotal), bdryVal_temp(etotal)
    real(8),DIMENSION(:),ALLOCATABLE:: fncn_s, gamma_approx
    real(8),DIMENSION(:,:),ALLOCATABLE:: dfncn
    real(8),DIMENSION(:,:),ALLOCATABLE:: dfncn_s, x_s
    real(8),DIMENSION(:,:,:),ALLOCATABLE:: gamma_mat_approx_inv
    

    ALLOCATE(dfncn(SPH_dim,ntotal))
    dfncn=0.D0
    
    ALLOCATE(fncn_s(etotal), dfncn_s(SPH_dim,etotal))
    fncn_s=0.D0
    dfncn_s=0.D0
    G=1.D0
    
    if (Allocated(pBC_edges)) then
         call PeriodicParameter(fncn(:))
    endif
    
    ! This evaluates function value at the boundary. This can be further generalized
    do s= 1, etotal
       fncn_s(s)=fncn(nedge_rel_edge(s))
    enddo
    
    ! Initialize laplacian as zero
    lap_fncn=0.D0
    
    if(oprtrTyype .eq. 1) then
    
        !---------------------------------------------------- KCBI Formulation ------------------------------------------------   
        ! Determine function's first order derivative approximation using KCBI formulation
        call GradientKCBI( dfncn, fncn, fncn_s, gamma_cont, SPH_dim, niac, pair_i, pair_j, dwdx, eniac, epair_a, epair_s, del_gamma_as, mass, rho, itype, ntotal, etotal)

        ! For periodic condition, update dfncn
        if (Allocated(pBC_edges)) then
            do d=1,SPH_dim
                call PeriodicParameter(dfncn(d,:))
            enddo
        endif
        
        ! Determine function's laplacian using KCBI as Divergence of gradient 𝜵∙𝜵(f_a) 
        ! This evaluates function first order derivative value at the boundary. This can be further generalized for numerical integration along edge...
        do s= 1, etotal
            dfncn_s(:,s)=dfncn(:,nedge_rel_edge(s))
        enddo
        
        call DivergenceKCBI( lap_fncn, dfncn, dfncn_s, gamma_cont, SPH_dim, niac, pair_i, pair_j, dwdx, eniac, epair_a, epair_s, del_gamma_as, mass, rho, itype, ntotal, etotal)
    
    elseif(oprtrTyype .eq. 2) then
        !---------------------------------------------------- PCBI Formulation ------------------------------------------------   
        ! Determine function's first order derivative approximation using PCBI formulation
        call GradientPCBI( dfncn, fncn, fncn_s, gamma_mat_inv, SPH_dim, niac, pair_i, pair_j, dwdx, eniac, epair_a, epair_s, del_gamma_as, mass, rho, itype, ntotal, etotal)

        ! For periodic condition, update dfncn
        if (Allocated(pBC_edges)) then
            do d=1,SPH_dim
                call PeriodicParameter(dfncn(d,:))
            enddo
        endif
        
        ! Determine function's laplacian using PCBI as Divergence of gradient 𝜵∙𝜵(f_a) 
        ! This evaluates function first order derivative value at the boundary. This can be further generalized for numerical integration along edge...
        do s= 1, etotal
            dfncn_s(:,s)=dfncn(:,nedge_rel_edge(s))
        enddo
        call DivergencePCBI( lap_fncn, dfncn, dfncn_s, gamma_mat_inv, SPH_dim, niac, pair_i, pair_j, dwdx, eniac, epair_a, epair_s, del_gamma_as, mass, rho, itype, ntotal, etotal)

     
    elseif(oprtrTyype .eq. 3) then
        !---------------------------------------------------- PCBI based weakly consistent Formulation I ------------------------------------------------   
       
        ALLOCATE( gamma_mat_approx_inv(SPH_dim,SPH_dim,ntotal)) 
        gamma_mat_approx_inv=0.D0
        do a=1,ntotal
            do d=1,SPH_dim
                    gamma_mat_approx_inv(d,d,a)= 1.D0/gamma_cont(a)
            enddo
        enddo
        
        ! Determine function's first order derivative approximation using WCBI1 formulation
        call GradientPCBI( dfncn, fncn, fncn_s, gamma_mat_approx_inv, SPH_dim, niac, pair_i, pair_j, dwdx, eniac, epair_a, epair_s, del_gamma_as, mass, rho, itype, ntotal, etotal)

        ! For periodic condition, update dfncn
        if (Allocated(pBC_edges)) then
            do d=1,SPH_dim
                call PeriodicParameter(dfncn(d,:))
            enddo
        endif
        
        ! Determine function's laplacian using WCBI1 as Divergence of gradient 𝜵∙𝜵(f_a) 
        ! This evaluates function first order derivative value at the boundary. This can be further generalized for numerical integration along edge...
        do s= 1, etotal
           dfncn_s(:,s)=dfncn(:,nedge_rel_edge(s))
        enddo
        call DivergencePCBI( lap_fncn, dfncn, dfncn_s, gamma_mat_approx_inv, SPH_dim, niac, pair_i, pair_j, dwdx, eniac, epair_a, epair_s, del_gamma_as, mass, rho, itype, ntotal, etotal)

        DEALLOCATE( gamma_mat_approx_inv )
    
    elseif(oprtrTyype .eq. 4) then
    !---------------------------------------------------- PCBI based weakly consistent Formulation 2------------------------------------------------   
    
        ALLOCATE( gamma_mat_approx_inv(SPH_dim,SPH_dim,ntotal)) 
        gamma_mat_approx_inv=0.D0
        do a=1,ntotal
            do d=1,SPH_dim
                    gamma_mat_approx_inv(d,d,a)= 1.D0/gamma_discrt(a)
            enddo
        enddo
    
        ! Determine function's first order derivative approximation using WCBI2 formulation
        call GradientPCBI( dfncn, fncn, fncn_s, gamma_mat_approx_inv, SPH_dim, niac, pair_i, pair_j, dwdx, eniac, epair_a, epair_s, del_gamma_as, mass, rho, itype, ntotal, etotal)

        ! For periodic condition, update dfncn
        if (Allocated(pBC_edges)) then
            do d=1,SPH_dim
                call PeriodicParameter(dfncn(d,:))
            enddo
        endif
        
        ! Determine function's laplacian using WCBI2 as Divergence of gradient 𝜵∙𝜵(f_a) 
        ! This evaluates function first order derivative value at the boundary. This can be further generalized for numerical integration along edge...
        do s= 1, etotal
           dfncn_s(:,s)=dfncn(:,nedge_rel_edge(s))
        enddo
        call DivergencePCBI( lap_fncn, dfncn, dfncn_s, gamma_mat_approx_inv, SPH_dim, niac, pair_i, pair_j, dwdx, eniac, epair_a, epair_s, del_gamma_as, mass, rho, itype, ntotal, etotal)

        DEALLOCATE( gamma_mat_approx_inv )
    
    
    elseif(oprtrTyype .eq. 5) then
    !---------------------------------------------------- PCBI based weakly consistent Formulation 3------------------------------------------------   
        
        ALLOCATE( gamma_mat_approx_inv(SPH_dim,SPH_dim,ntotal)) 
        gamma_mat_approx_inv=0.D0
        do a=1,ntotal
            do d=1,SPH_dim
                    gamma_mat_approx_inv(d,d,a)= 1.D0/gamma_mat(d,d,a)
            enddo
        enddo
    
        ! Determine function's first order derivative approximation using WCBI3 formulation
        call GradientPCBI( dfncn, fncn, fncn_s, gamma_mat_approx_inv, SPH_dim, niac, pair_i, pair_j, dwdx, eniac, epair_a, epair_s, del_gamma_as, mass, rho, itype, ntotal, etotal)

        ! For periodic condition, update dfncn
        if (Allocated(pBC_edges)) then
            do d=1,SPH_dim
                call PeriodicParameter(dfncn(d,:))
            enddo
        endif
        
        ! Determine function's laplacian using WCBI3 as Divergence of gradient 𝜵∙𝜵(f_a) 
        ! This evaluates function first order derivative value at the boundary. This can be further generalized for numerical integration along edge...
        do s= 1, etotal
           dfncn_s(:,s)=dfncn(:,nedge_rel_edge(s))
        enddo
        
        call DivergencePCBI( lap_fncn, dfncn, dfncn_s, gamma_mat_approx_inv, SPH_dim, niac, pair_i, pair_j, dwdx, eniac, epair_a, epair_s, del_gamma_as, mass, rho, itype, ntotal, etotal)
  
        DEALLOCATE( gamma_mat_approx_inv )
    
    
    elseif(oprtrTyype .eq. 6) then
    !---------------------------------------------------- CSPM Formulation ------------------------------------------------   

        !ALLOCATE( gamma_mat_approx_inv(SPH_dim,SPH_dim,ntotal)) 
        ! Determine function's first order derivative approximation using CSPM formulation
        call GradientCSPM( dfncn, fncn, xi1_mat_inv, SPH_dim, niac, pair_i, pair_j, dwdx, mass, rho, itype, ntotal)

        ! For periodic condition, update dfncn
        if (Allocated(pBC_edges)) then
            do d=1,SPH_dim
                call PeriodicParameter(dfncn(d,:))
            enddo
        endif
        
        ! Determine function's laplacian using CSPM as Divergence of gradient 𝜵∙𝜵(f_a) 
        call DivergenceCSPM(lap_fncn, dfncn, xi1_mat_inv, SPH_dim, niac, pair_i, pair_j, dwdx, mass, rho, itype, ntotal)

    
    elseif(oprtrTyype .eq. 7) then
     !---------------------------------------------------- KCBI-CSPM Formulation ------------------------------------------------   
    
        !ALLOCATE( gamma_mat_approx_inv(SPH_dim,SPH_dim,ntotal)) 
        ! Determine function's first order derivative approximation using KCBI-CSPM formulation
        call GradientKCBI_CSPM(dfncn, fncn, xi_cont_mat_inv, SPH_dim, niac, pair_i, pair_j, dwdx, eniac, epair_a, del_gamma_as, mass, rho, itype, ntotal, etotal)


        ! For periodic condition, update dfncn
        if (Allocated(pBC_edges)) then
            do d=1,SPH_dim
                call PeriodicParameter(dfncn(d,:))
            enddo
        endif
        
        ! Determine function's laplacian using KCBI-CSPM as Divergence of gradient 𝜵∙𝜵(f_a) 
        call DivergenceKCBI_CSPM(lap_fncn, dfncn, xi_cont_mat_inv, SPH_dim, niac, pair_i, pair_j, dwdx, eniac, epair_a, del_gamma_as, mass, rho, itype, ntotal, etotal)


    elseif(oprtrTyype .eq. 8) then
        !---------------------------------------------------- BIL-USAW Formulation (using analytical gradients)------------------------------------------------   
    
        !Determine function's first order derivative approximation using PCBI formulation
        call GradientPCBI( dfncn, fncn, fncn_s, gamma_mat_inv, SPH_dim, niac, pair_i, pair_j, dwdx, eniac, epair_a, epair_s, del_gamma_as, mass, rho, itype, ntotal, etotal)

        ! For periodic condition, update dfncn
        if (Allocated(pBC_edges)) then
            do d=1,SPH_dim
                call PeriodicParameter(dfncn(d,:))
            enddo
        endif

        ! Determine function's laplacian  𝜵∙kk 𝜵(f_a) using USAW formulation    
        ! we set kk as 1 in 𝜵∙kk 𝜵(f_a)
        kk=1.D0
        
        call LaplacianSPHTraditional(lap_fncn, fncn, kk, x(1:SPH_dim,1:ntotal), &
                 & gamma_cont, SPH_dim, niac, pair_i, pair_j, dwdx, mass, rho, itype, ntotal)        
        
        do a=1,ntotal
            dfncn(:,a)= kk(a)*dfncn(:,a)
        enddo
        ! Boundary vales of fucntions are its analytical values at the boundary
        do s= 1, etotal
           dfncn_s(:,s)=dfncn(:,nedge_rel_edge(s))
        enddo

        call BILUSAW(lap_fncn, dfncn, dfncn_s, gamma_cont, SPH_dim, eniac, epair_a, epair_s, &
                        & del_gamma_as, etype, etotal, ntotal, etype_periodic, &
                        & etype_thermal_neumann,etype_thermal_dirichlet,bdryVal_temp,etotal)

        
    elseif(oprtrTyype .eq. 9) then    
    !---------------------------------------------------- BIL1 Formulation (using analytical gradients)------------------------------------------------   

        !Determine function's first order derivative approximation using PCBI formulation
        call GradientPCBI( dfncn, fncn, fncn_s, gamma_mat_inv, SPH_dim, niac, pair_i, pair_j, dwdx, eniac, epair_a, epair_s, del_gamma_as, mass, rho, itype, ntotal, etotal)

        ! For periodic condition, update dfncn
        if (Allocated(pBC_edges)) then
            do d=1,SPH_dim
                call PeriodicParameter(dfncn(d,:))
            enddo
        endif
        
        ! Determine function's laplacian  𝜵∙kk 𝜵(f_a) using BIL1 formulation    
        ! we set kk as 1 in 𝜵∙kk 𝜵(f_a)
        kk=1.D0
        
        call LaplacianSPHTraditional(lap_fncn, fncn, kk, x(1:SPH_dim,1:ntotal), &
                 & gamma_cont, SPH_dim, niac, pair_i, pair_j, dwdx, mass, rho, itype, ntotal)  

        
        call BIL1(lap_fncn, dfncn, gamma_cont, SPH_dim, eniac, epair_a, epair_s, &
                        & del_gamma_as, etype, etotal, ntotal, etype_periodic, &
                        & etype_thermal_neumann,etype_thermal_dirichlet,bdryVal_temp, 2.D0)
        
    elseif(oprtrTyype .eq. 10) then 
    !---------------------------------------------------- BIL2 Formulation (using analytical gradients)------------------------------------------------   

        !Determine function's first order derivative approximation using PCBI formulation
        call GradientPCBI( dfncn, fncn, fncn_s, gamma_mat_inv, SPH_dim, niac, pair_i, pair_j, dwdx, eniac, epair_a, epair_s, del_gamma_as, mass, rho, itype, ntotal, etotal)

        ! For periodic condition, update dfncn
        if (Allocated(pBC_edges)) then
            do d=1,SPH_dim
                call PeriodicParameter(dfncn(d,:))
            enddo
        endif
        
        ! Determine function's laplacian  𝜵∙kk 𝜵(f_a) using BIL2 formulation    
        ! we set kk as 1 in 𝜵∙kk 𝜵(f_a)
        kk=1.D0
        
        call LaplacianSPHTraditional(lap_fncn, fncn, kk, x(1:SPH_dim,1:ntotal), &
                 & gamma_cont, SPH_dim, niac, pair_i, pair_j, dwdx, mass, rho, itype, ntotal)  

        ! Boundary values of fucntions are its analytical values at the boundary
        do s= 1, etotal
           dfncn_s(:,s)=dfncn(:,nedge_rel_edge(s))
        enddo
        
        call BIL2(lap_fncn, dfncn_s, gamma_cont, SPH_dim, eniac, epair_a, epair_s, &
                        & del_gamma_as, etype, etotal, ntotal, etype_periodic, &
                        & etype_thermal_neumann,etype_thermal_dirichlet,bdryVal_temp,etotal, 2.D0)

    elseif(oprtrTyype .eq. 11) then 
    
    !---------------------------------------------------- BIL-Macia Formulation------------------------------------------------   

        ALLOCATE(x_s(SPH_dim,etotal))       
        x_s=0.D0
        
        ! Determine function's laplacian  𝜵∙ 𝜵(f_a) using BIL-Macia formulation 
        ! Boundary values of fucntions are its analytical values at the boundary
        do s= 1, etotal
           x_s(:,s)=x(:,nedge_rel_edge(s))
        enddo
        
        ! Determine function's laplacian  𝜵∙kk 𝜵(f_a) using BIL formulation    
        ! we set kk as 1 in 𝜵∙kk 𝜵(f_a)
        kk=1.D0
        
        call LaplacianSPHTraditional(lap_fncn, fncn, kk, x(1:SPH_dim,1:ntotal), &
                     & gamma_cont, SPH_dim, niac, pair_i, pair_j, dwdx, mass, rho, itype, ntotal)  

        ! Boundary vales of fucntions are its analytical values at the boundary
        do s= 1, etotal
           kk_s(s)=1.D0
        enddo
    
        call BILMacia(lap_fncn, fncn, fncn_s, x(1:SPH_dim,1:ntotal), x_s, gamma_cont, SPH_dim, eniac,&
                & epair_a, epair_s, del_gamma_as, etype, etotal, ntotal, kk_s)
        
        deallocate(x_s)
  
    elseif(oprtrTyype .eq. 12) then 
    
    !---------------------------------------------------- BIL-Macia Formulation with 2X boundary------------------------------------------------   

        ALLOCATE(x_s(SPH_dim,etotal))       
        x_s=0.D0
        
        ! Determine function's laplacian  𝜵∙ 𝜵(f_a) using BIL-Macia formulation 
        ! Boundary values of fucntions are its analytical values at the boundary
        do s= 1, etotal
           x_s(:,s)=x(:,nedge_rel_edge(s))
        enddo
  
        ! Determine function's laplacian  𝜵∙kk 𝜵(f_a) using BIL formulation    
        ! we set kk as 1 in 𝜵∙kk 𝜵(f_a)
        kk=1.D0
        
        call LaplacianSPHTraditional(lap_fncn, fncn, kk, x(1:SPH_dim,1:ntotal), &
                     & gamma_cont, SPH_dim, niac, pair_i, pair_j, dwdx, mass, rho, itype, ntotal)  

        ! Boundary vales of fucntions are its analytical values at the boundary
        do s= 1, etotal
           kk_s(s)=2.D0
        enddo
    
        call BILMacia(lap_fncn, fncn, fncn_s, x(1:SPH_dim,1:ntotal), x_s, gamma_cont, SPH_dim, eniac,&
                & epair_a, epair_s, del_gamma_as, etype, etotal, ntotal, kk_s)
        
        deallocate(x_s)
    
        elseif(oprtrTyype .eq. 13) then    
    !---------------------------------------------------- BIL1 Formulation (using discrete gamma for laplacian)------------------------------------------------   
        !Determine function's first order derivative approximation using PCBI formulation
        call GradientPCBI( dfncn, fncn, fncn_s, gamma_mat_inv, SPH_dim, niac, pair_i, pair_j, dwdx, eniac, epair_a, epair_s, del_gamma_as, mass, rho, itype, ntotal, etotal)

        ! For periodic condition, update dfncn
        if (Allocated(pBC_edges)) then
            do d=1,SPH_dim
                call PeriodicParameter(dfncn(d,:))
            enddo
        endif
        
        ! Determine function's laplacian  𝜵∙kk 𝜵(f_a) using BIL1 formulation    
        ! we set kk as 1 in 𝜵∙kk 𝜵(f_a)
        
        call LaplacianSPHTraditional(lap_fncn, fncn, G, x(1:SPH_dim,1:ntotal), &
                 & gamma_discrt, SPH_dim, niac, pair_i, pair_j, dwdx, mass, rho, itype, ntotal)  

        
        call BIL1(lap_fncn, dfncn, gamma_discrt, SPH_dim, eniac, epair_a, epair_s, &
                        & del_gamma_as, etype, etotal, ntotal, etype_periodic, &
                        & etype_thermal_neumann,etype_thermal_dirichlet,bdryVal_temp, 2.D0)
        
     elseif(oprtrTyype .eq. 14) then    
    !---------------------------------------------------- BIL2 Formulation (using discrete gamma for laplacian) ------------------------------------------------   

              !Determine function's first order derivative approximation using PCBI formulation
        call GradientPCBI( dfncn, fncn, fncn_s, gamma_mat_inv, SPH_dim, niac, pair_i, pair_j, dwdx, eniac, epair_a, epair_s, del_gamma_as, mass, rho, itype, ntotal, etotal)

        ! For periodic condition, update dfncn
        if (Allocated(pBC_edges)) then
            do d=1,SPH_dim
                call PeriodicParameter(dfncn(d,:))
            enddo
        endif
        
        ! Determine function's laplacian  𝜵∙kk 𝜵(f_a) using BIL2 formulation    
        ! we set kk as 1 in 𝜵∙kk 𝜵(f_a)

        call LaplacianSPHTraditional(lap_fncn, fncn, G, x(1:SPH_dim,1:ntotal), &
                 & gamma_discrt, SPH_dim, niac, pair_i, pair_j, dwdx, mass, rho, itype, ntotal)  

        ! Boundary values of fucntions are its analytical values at the boundary
        do s= 1, etotal
           dfncn_s(:,s)=dfncn(:,nedge_rel_edge(s))
        enddo
        
        call BIL2(lap_fncn, dfncn_s, gamma_discrt, SPH_dim, eniac, epair_a, epair_s, &
                        & del_gamma_as, etype, etotal, ntotal, etype_periodic, &
                        & etype_thermal_neumann,etype_thermal_dirichlet,bdryVal_temp,etotal, 2.D0)

    elseif(oprtrTyype .eq. 15) then    
    !---------------------------------------------------- BIL1 Corrected Formulation------------------------------------------------   

        !Determine function's first order derivative approximation using PCBI formulation
        call GradientPCBI( dfncn, fncn, fncn_s, gamma_mat_inv, SPH_dim, niac, pair_i, pair_j, dwdx, eniac, epair_a, epair_s, del_gamma_as, mass, rho, itype, ntotal, etotal)

        ! For periodic condition, update dfncn
        if (Allocated(pBC_edges)) then
            do d=1,SPH_dim
                call PeriodicParameter(dfncn(d,:))
            enddo
        endif

        ! Determine function's laplacian  𝜵∙ 𝜵(f_a) using BIL1 formulation
        call LaplacianSPHTraditionalCorrected(lap_fncn, fncn, G, x(1:SPH_dim,1:ntotal), &
            & gamma_mat_inv, SPH_dim, niac, pair_i, pair_j, dwdx, mass, rho, itype, ntotal)
        
        call BIL1Corrected(lap_fncn, dfncn, gamma_mat_inv, SPH_dim, eniac,&
            & epair_a, epair_s, del_gamma_as, etype, etotal, ntotal, 2.D0)

        
     elseif(oprtrTyype .eq. 16) then    
    !----------------------------------------------------  BIL2 Corrected Formulation ------------------------------------------------   

      !Determine function's first order derivative approximation using PCBI formulation
        call GradientPCBI( dfncn, fncn, fncn_s, gamma_mat_inv, SPH_dim, niac, pair_i, pair_j,&
            & dwdx, eniac, epair_a, epair_s, del_gamma_as, mass, rho, itype, ntotal, etotal)

        ! For periodic condition, update dfncn
        if (Allocated(pBC_edges)) then
            do d=1,SPH_dim
                call PeriodicParameter(dfncn(d,:))
            enddo
        endif

        ! Determine function's laplacian  𝜵∙ 𝜵(f_a) using BIL2 formulation
        call LaplacianSPHTraditionalCorrected(lap_fncn, fncn, G, x(1:SPH_dim,1:ntotal), &
            & gamma_mat_inv, SPH_dim, niac, pair_i, pair_j, dwdx, mass, rho, itype, ntotal)
        ! Boundary values of fucntions are its analytical values at the boundary
        do s= 1, etotal
           dfncn_s(:,s)=dfncn(:,nedge_rel_edge(s))
        enddo
        
        call BIL2Corrected(lap_fncn, dfncn_s, gamma_mat_inv, SPH_dim, eniac, &
            & epair_a, epair_s, del_gamma_as, etype, etotal, ntotal, 2.D0)   

        
        elseif(oprtrTyype .eq. 17) then    
    !---------------------------------------------------- BIL1 Corrected Formulation inside only------------------------------------------------   

        !Determine function's first order derivative approximation using PCBI formulation
        call GradientPCBI( dfncn, fncn, fncn_s, gamma_mat_inv, SPH_dim, niac, pair_i, pair_j, &
            & dwdx, eniac, epair_a, epair_s, del_gamma_as, mass, rho, itype, ntotal, etotal)

        ! For periodic condition, update dfncn
        if (Allocated(pBC_edges)) then
            do d=1,SPH_dim
                call PeriodicParameter(dfncn(d,:))
            enddo
        endif
        
        ALLOCATE( gamma_mat_approx_inv(SPH_dim,SPH_dim,ntotal)) 
        gamma_mat_approx_inv=0.D0
        do a=1,ntotal
            if((gamma_cont(a) .lt. 1.D0) .and. (gamma_cont(a) .gt. 0.D0)) then
                do d=1,SPH_dim
                    gamma_mat_approx_inv(d,d,a)= 1.D0/gamma_cont(a)
                enddo
            else
                gamma_mat_approx_inv(:,:,a)=gamma_mat_inv(:,:,a)
            endif            
        enddo

        ! Determine function's laplacian  𝜵∙ 𝜵(f_a) using BIL1 formulation
        call LaplacianSPHTraditionalCorrected(lap_fncn, fncn, G, x(1:SPH_dim,1:ntotal),&
            & gamma_mat_approx_inv, SPH_dim, niac, pair_i, pair_j, dwdx, mass, rho, itype, ntotal)
        
        call BIL1Corrected(lap_fncn, dfncn, gamma_mat_approx_inv, SPH_dim, eniac,&
            & epair_a, epair_s, del_gamma_as, etype, etotal, ntotal, 2.D0)

        deallocate(gamma_mat_approx_inv)
        
     elseif(oprtrTyype .eq. 18) then    
    !---------------------------------------------------- BIL2 Corrected Formulation inside only------------------------------------------------   

      !Determine function's first order derivative approximation using PCBI formulation
        call GradientPCBI( dfncn, fncn, fncn_s, gamma_mat_inv, SPH_dim, niac, pair_i, pair_j,&
            & dwdx, eniac, epair_a, epair_s, del_gamma_as, mass, rho, itype, ntotal, etotal)

        ! For periodic condition, update dfncn
        if (Allocated(pBC_edges)) then
            do d=1,SPH_dim
                call PeriodicParameter(dfncn(d,:))
            enddo
        endif
        
        ALLOCATE( gamma_mat_approx_inv(SPH_dim,SPH_dim,ntotal)) 
        gamma_mat_approx_inv=gamma_mat_inv
        do a=1,ntotal
            if((gamma_cont(a) .lt. 1.D0) .and. (gamma_cont(a) .gt. 0.D0)) then
                do d=1,SPH_dim
                        gamma_mat_approx_inv(d,d,a)= 1.D0/gamma_cont(a)
                enddo
            endif
        enddo

        ! Determine function's laplacian  𝜵∙ 𝜵(f_a) using BIL2 formulation
        call LaplacianSPHTraditionalCorrected(lap_fncn, fncn, G, x(1:SPH_dim,1:ntotal), &
            & gamma_mat_approx_inv, SPH_dim, niac, pair_i, pair_j, dwdx, mass, rho, itype, ntotal)
        ! Boundary values of fucntions are its analytical values at the boundary
        do s= 1, etotal
           dfncn_s(:,s)=dfncn(:,nedge_rel_edge(s))
        enddo
        
        call BIL2Corrected(lap_fncn, dfncn_s, gamma_mat_approx_inv, SPH_dim, eniac,&
            & epair_a, epair_s, del_gamma_as, etype, etotal, ntotal, 2.D0)   
    
        deallocate(gamma_mat_approx_inv)
        
    elseif(oprtrTyype .eq. 19) then    
    !---------------------------------------------------- BIL1 Corrected Formulation, but gamma_cont for boundary integral ------------------------------------------------   

        !Determine function's first order derivative approximation using PCBI formulation
        call GradientPCBI( dfncn, fncn, fncn_s, gamma_mat_inv, SPH_dim, niac, pair_i, pair_j,&
            & dwdx, eniac, epair_a, epair_s, del_gamma_as, mass, rho, itype, ntotal, etotal)

        ! For periodic condition, update dfncn
        if (Allocated(pBC_edges)) then
            do d=1,SPH_dim
                call PeriodicParameter(dfncn(d,:))
            enddo
        endif
        
        ! Determine function's laplacian  𝜵∙ 𝜵(f_a) using BIL1 formulation
        call LaplacianSPHTraditionalCorrected(lap_fncn, fncn, G, x(1:SPH_dim,1:ntotal),&
            & gamma_mat_inv, SPH_dim, niac, pair_i, pair_j, dwdx, mass, rho, itype, ntotal)
        
        call BIL1(lap_fncn, dfncn, gamma_cont, SPH_dim, eniac, epair_a, epair_s, &
                        & del_gamma_as, etype, etotal, ntotal, etype_periodic, &
                        & etype_thermal_neumann,etype_thermal_dirichlet,bdryVal_temp, 2.D0)
        
     elseif(oprtrTyype .eq. 20) then    
    !---------------------------------------------------- BIL2 Corrected, but gamma_cont for boundary integral ------------------------------------------------   

      !Determine function's first order derivative approximation using PCBI formulation
        call GradientPCBI( dfncn, fncn, fncn_s, gamma_mat_inv, SPH_dim, niac,&
            & pair_i, pair_j, dwdx, eniac, epair_a, epair_s, del_gamma_as, mass, rho, itype, ntotal, etotal)

        ! For periodic condition, update dfncn
        if (Allocated(pBC_edges)) then
            do d=1,SPH_dim
                call PeriodicParameter(dfncn(d,:))
            enddo
        endif
        
        ! Determine function's laplacian  𝜵∙ 𝜵(f_a) using BIL2 formulation
        call LaplacianSPHTraditionalCorrected(lap_fncn, fncn, G, x(1:SPH_dim,1:ntotal),&
            & gamma_mat_inv, SPH_dim, niac, pair_i, pair_j, dwdx, mass, rho, itype, ntotal)
      
        ! Boundary values of fucntions are its analytical values at the boundary
        do s= 1, etotal
           dfncn_s(:,s)=dfncn(:,nedge_rel_edge(s))
        enddo
        
        call BIL2(lap_fncn, dfncn_s, gamma_cont, SPH_dim, eniac, epair_a, epair_s, &
                        & del_gamma_as, etype, etotal, ntotal, etype_periodic, &
                        & etype_thermal_neumann,etype_thermal_dirichlet,bdryVal_temp,etotal, 2.D0)

        
    elseif(oprtrTyype .eq. 21) then    
    !---------------------------------------------------- BIL1 Corrected Formulation with Xi------------------------------------------------   

        !Determine function's first order derivative approximation using PCBI formulation
        call GradientPCBI( dfncn, fncn, fncn_s, gamma_mat_inv, SPH_dim, niac, pair_i, pair_j,&
            & dwdx, eniac, epair_a, epair_s, del_gamma_as, mass, rho, itype, ntotal, etotal)

        ! For periodic condition, update dfncn
        if (Allocated(pBC_edges)) then
            do d=1,SPH_dim
                call PeriodicParameter(dfncn(d,:))
            enddo
        endif

        ! Determine function's laplacian  𝜵∙ 𝜵(f_a) using BIL1 formulation
        call LaplacianSPHTraditionalCorrected(lap_fncn, fncn, G, x(1:SPH_dim,1:ntotal),&
            & xi1_mat_inv, SPH_dim, niac, pair_i, pair_j, dwdx, mass, rho, itype, ntotal)
        
        call BIL1(lap_fncn, dfncn, gamma_cont, SPH_dim, eniac, epair_a, epair_s, &
                        & del_gamma_as, etype, etotal, ntotal, etype_periodic, &
                        & etype_thermal_neumann,etype_thermal_dirichlet,bdryVal_temp, 2.D0)
        
     elseif(oprtrTyype .eq. 22) then    
    !---------------------------------------------------- BIL2 Corrected Formulation with Xi------------------------------------------------   

      !Determine function's first order derivative approximation using PCBI formulation
        call GradientPCBI( dfncn, fncn, fncn_s, gamma_mat_inv, SPH_dim, niac,&
            & pair_i, pair_j, dwdx, eniac, epair_a, epair_s, del_gamma_as, mass, rho, itype, ntotal, etotal)

        ! For periodic condition, update dfncn
        if (Allocated(pBC_edges)) then
            do d=1,SPH_dim
                call PeriodicParameter(dfncn(d,:))
            enddo
        endif


        ! Determine function's laplacian  𝜵∙ 𝜵(f_a) using BIL2 formulation
        call LaplacianSPHTraditionalCorrected(lap_fncn, fncn, G, x(1:SPH_dim,1:ntotal),&
            & xi1_mat_inv, SPH_dim, niac, pair_i, pair_j, dwdx, mass, rho, itype, ntotal)
      
        ! Boundary values of fucntions are its analytical values at the boundary
        do s= 1, etotal
           dfncn_s(:,s)=dfncn(:,nedge_rel_edge(s))
        enddo
        
        call BIL2(lap_fncn, dfncn_s, gamma_cont, SPH_dim, eniac, epair_a, epair_s, &
                        & del_gamma_as, etype, etotal, ntotal, etype_periodic, &
                        & etype_thermal_neumann,etype_thermal_dirichlet,bdryVal_temp,etotal, 2.D0)
    
    elseif(oprtrTyype .eq. 23) then    
    !---------------------------------------------------- BIL1 Corrected Formulation Xi With γ_a------------------------------------------------   

        !Determine function's first order derivative approximation using PCBI formulation
        call GradientPCBI( dfncn, fncn, fncn_s, gamma_mat_inv, SPH_dim, niac, pair_i, pair_j,&
            & dwdx, eniac, epair_a, epair_s, del_gamma_as, mass, rho, itype, ntotal, etotal)

        ! For periodic condition, update dfncn
        if (Allocated(pBC_edges)) then
            do d=1,SPH_dim
                call PeriodicParameter(dfncn(d,:))
            enddo
        endif
        
        ALLOCATE( gamma_mat_approx_inv(SPH_dim,SPH_dim,ntotal)) 
        gamma_mat_approx_inv=0.D0
        do a=1,ntotal
            gamma_mat_approx_inv(:,:,a)=xi1_mat_inv(:,:,a)/gamma_cont(a)
        enddo

        ! Determine function's laplacian  𝜵∙ 𝜵(f_a) using BIL1 formulation
        call LaplacianSPHTraditionalCorrected(lap_fncn, fncn, G, x(1:SPH_dim,1:ntotal),&
            & gamma_mat_approx_inv, SPH_dim, niac, pair_i, pair_j, dwdx, mass, rho, itype, ntotal)
        
        call BIL1(lap_fncn, dfncn, gamma_cont, SPH_dim, eniac, epair_a, epair_s, &
                        & del_gamma_as, etype, etotal, ntotal, etype_periodic, &
                        & etype_thermal_neumann,etype_thermal_dirichlet,bdryVal_temp, 2.D0)

        deallocate(gamma_mat_approx_inv)
        
     elseif(oprtrTyype .eq. 24) then    
    !---------------------------------------------------- BIL2 Corrected Formulation Xi With γ_a------------------------------------------------   

      !Determine function's first order derivative approximation using PCBI formulation
        call GradientPCBI( dfncn, fncn, fncn_s, gamma_mat_inv, SPH_dim, niac,&
            & pair_i, pair_j, dwdx, eniac, epair_a, epair_s, del_gamma_as, mass, rho, itype, ntotal, etotal)

        ! For periodic condition, update dfncn
        if (Allocated(pBC_edges)) then
            do d=1,SPH_dim
                call PeriodicParameter(dfncn(d,:))
            enddo
        endif
        
        ALLOCATE( gamma_mat_approx_inv(SPH_dim,SPH_dim,ntotal)) 
        gamma_mat_approx_inv=0.D0
        do a=1,ntotal
            gamma_mat_approx_inv(:,:,a)=xi1_mat_inv(:,:,a)/gamma_cont(a)
        enddo

        ! Determine function's laplacian  𝜵∙ 𝜵(f_a) using BIL2 formulation
        call LaplacianSPHTraditionalCorrected(lap_fncn, fncn, G, x(1:SPH_dim,1:ntotal),&
            & gamma_mat_approx_inv, SPH_dim, niac, pair_i, pair_j, dwdx, mass, rho, itype, ntotal)
      
        ! Boundary values of fucntions are its analytical values at the boundary
        do s= 1, etotal
           dfncn_s(:,s)=dfncn(:,nedge_rel_edge(s))
        enddo
        
        call BIL2(lap_fncn, dfncn_s, gamma_cont, SPH_dim, eniac, epair_a, epair_s, &
                        & del_gamma_as, etype, etotal, ntotal, etype_periodic, &
                        & etype_thermal_neumann,etype_thermal_dirichlet,bdryVal_temp,etotal, 2.D0)
    
        deallocate(gamma_mat_approx_inv)    
    
    elseif(oprtrTyype .eq. 25) then
        !---------------------------------------------------- BIL-USAW Formulation corrected------------------------------------------------   
    
        ! Determine function's first order derivative approximation using PCBI formulation
         call GradientPCBI( dfncn, fncn, fncn_s, gamma_mat_inv, SPH_dim, niac,&
            & pair_i, pair_j, dwdx, eniac, epair_a, epair_s, del_gamma_as, mass, rho, itype, ntotal, etotal)
        
        ! For periodic condition, update dfncn
        if (Allocated(pBC_edges)) then
            do d=1,SPH_dim
                call PeriodicParameter(dfncn(d,:))
            enddo
        endif
        
        ! Determine function's laplacian  𝜵∙kk 𝜵(f_a) using USAW formulation    
        ! we set kk as 1 in 𝜵∙kk 𝜵(f_a)
        kk=1.D0
        

        call LaplacianSPHTraditionalCorrected(lap_fncn, fncn, kk, x(1:SPH_dim,1:ntotal),&
            & gamma_mat_inv, SPH_dim, niac, pair_i, pair_j, dwdx, mass, rho, itype, ntotal)
        
        do a=1,ntotal
            dfncn(:,a)= kk(a)*dfncn(:,a)
        enddo
        ! Boundary vales of fucntions are its analytical values at the boundary
        do s= 1, etotal
           dfncn_s(:,s)=dfncn(:,nedge_rel_edge(s))
        enddo
        
        call BIL1Corrected(lap_fncn, dfncn, gamma_mat_inv, SPH_dim, eniac,&
            & epair_a, epair_s, del_gamma_as, etype, etotal, ntotal, 1.D0)
        call BIL2Corrected(lap_fncn, dfncn_s, gamma_mat_inv, SPH_dim, eniac,&
            & epair_a, epair_s, del_gamma_as, etype, etotal, ntotal, 1.D0)        
        
        
    elseif(oprtrTyype .eq. 26) then
        !---------------------------------------------------- BIL-USAW Formulation corrected inside only------------------------------------------------   
    
        ! Determine function's first order derivative approximation using PCBI formulation
         call GradientPCBI( dfncn, fncn, fncn_s, gamma_mat_inv, SPH_dim, niac,&
            & pair_i, pair_j, dwdx, eniac, epair_a, epair_s, del_gamma_as, mass, rho, itype, ntotal, etotal)
        
        ! For periodic condition, update dfncn
        if (Allocated(pBC_edges)) then
            do d=1,SPH_dim
                call PeriodicParameter(dfncn(d,:))
            enddo
        endif
        
        
        ALLOCATE( gamma_mat_approx_inv(SPH_dim,SPH_dim,ntotal)) 
        gamma_mat_approx_inv=0.D0
        do a=1,ntotal
            if((gamma_cont(a) .lt. 1.D0) .and. (gamma_cont(a) .gt. 0.D0)) then
                do d=1,SPH_dim
                    gamma_mat_approx_inv(d,d,a)= 1.D0/gamma_cont(a)
                enddo
            else
                gamma_mat_approx_inv(:,:,a)=gamma_mat_inv(:,:,a)
            endif            
        enddo 
        
        ! Determine function's laplacian  𝜵∙kk 𝜵(f_a) using USAW formulation    
        ! we set kk as 1 in 𝜵∙kk 𝜵(f_a)
        kk=1.D0
        
        call LaplacianSPHTraditionalCorrected(lap_fncn, fncn, kk, x(1:SPH_dim,1:ntotal),&
            & gamma_mat_approx_inv, SPH_dim, niac, pair_i, pair_j, dwdx, mass, rho, itype, ntotal)
        
        do a=1,ntotal
            dfncn(:,a)= kk(a)*dfncn(:,a)
        enddo
        ! Boundary vales of fucntions are its analytical values at the boundary
        do s= 1, etotal
           dfncn_s(:,s)=dfncn(:,nedge_rel_edge(s))
        enddo
        
        call BIL1(lap_fncn, dfncn, gamma_cont, SPH_dim, eniac, epair_a, epair_s, &
                        & del_gamma_as, etype, etotal, ntotal, etype_periodic, &
                        & etype_thermal_neumann,etype_thermal_dirichlet,bdryVal_temp, 1.D0)
        call BIL2(lap_fncn, dfncn_s, gamma_cont, SPH_dim, eniac, epair_a, epair_s, &
                        & del_gamma_as, etype, etotal, ntotal, etype_periodic, &
                        & etype_thermal_neumann,etype_thermal_dirichlet,bdryVal_temp,etotal, 1.D0)
        deallocate(gamma_mat_approx_inv)
        
    elseif(oprtrTyype .eq. 27) then
        !---------------------------------------------------- BIL-USAW Formulation corrected, but gamma_cont for boundary integral------------------------------------------------   
    
        ! Determine function's first order derivative approximation using PCBI formulation
         call GradientPCBI( dfncn, fncn, fncn_s, gamma_mat_inv, SPH_dim, niac,&
            & pair_i, pair_j, dwdx, eniac, epair_a, epair_s, del_gamma_as, mass, rho, itype, ntotal, etotal)
        
        ! For periodic condition, update dfncn
        if (Allocated(pBC_edges)) then
            do d=1,SPH_dim
                call PeriodicParameter(dfncn(d,:))
            enddo
        endif
        
        
        ! Determine function's laplacian  𝜵∙kk 𝜵(f_a) using USAW formulation    
        ! we set kk as 1 in 𝜵∙kk 𝜵(f_a)
        kk=1.D0
        
        call LaplacianSPHTraditionalCorrected(lap_fncn, fncn, kk, x(1:SPH_dim,1:ntotal),&
            & gamma_mat_inv, SPH_dim, niac, pair_i, pair_j, dwdx, mass, rho, itype, ntotal)
        
        do a=1,ntotal
            dfncn(:,a)= kk(a)*dfncn(:,a)
        enddo
        ! Boundary vales of fucntions are its analytical values at the boundary
        do s= 1, etotal
           dfncn_s(:,s)=dfncn(:,nedge_rel_edge(s))
        enddo

        call BIL1(lap_fncn, dfncn, gamma_cont, SPH_dim, eniac, epair_a, epair_s, &
                        & del_gamma_as, etype, etotal, ntotal, etype_periodic, &
                        & etype_thermal_neumann,etype_thermal_dirichlet,bdryVal_temp, 1.D0)
        call BIL2(lap_fncn, dfncn_s, gamma_cont, SPH_dim, eniac, epair_a, epair_s, &
                        & del_gamma_as, etype, etotal, ntotal, etype_periodic, &
                        & etype_thermal_neumann,etype_thermal_dirichlet,bdryVal_temp,etotal, 1.D0)
        
     elseif(oprtrTyype .eq. 28) then
        !---------------------------------------------------- BIL-USAW Formulation corrected with gamma_discrete ------------------------------------------------   
    
        ! Determine function's first order derivative approximation using PCBI formulation
         call GradientPCBI( dfncn, fncn, fncn_s, gamma_mat_inv, SPH_dim, niac,&
            & pair_i, pair_j, dwdx, eniac, epair_a, epair_s, del_gamma_as, mass, rho, itype, ntotal, etotal)
        
        ! For periodic condition, update dfncn
        if (Allocated(pBC_edges)) then
            do d=1,SPH_dim
                call PeriodicParameter(dfncn(d,:))
            enddo
        endif
        
        
        
        ! Determine function's laplacian  𝜵∙kk 𝜵(f_a) using USAW formulation    
        ! we set kk as 1 in 𝜵∙kk 𝜵(f_a)
        kk=1.D0
        
        call LaplacianSPHTraditional(lap_fncn, fncn, kk, x(1:SPH_dim,1:ntotal), &
                 & gamma_discrt, SPH_dim, niac, pair_i, pair_j, dwdx, mass, rho, itype, ntotal) 

        
        do a=1,ntotal
            dfncn(:,a)= kk(a)*dfncn(:,a)
        enddo
        ! Boundary vales of fucntions are its analytical values at the boundary
        do s= 1, etotal
           dfncn_s(:,s)=dfncn(:,nedge_rel_edge(s))
        enddo
        
        call BIL1(lap_fncn, dfncn, gamma_discrt, SPH_dim, eniac, epair_a, epair_s, &
                        & del_gamma_as, etype, etotal, ntotal, etype_periodic, &
                        & etype_thermal_neumann,etype_thermal_dirichlet,bdryVal_temp, 1.D0)
        call BIL2(lap_fncn, dfncn_s, gamma_discrt, SPH_dim, eniac, epair_a, epair_s, &
                        & del_gamma_as, etype, etotal, ntotal, etype_periodic, &
                        & etype_thermal_neumann,etype_thermal_dirichlet,bdryVal_temp,etotal, 1.D0)    
        
        
    elseif(oprtrTyype .eq. 29) then
        !---------------------------------------------------- BIL-USAW Formulation corrected inside only with gamma_discrete------------------------------------------------   
    
        ! Determine function's first order derivative approximation using PCBI formulation
         call GradientPCBI( dfncn, fncn, fncn_s, gamma_mat_inv, SPH_dim, niac,&
            & pair_i, pair_j, dwdx, eniac, epair_a, epair_s, del_gamma_as, mass, rho, itype, ntotal, etotal)
        
        ! For periodic condition, update dfncn
        if (Allocated(pBC_edges)) then
            do d=1,SPH_dim
                call PeriodicParameter(dfncn(d,:))
            enddo
        endif
        
        ALLOCATE( gamma_approx(ntotal)) 
        gamma_approx=0.D0
        do a=1,ntotal
            if((gamma_cont(a) .lt. 1.D0) .and. (gamma_cont(a) .gt. 0.D0)) then
                do d=1,SPH_dim
                    gamma_approx(a)= gamma_cont(a)
                enddo
            else
                gamma_approx(a)= gamma_discrt(a)
            endif            
        enddo 
        
        ! Determine function's laplacian  𝜵∙kk 𝜵(f_a) using USAW formulation    
        ! we set kk as 1 in 𝜵∙kk 𝜵(f_a)
        kk=1.D0
        
        call LaplacianSPHTraditional(lap_fncn, fncn, kk, x(1:SPH_dim,1:ntotal), &
                 & gamma_approx, SPH_dim, niac, pair_i, pair_j, dwdx, mass, rho, itype, ntotal)
        
        do a=1,ntotal
            dfncn(:,a)= kk(a)*dfncn(:,a)
        enddo
        ! Boundary vales of fucntions are its analytical values at the boundary
        do s= 1, etotal
           dfncn_s(:,s)=dfncn(:,nedge_rel_edge(s))
        enddo
        
        call BIL1(lap_fncn, dfncn, gamma_cont, SPH_dim, eniac, epair_a, epair_s, &
                        & del_gamma_as, etype, etotal, ntotal, etype_periodic, &
                        & etype_thermal_neumann,etype_thermal_dirichlet,bdryVal_temp, 1.D0)
        call BIL2(lap_fncn, dfncn_s, gamma_cont, SPH_dim, eniac, epair_a, epair_s, &
                        & del_gamma_as, etype, etotal, ntotal, etype_periodic, &
                        & etype_thermal_neumann,etype_thermal_dirichlet,bdryVal_temp,etotal, 1.D0)
        
        deallocate(gamma_approx)

        
    elseif(oprtrTyype .eq. 30) then
        !---------------------------------------------------- BIL-USAW Formulation corrected with gamma_discrete, but gamma_cont for boundary integral------------------------------------------------   
    
        ! Determine function's first order derivative approximation using PCBI formulation
         call GradientPCBI( dfncn, fncn, fncn_s, gamma_mat_inv, SPH_dim, niac,&
            & pair_i, pair_j, dwdx, eniac, epair_a, epair_s, del_gamma_as, mass, rho, itype, ntotal, etotal)
        
        ! For periodic condition, update dfncn
        if (Allocated(pBC_edges)) then
            do d=1,SPH_dim
                call PeriodicParameter(dfncn(d,:))
            enddo
        endif
        
        
        ! Determine function's laplacian  𝜵∙kk 𝜵(f_a) using USAW formulation    
        ! we set kk as 1 in 𝜵∙kk 𝜵(f_a)
        kk=1.D0
        
        call LaplacianSPHTraditional(lap_fncn, fncn, kk, x(1:SPH_dim,1:ntotal), &
                 & gamma_discrt, SPH_dim, niac, pair_i, pair_j, dwdx, mass, rho, itype, ntotal)
        
        do a=1,ntotal
            dfncn(:,a)= kk(a)*dfncn(:,a)
        enddo
        ! Boundary vales of fucntions are its analytical values at the boundary
        do s= 1, etotal
           dfncn_s(:,s)=dfncn(:,nedge_rel_edge(s))
        enddo

        call BIL1(lap_fncn, dfncn, gamma_cont, SPH_dim, eniac, epair_a, epair_s, &
                        & del_gamma_as, etype, etotal, ntotal, etype_periodic, &
                        & etype_thermal_neumann,etype_thermal_dirichlet,bdryVal_temp, 1.D0)
        call BIL2(lap_fncn, dfncn_s, gamma_cont, SPH_dim, eniac, epair_a, epair_s, &
                        & del_gamma_as, etype, etotal, ntotal, etype_periodic, &
                        & etype_thermal_neumann,etype_thermal_dirichlet,bdryVal_temp,etotal, 1.D0)
    
    elseif(oprtrTyype .eq. 31) then
        !---------------------------------------------------- BIL-MEA Formulation (using analytical gradients)------------------------------------------------   
    
        !Determine function's first order derivative approximation using PCBI formulation
        call GradientPCBI( dfncn, fncn, fncn_s, gamma_mat_inv, SPH_dim, niac, pair_i, pair_j, dwdx, eniac, epair_a, epair_s, del_gamma_as, mass, rho, itype, ntotal, etotal)

        ! For periodic condition, update dfncn
        if (Allocated(pBC_edges)) then
            do d=1,SPH_dim
                call PeriodicParameter(dfncn(d,:))
            enddo
        endif

        ! Determine function's laplacian  𝜵∙kk 𝜵(f_a) using USAW formulation    
        ! we set kk as 1 in 𝜵∙kk 𝜵(f_a)
        kk=1.D0
        
        call LaplacianSPHTraditional(lap_fncn, fncn, kk, x(1:SPH_dim,1:ntotal), &
                 & gamma_cont, SPH_dim, niac, pair_i, pair_j, dwdx, mass, rho, itype, ntotal)        
        
        ALLOCATE(x_s(SPH_dim,etotal))
        do s= 1, etotal
           x_s(:,s)=x(:,nedge_rel_edge(s))
           kk_s(s)=1.D0
        enddo

        call BILMEA(lap_fncn, fncn, fncn_s, x(1:SPH_dim,1:ntotal), x_s, gamma_cont, SPH_dim, eniac, epair_a, epair_s, &
                &  del_gamma_as, etype, etotal, ntotal,  etype_periodic, etype_thermal_neumann,etype_thermal_dirichlet, &
                & bdryVal_temp,etotal, surf_norm, kk, kk_s)
      
        deallocate(x_s)
        
    elseif(oprtrTyype .eq. 32) then
        !---------------------------------------------------- BIL-MEA Formulation (using disrete gamma)------------------------------------------------   
    
        !Determine function's first order derivative approximation using PCBI formulation
        call GradientPCBI( dfncn, fncn, fncn_s, gamma_mat_inv, SPH_dim, niac, pair_i, pair_j, dwdx, eniac, epair_a, epair_s, del_gamma_as, mass, rho, itype, ntotal, etotal)

        ! For periodic condition, update dfncn
        if (Allocated(pBC_edges)) then
            do d=1,SPH_dim
                call PeriodicParameter(dfncn(d,:))
            enddo
        endif

        ! Determine function's laplacian  𝜵∙kk 𝜵(f_a) using USAW formulation    
        ! we set kk as 1 in 𝜵∙kk 𝜵(f_a)
        kk=1.D0
        
        call LaplacianSPHTraditional(lap_fncn, fncn, kk, x(1:SPH_dim,1:ntotal), &
                 & gamma_discrt, SPH_dim, niac, pair_i, pair_j, dwdx, mass, rho, itype, ntotal)        
        
        ALLOCATE(x_s(SPH_dim,etotal))
        do s= 1, etotal
           x_s(:,s)=x(:,nedge_rel_edge(s))
           kk_s(s)=1.D0
        enddo

        call BILMEA(lap_fncn, fncn, fncn_s, x(1:SPH_dim,1:ntotal), x_s, gamma_discrt, SPH_dim, eniac, epair_a, epair_s, &
                &  del_gamma_as, etype, etotal, ntotal,  etype_periodic, etype_thermal_neumann,etype_thermal_dirichlet, &
                & bdryVal_temp,etotal, surf_norm, kk, kk_s)
      
        deallocate(x_s)
        
        
    endif

    
    
 !---------------------------------------------------------------------------------------------------------------------------


deallocate(dfncn, dfncn_s)
deallocate(fncn_s)

    
end       
