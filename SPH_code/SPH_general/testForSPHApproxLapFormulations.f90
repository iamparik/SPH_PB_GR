!
!  SUBROUTINE: testForSPHApproxFormulations
!
!  PURPOSE:  Subroutine that calls various SPH operators 
!           to find interpolations of different functions
!
!   CREATED:        05/02/2022       by  PARIKSHIT BOREGOWDA
!   Last Modified:  01/28/2023       by  PARIKSHIT BOREGOWDA 
!****************************************************************************
    
subroutine testForSPHApproxLapFormulations 
    use config_parameter, only: SPH_dim, dataOutputPath, &
            & etype_periodic, etype_SolidWall1, etype_SolidWall2 
    use config_geometry, only:  dx_r
    use particle_data ,   only: rho, mass, itype, x, nreal,nedge,nghost,ntotal,&
            & niac, pair_i, pair_j, w, w_aa, dwdx,     &
            & eniac, epair_a, epair_s,nedge_rel_edge, edge, &
            & gamma_cont, gamma_mat, gamma_discrt, del_gamma_as, &
            & gamma_mat_inv, xi1_mat_inv, xi_cont_mat_inv, etotal,etype, &
            & nperiodic,pBC_duplicate_pair, pBC_edges, &
            & hsml_factor_name, hsml_factor_code
    implicit none
    
    integer(4) a, k, d, s
    real(8) kk(ntotal), kk_s(etotal)
    real(8),DIMENSION(:),ALLOCATABLE:: fncn, fncn0, fncn_s0, lap_fncn, lap_fncn0
    real(8),DIMENSION(:,:),ALLOCATABLE:: dfncn, bdryVal
    real(8),DIMENSION(:,:),ALLOCATABLE:: dfncn0, dfncn_s, x_s
    real(8),DIMENSION(:,:,:),ALLOCATABLE:: gamma_mat_approx_inv
    real(8),DIMENSION(:),ALLOCATABLE:: gamma_approx
    real(8)::RE(2,30),REB(2,30),REI(2,30)
    LOGICAL :: relErrorLogical
    real(8) G_s(etotal), G(ntotal)
    
    ALLOCATE(fncn(ntotal),fncn0(ntotal), lap_fncn(ntotal), lap_fncn0(ntotal))
    ALLOCATE( dfncn(SPH_dim,ntotal), dfncn0(SPH_dim,ntotal))
    fncn=0.D0
    dfncn=0.D0
    lap_fncn=0.D0
    G=1.D0
    
    relErrorLogical = .TRUE.
    ! Initialize the test function to be whatever is to be tested
    ! This can be changed to calculate different functions
    

    do a =1,ntotal
       fncn0(a)=(x(1,a)**6+x(2,a)**6)
       do d=1,SPH_dim
            dfncn0(d,a)=6.D0*(x(d,a)**5)
       enddo
       lap_fncn0(a)=30.D0*(x(1,a)**4+x(2,a)**4)
    enddo
    
    !do a =1,ntotal
    !   fncn0(a)=(x(1,a)**2+x(2,a)**2) 
    !   do d=1,SPH_dim
    !        dfncn0(d,a)=2.D0*x(d,a)
    !   enddo
    !   lap_fncn0(a)=4.D0
    !enddo

    
    
    ALLOCATE(fncn_s0(etotal), dfncn_s(SPH_dim,etotal),x_s(SPH_dim,etotal))
    fncn_s0=0.D0
    dfncn_s=0.D0
    x_s=0.D0
    
    ! This evaluates function value at the boundary. This can be further generalized
    do s= 1, etotal
       fncn_s0(s)=fncn0(nedge_rel_edge(s))
       x_s(:,s)=x(:,nedge_rel_edge(s))
    enddo
    
    if (Allocated(pBC_edges)) then
        do a = 1, nperiodic
            fncn0(pBC_duplicate_pair(2,a))=fncn0(pBC_duplicate_pair(1,a))
            dfncn0(:,pBC_duplicate_pair(2,a))=dfncn0(:,pBC_duplicate_pair(1,a))
            lap_fncn0(pBC_duplicate_pair(2,a))=lap_fncn0(pBC_duplicate_pair(1,a))
        enddo
    endif
    
 !---------------------------------------------------- KCBI Formulation ------------------------------------------------   
    write(*,*) "********************************************************************************************************* "
    write(*,*) "KCBI approximations for given function, BIL CASE = 01....  "
    ! Determine function approximation using KCBI formulation
    fncn=0.D0
    call zerothOrderFunction(fncn0, fncn, gamma_cont, niac, pair_i, pair_j, w, w_aa, mass, rho, itype, ntotal)
    
    !ALLOCATE( gamma_mat_approx_inv(SPH_dim,SPH_dim,ntotal)) 
    ! Determine function's first order derivative approximation using KCBI formulation
    dfncn=0.D0
    call GradientKCBI( dfncn, fncn0, fncn_s0, gamma_cont, SPH_dim, niac, pair_i, pair_j, dwdx, eniac, epair_a, epair_s, del_gamma_as, mass, rho, itype, ntotal, etotal)

    ! Determine function's laplacian using KCBI as Divergence of gradient 𝜵∙𝜵(f_a) 
    lap_fncn=0.D0
    ! This evaluates function first order derivative value at the boundary. This can be further generalized for numerical integration along edge...
    do s= 1, etotal
       !dfncn_s(:,s)=dfncn(:,nedge_rel_edge(s))
       dfncn_s(:,s)=dfncn0(:,nedge_rel_edge(s))
    enddo
    call DivergenceKCBI( lap_fncn, dfncn, dfncn_s, gamma_cont, SPH_dim, niac, pair_i, pair_j, dwdx, eniac, epair_a, epair_s, del_gamma_as, mass, rho, itype, ntotal, etotal)

    !relative error
    call relErrorFormulation(lap_fncn0,lap_fncn,RE(:,1),REB(:,1),REI(:,1),relErrorLogical)
    
    ! Output of the KCBI formulation
    call outputFncnApprx(int(1),x, fncn0, fncn, dfncn0, dfncn, lap_fncn0, lap_fncn, nedge_rel_edge(1:etotal),edge(SPH_dim,1:etotal), etotal,SPH_dim,ntotal,itype(1:ntotal),etype(1:etotal),dataOutputPath)
    
    
  !---------------------------------------------------- PCBI Formulation ------------------------------------------------   
    write(*,*) "********************************************************************************************************* "
    write(*,*) "PCBI approximations for given function, BIL CASE = 02 ....  "
    ! Determine function approximation using PCBI formulation
    fncn=0.D0
    call zerothOrderFunction(fncn0, fncn, gamma_discrt, niac, pair_i, pair_j, w, w_aa, mass, rho, itype, ntotal)
    
    !ALLOCATE( gamma_mat_approx_inv(SPH_dim,SPH_dim,ntotal)) 
    ! Determine function's first order derivative approximation using PCBI formulation
    dfncn=0.D0
    call GradientPCBI( dfncn, fncn0, fncn_s0, gamma_mat_inv, SPH_dim, niac, pair_i, pair_j, dwdx, eniac, epair_a, epair_s, del_gamma_as, mass, rho, itype, ntotal, etotal)

    ! Determine function's laplacian using PCBI as Divergence of gradient 𝜵∙𝜵(f_a) 
    lap_fncn=0.D0
    ! This evaluates function first order derivative value at the boundary. This can be further generalized for numerical integration along edge...
    do s= 1, etotal
       !dfncn_s(:,s)=dfncn(:,nedge_rel_edge(s))
       dfncn_s(:,s)=dfncn0(:,nedge_rel_edge(s))
    enddo
    call DivergencePCBI( lap_fncn, dfncn, dfncn_s, gamma_mat_inv, SPH_dim, niac, pair_i, pair_j, dwdx, eniac, epair_a, epair_s, del_gamma_as, mass, rho, itype, ntotal, etotal)

    !relative error
    !call relErrorFormulation(lap_fncn0,lap_fncn,int(2))
    call relErrorFormulation(lap_fncn0,lap_fncn,RE(:,2),REB(:,2),REI(:,2),relErrorLogical)
    
    
    ! Output of the PCBI formulation
    call outputFncnApprx(int(2),x, fncn0, fncn, dfncn0, dfncn, lap_fncn0, lap_fncn, nedge_rel_edge(1:etotal),edge(SPH_dim,1:etotal), etotal,SPH_dim,ntotal,itype(1:ntotal),etype(1:etotal),dataOutputPath)
     
    
    !---------------------------------------------------- PCBI based weakly consistent Formulation I ------------------------------------------------   
    write(*,*) "********************************************************************************************************* "
    write(*,*) "WCBI1 approximations for given function, BIL CASE = 03 ....  "
    ! Determine function approximation using WCBI1 formulation
    fncn=0.D0
    call zerothOrderFunction(fncn0, fncn, gamma_discrt, niac, pair_i, pair_j, w, w_aa, mass, rho, itype, ntotal)
    
    ALLOCATE( gamma_mat_approx_inv(SPH_dim,SPH_dim,ntotal)) 
    gamma_mat_approx_inv=0.D0
    do a=1,ntotal
        do d=1,SPH_dim
                gamma_mat_approx_inv(d,d,a)= 1.D0/gamma_cont(a)
        enddo
    enddo
    
    !ALLOCATE( gamma_mat_approx_inv(SPH_dim,SPH_dim,ntotal)) 
    ! Determine function's first order derivative approximation using WCBI1 formulation
    dfncn=0.D0
    call GradientPCBI( dfncn, fncn0, fncn_s0, gamma_mat_approx_inv, SPH_dim, niac, pair_i, pair_j, dwdx, eniac, epair_a, epair_s, del_gamma_as, mass, rho, itype, ntotal, etotal)

    ! Determine function's laplacian using WCBI1 as Divergence of gradient 𝜵∙𝜵(f_a) 
    lap_fncn=0.D0
    ! This evaluates function first order derivative value at the boundary. This can be further generalized for numerical integration along edge...
    do s= 1, etotal
       !dfncn_s(:,s)=dfncn(:,nedge_rel_edge(s))
       dfncn_s(:,s)=dfncn0(:,nedge_rel_edge(s))
    enddo
    call DivergencePCBI( lap_fncn, dfncn, dfncn_s, gamma_mat_approx_inv, SPH_dim, niac, pair_i, pair_j, dwdx, eniac, epair_a, epair_s, del_gamma_as, mass, rho, itype, ntotal, etotal)

    !relative error
    call relErrorFormulation(lap_fncn0,lap_fncn,RE(:,3),REB(:,3),REI(:,3),relErrorLogical)  
        
     DEALLOCATE( gamma_mat_approx_inv )
    
    ! Output of the WCBI1 formulation
    call outputFncnApprx(int(3),x, fncn0, fncn, dfncn0, dfncn, lap_fncn0, lap_fncn, nedge_rel_edge(1:etotal),edge(SPH_dim,1:etotal), etotal,SPH_dim,ntotal,itype(1:ntotal),etype(1:etotal),dataOutputPath)
  
!---------------------------------------------------- PCBI based weakly consistent Formulation 2------------------------------------------------   
    write(*,*) "********************************************************************************************************* "
    write(*,*) "WCBI2 approximations for given function, BIL CASE = 04 ....  "
    ! Determine function approximation using WCBI2 formulation
    fncn=0.D0
    call zerothOrderFunction(fncn0, fncn, gamma_discrt, niac, pair_i, pair_j, w, w_aa, mass, rho, itype, ntotal)
    
    ALLOCATE( gamma_mat_approx_inv(SPH_dim,SPH_dim,ntotal)) 
    gamma_mat_approx_inv=0.D0
    do a=1,ntotal
        do d=1,SPH_dim
                gamma_mat_approx_inv(d,d,a)= 1.D0/gamma_discrt(a)
        enddo
    enddo
    
    !ALLOCATE( gamma_mat_approx_inv(SPH_dim,SPH_dim,ntotal)) 
    ! Determine function's first order derivative approximation using WCBI2 formulation
    dfncn=0.D0
    call GradientPCBI( dfncn, fncn0, fncn_s0, gamma_mat_approx_inv, SPH_dim, niac, pair_i, pair_j, dwdx, eniac, epair_a, epair_s, del_gamma_as, mass, rho, itype, ntotal, etotal)

    ! Determine function's laplacian using WCBI2 as Divergence of gradient 𝜵∙𝜵(f_a) 
    lap_fncn=0.D0
    ! This evaluates function first order derivative value at the boundary. This can be further generalized for numerical integration along edge...
    do s= 1, etotal
       !dfncn_s(:,s)=dfncn(:,nedge_rel_edge(s))
       dfncn_s(:,s)=dfncn0(:,nedge_rel_edge(s))
    enddo
    call DivergencePCBI( lap_fncn, dfncn, dfncn_s, gamma_mat_approx_inv, SPH_dim, niac, pair_i, pair_j, dwdx, eniac, epair_a, epair_s, del_gamma_as, mass, rho, itype, ntotal, etotal)

    !relative error
    call relErrorFormulation(lap_fncn0,lap_fncn,RE(:,4),REB(:,4),REI(:,4),relErrorLogical)
        
     DEALLOCATE( gamma_mat_approx_inv )
    
    ! Output of the WCBI2 formulation
    call outputFncnApprx(int(4),x, fncn0, fncn, dfncn0, dfncn, lap_fncn0, lap_fncn, nedge_rel_edge(1:etotal),edge(SPH_dim,1:etotal), etotal,SPH_dim,ntotal,itype(1:ntotal),etype(1:etotal),dataOutputPath)

!---------------------------------------------------- PCBI based weakly consistent Formulation 3------------------------------------------------   
    write(*,*) "********************************************************************************************************* "
    write(*,*) "WCBI3 approximations for given function, BIL CASE = 05 ....  "
    ! Determine function approximation using WCBI3 formulation
    fncn=0.D0
    call zerothOrderFunction(fncn0, fncn, gamma_discrt, niac, pair_i, pair_j, w, w_aa, mass, rho, itype, ntotal)
    
    ALLOCATE( gamma_mat_approx_inv(SPH_dim,SPH_dim,ntotal)) 
    gamma_mat_approx_inv=0.D0
    do a=1,ntotal
        do d=1,SPH_dim
                gamma_mat_approx_inv(d,d,a)= 1.D0/gamma_mat(d,d,a)
        enddo
    enddo
    
    !ALLOCATE( gamma_mat_approx_inv(SPH_dim,SPH_dim,ntotal)) 
    ! Determine function's first order derivative approximation using WCBI3 formulation
    dfncn=0.D0
    call GradientPCBI( dfncn, fncn0, fncn_s0, gamma_mat_approx_inv, SPH_dim, niac, pair_i, pair_j, dwdx, eniac, epair_a, epair_s, del_gamma_as, mass, rho, itype, ntotal, etotal)

    ! Determine function's laplacian using WCBI3 as Divergence of gradient 𝜵∙𝜵(f_a) 
    lap_fncn=0.D0
    ! This evaluates function first order derivative value at the boundary. This can be further generalized for numerical integration along edge...
    do s= 1, etotal
       !dfncn_s(:,s)=dfncn(:,nedge_rel_edge(s))
       dfncn_s(:,s)=dfncn0(:,nedge_rel_edge(s))
    enddo
    call DivergencePCBI( lap_fncn, dfncn, dfncn_s, gamma_mat_approx_inv, SPH_dim, niac, pair_i, pair_j, dwdx, eniac, epair_a, epair_s, del_gamma_as, mass, rho, itype, ntotal, etotal)

    !relative error
    call relErrorFormulation(lap_fncn0,lap_fncn,RE(:,5),REB(:,5),REI(:,5),relErrorLogical)   
        
     DEALLOCATE( gamma_mat_approx_inv )
    
    ! Output of the WCBI3 formulation
    call outputFncnApprx(int(5),x, fncn0, fncn, dfncn0, dfncn, lap_fncn0, lap_fncn, nedge_rel_edge(1:etotal),edge(SPH_dim,1:etotal), etotal,SPH_dim,ntotal,itype(1:ntotal),etype(1:etotal),dataOutputPath)
    
    
!---------------------------------------------------- CSPM Formulation ------------------------------------------------   
    write(*,*) "********************************************************************************************************* "
    write(*,*) "CSPM approximations for given function, BIL CASE = 06....  "
    ! Determine function approximation using CSPM formulation
    fncn=0.D0
    call zerothOrderFunction(fncn0, fncn, gamma_discrt, niac, pair_i, pair_j, w, w_aa, mass, rho, itype, ntotal)
    
    !ALLOCATE( gamma_mat_approx_inv(SPH_dim,SPH_dim,ntotal)) 
    ! Determine function's first order derivative approximation using CSPM formulation
    dfncn=0.D0
    call GradientCSPM( dfncn, fncn0, xi1_mat_inv, SPH_dim, niac, pair_i, pair_j, dwdx, mass, rho, itype, ntotal)

    ! Determine function's laplacian using CSPM as Divergence of gradient 𝜵∙𝜵(f_a) 
    lap_fncn=0.D0
  
    call DivergenceCSPM(lap_fncn, dfncn, xi1_mat_inv, SPH_dim, niac, pair_i, pair_j, dwdx, mass, rho, itype, ntotal)

    !relative error
    call relErrorFormulation(lap_fncn0,lap_fncn,RE(:,6),REB(:,6),REI(:,6),relErrorLogical)
    
    ! Output of the CSPM formulation
    call outputFncnApprx(int(6),x, fncn0, fncn, dfncn0, dfncn, lap_fncn0, lap_fncn, nedge_rel_edge(1:etotal),edge(SPH_dim,1:etotal), etotal,SPH_dim,ntotal,itype(1:ntotal),etype(1:etotal),dataOutputPath)
    
 !---------------------------------------------------- KCBI-CSPM Formulation ------------------------------------------------   
    write(*,*) "********************************************************************************************************* "
    write(*,*) "CSPM approximations for given function, BIL CASE = 07 ....  "
    ! Determine function approximation using KCBI-CSPM formulation
    fncn=0.D0
    call zerothOrderFunction(fncn0, fncn, gamma_discrt, niac, pair_i, pair_j, w, w_aa, mass, rho, itype, ntotal)
    
    !ALLOCATE( gamma_mat_approx_inv(SPH_dim,SPH_dim,ntotal)) 
    ! Determine function's first order derivative approximation using KCBI-CSPM formulation
    dfncn=0.D0
    call GradientKCBI_CSPM(dfncn, fncn0, xi1_mat_inv, SPH_dim, niac, pair_i, pair_j, dwdx, eniac, epair_a, del_gamma_as, mass, rho, itype, ntotal, etotal)


    ! Determine function's laplacian using KCBI-CSPM as Divergence of gradient 𝜵∙𝜵(f_a) 
    lap_fncn=0.D0
  
    call DivergenceKCBI_CSPM(lap_fncn, dfncn, xi1_mat_inv, SPH_dim, niac, pair_i, pair_j, dwdx, eniac, epair_a, del_gamma_as, mass, rho, itype, ntotal, etotal)

    !relative error
    call relErrorFormulation(lap_fncn0,lap_fncn,RE(:,7),REB(:,7),REI(:,7),relErrorLogical)
    
    ! Output of the KCBI-CSPM formulation
    call outputFncnApprx(int(7),x, fncn0, fncn, dfncn0, dfncn, lap_fncn0, lap_fncn, nedge_rel_edge(1:etotal),edge(SPH_dim,1:etotal), etotal,SPH_dim,ntotal,itype(1:ntotal),etype(1:etotal),dataOutputPath)

!---------------------------------------------------- BIL-USAW Formulation (using analytical gradients)------------------------------------------------   
    write(*,*) "********************************************************************************************************* "
    write(*,*) "BIL-USAW Formulation BIL CASE = 08 ....  "
    ! Values of functions are used as analytical values
    fncn=0.D0
    fncn=fncn0
    
    !Values of gradients are analytical values at boundary and for particles close to boundary
    dfncn=0.D0
    dfncn=dfncn0

    ! Determine function's laplacian  𝜵∙kk 𝜵(f_a) using USAW formulation
    lap_fncn=0.D0
    
    ! Boundary vales of fucntions are its analytical values at the boundary
    do s= 1, etotal
       dfncn_s(:,s)=dfncn(:,nedge_rel_edge(s))
    enddo
    
    ! Determine function's laplacian  𝜵∙kk 𝜵(f_a) using USAW formulation    
    ! we set kk as 1 in 𝜵∙kk 𝜵(f_a)
    kk=1.D0
        
    call LaplacianSPHTraditional(lap_fncn, fncn, kk, x(1:SPH_dim,1:ntotal), &
                 & gamma_cont, SPH_dim, niac, pair_i, pair_j, dwdx, mass, rho, itype, ntotal)
    
    allocate(bdryVal(SPH_dim,etotal))
    bdryVal=0.D0
    call BIL1(lap_fncn, dfncn, gamma_cont, SPH_dim, eniac, epair_a, epair_s, &
                        & del_gamma_as, etype, etotal, ntotal, etype_periodic, &
                        & etype_SolidWall2,etype_SolidWall1,bdryVal, 1.D0)
    call BIL2(lap_fncn, dfncn_s, gamma_cont, SPH_dim, eniac, epair_a, epair_s, &
                        & del_gamma_as, etype, etotal, ntotal, etype_periodic, &
                        & etype_SolidWall2,etype_SolidWall1,bdryVal,etotal, 1.D0)
    deallocate(bdryval)
    
    !relative error
    call relErrorFormulation(lap_fncn0,lap_fncn,RE(:,8),REB(:,8),REI(:,8),relErrorLogical)
    
    ! Output of the KCBI-CSPM formulation
    call outputFncnApprx(int(8),x, fncn0, fncn, dfncn0, dfncn, lap_fncn0, lap_fncn, nedge_rel_edge(1:etotal),edge(SPH_dim,1:etotal), etotal,SPH_dim,ntotal,itype(1:ntotal),etype(1:etotal),dataOutputPath)

!---------------------------------------------------- BIL1 Formulation (using analytical gradients)------------------------------------------------   
    write(*,*) "********************************************************************************************************* "
    write(*,*) "BIL1 Formulation BIL CASE = 09 ....  "
    
    ! Values of functions are used as analytical values
    fncn=0.D0
    fncn=fncn0
    
    !Values of gradients are analytical values at boundary and for particles close to boundary
    dfncn=0.D0
    dfncn=dfncn0

    ! Determine function's laplacian  𝜵∙ 𝜵(f_a) using BIL1 formulation where β is assumed to be negligible in taylor series to obtan lap
    lap_fncn=0.D0
  
    ! Determine function's laplacian  𝜵∙kk 𝜵(f_a) using BIL1 formulation    
        ! we set kk as 1 in 𝜵∙kk 𝜵(f_a)
    kk=1.D0
        
    call LaplacianSPHTraditional(lap_fncn, fncn, kk, x(1:SPH_dim,1:ntotal), &
                 & gamma_cont, SPH_dim, niac, pair_i, pair_j, dwdx, mass, rho, itype, ntotal)  

    
    call BIL1(lap_fncn, dfncn, gamma_cont, SPH_dim, eniac, epair_a, epair_s, &
                        & del_gamma_as, etype, etotal, ntotal, etype_periodic, &
                        & etype_SolidWall2,etype_SolidWall1,bdryVal, 2.D0)

    !relative error
    call relErrorFormulation(lap_fncn0,lap_fncn,RE(:,9),REB(:,9),REI(:,9),relErrorLogical)
    
    ! Output of the BIL1 formulation
    call outputFncnApprx(int(9),x, fncn0, fncn, dfncn0, dfncn, lap_fncn0, lap_fncn, nedge_rel_edge(1:etotal),edge(SPH_dim,1:etotal), etotal,SPH_dim,ntotal,itype(1:ntotal),etype(1:etotal),dataOutputPath)

!---------------------------------------------------- BIL2 Formulation (using analytical gradients)------------------------------------------------   
    write(*,*) "********************************************************************************************************* "
    write(*,*) "BIL2 Formulation (using analytical gradients) to find Laplacian BIL CASE =10....   "
    
    ! Values of functions are used as analytical values
    fncn=0.D0
    fncn=fncn0
    
    !Values of gradients are analytical values at boundary and for particles close to boundary
    dfncn=0.D0
    dfncn=dfncn0

    ! Determine function's laplacian  𝜵∙ 𝜵(f_a) using BIL2 formulation 
    lap_fncn=0.D0

    ! Boundary vales of fucntions are its analytical values at the boundary
    do s= 1, etotal
       dfncn_s(:,s)=dfncn(:,nedge_rel_edge(s))
    enddo
  
    kk=1.D0
        
    call LaplacianSPHTraditional(lap_fncn, fncn, kk, x(1:SPH_dim,1:ntotal), &
                 & gamma_cont, SPH_dim, niac, pair_i, pair_j, dwdx, mass, rho, itype, ntotal)
    
    allocate(bdryVal(SPH_dim,etotal))
    bdryVal=0.D0
    call BIL2(lap_fncn, dfncn_s, gamma_cont, SPH_dim, eniac, epair_a, epair_s, &
                        & del_gamma_as, etype, etotal, ntotal, etype_periodic, &
                        & etype_SolidWall2,etype_SolidWall1,bdryVal,etotal, 2.D0)
    deallocate(bdryval)
    !relative error
    call relErrorFormulation(lap_fncn0,lap_fncn,RE(:,10),REB(:,10),REI(:,10),relErrorLogical)
    
    ! Output of the BIL2 formulation
    call outputFncnApprx(int(10),x, fncn0, fncn, dfncn0, dfncn, lap_fncn0, lap_fncn, nedge_rel_edge(1:etotal),edge(SPH_dim,1:etotal), etotal,SPH_dim,ntotal,itype(1:ntotal),etype(1:etotal),dataOutputPath)

!---------------------------------------------------- BIL-Macia Formulation ------------------------------------------------   
    write(*,*) "********************************************************************************************************* "
    write(*,*) "BIL-Macia Formulation (using analytical gradients) to find Laplacian BIL CASE =11....   "
    
    ! Values of functions are used as analytical values
    fncn=0.D0
    fncn=fncn0
    
    !Values of gradients are analytical values at boundary and for particles close to boundary
    dfncn=0.D0
    dfncn=dfncn0

    ! Determine function's laplacian  𝜵∙ 𝜵(f_a) using BIL1 formulation where β is assumed to be negligible in taylor series to obtan lap
    lap_fncn=0.D0
  
    ! Determine function's laplacian  𝜵∙kk 𝜵(f_a) using BIL formulation    
        ! we set kk as 1 in 𝜵∙kk 𝜵(f_a)
    kk=1.D0
        
    call LaplacianSPHTraditional(lap_fncn, fncn, kk, x(1:SPH_dim,1:ntotal), &
                 & gamma_cont, SPH_dim, niac, pair_i, pair_j, dwdx, mass, rho, itype, ntotal)  

    ! Boundary vales of fucntions are its analytical values at the boundary
    do s= 1, etotal
       kk_s(s)=1.D0
    enddo
    
    call BILMacia(lap_fncn, fncn0, fncn_s0, x, x_s, gamma_cont, SPH_dim, eniac,&
            & epair_a, epair_s, del_gamma_as, etype, etotal, ntotal, kk_s)
    
    !relative error
    call relErrorFormulation(lap_fncn0,lap_fncn,RE(:,11),REB(:,11),REI(:,11),relErrorLogical)
    
    ! Output of the BIL2 formulation
    call outputFncnApprx(int(11),x, fncn0, fncn, dfncn0, dfncn, lap_fncn0, lap_fncn, nedge_rel_edge(1:etotal),edge(SPH_dim,1:etotal), etotal,SPH_dim,ntotal,itype(1:ntotal),etype(1:etotal),dataOutputPath)

!---------------------------------------------------- BIL-Macia Formulation with 2X boundary ------------------------------------------------   
    write(*,*) "********************************************************************************************************* "
    write(*,*) "BIL-Macia Formulation BIL CASE =12....   "
    
    ! Values of functions are used as analytical values
    fncn=0.D0
    fncn=fncn0
    
    !Values of gradients are analytical values at boundary and for particles close to boundary
    dfncn=0.D0
    dfncn=dfncn0

    ! Determine function's laplacian  𝜵∙ 𝜵(f_a) using BIL1 formulation where β is assumed to be negligible in taylor series to obtan lap
    lap_fncn=0.D0
  
    ! Determine function's laplacian  𝜵∙kk 𝜵(f_a) using BIL1 formulation    
        ! we set kk as 1 in 𝜵∙kk 𝜵(f_a)
    kk=1.D0
        
    call LaplacianSPHTraditional(lap_fncn, fncn, kk, x(1:SPH_dim,1:ntotal), &
                 & gamma_cont, SPH_dim, niac, pair_i, pair_j, dwdx, mass, rho, itype, ntotal)  

    ! Boundary vales of fucntions are its analytical values at the boundary
    do s= 1, etotal
       kk_s(s)=2.D0
    enddo
    
    call BILMacia(lap_fncn, fncn0, fncn_s0, x, x_s, gamma_cont, SPH_dim, eniac,&
            & epair_a, epair_s, del_gamma_as, etype, etotal, ntotal, kk_s)
    
    !relative error
    call relErrorFormulation(lap_fncn0,lap_fncn,RE(:,12),REB(:,12),REI(:,12),relErrorLogical)
    
    ! Output of the BIL2 formulation
    call outputFncnApprx(int(12),x, fncn0, fncn, dfncn0, dfncn, lap_fncn0, lap_fncn, nedge_rel_edge(1:etotal),edge(SPH_dim,1:etotal), etotal,SPH_dim,ntotal,itype(1:ntotal),etype(1:etotal),dataOutputPath)

!---------------------------------------------------- BIL1 Formulation discrete gamma------------------------------------------------   
    write(*,*) "********************************************************************************************************* "
    write(*,*) "BIL1 Formulation gamma dsicrete BIL CASE =13....  "
    
    ! Values of functions are used as analytical values
    fncn=0.D0
    fncn=fncn0
    
    !Values of gradients are analytical values at boundary and for particles close to boundary
    dfncn=0.D0
    dfncn=dfncn0

    ! Determine function's laplacian  𝜵∙ 𝜵(f_a) using BIL1 formulation where β is assumed to be negligible in taylor series to obtan lap
    lap_fncn=0.D0
  
    ! Determine function's laplacian  𝜵∙kk 𝜵(f_a) using BIL1 formulation    
        ! we set kk as 1 in 𝜵∙kk 𝜵(f_a)
    kk=1.D0
        
    call LaplacianSPHTraditional(lap_fncn, fncn, kk, x(1:SPH_dim,1:ntotal), &
                 & gamma_discrt, SPH_dim, niac, pair_i, pair_j, dwdx, mass, rho, itype, ntotal)  

    
    call BIL1(lap_fncn, dfncn, gamma_discrt, SPH_dim, eniac, epair_a, epair_s, &
                        & del_gamma_as, etype, etotal, ntotal, etype_periodic, &
                        & etype_SolidWall2,etype_SolidWall1,bdryVal, 2.D0)

    !relative error
    call relErrorFormulation(lap_fncn0,lap_fncn,RE(:,13),REB(:,13),REI(:,13),relErrorLogical)
    
    ! Output of the BIL1 formulation
    call outputFncnApprx(int(13),x, fncn0, fncn, dfncn0, dfncn, lap_fncn0, lap_fncn, nedge_rel_edge(1:etotal),edge(SPH_dim,1:etotal), etotal,SPH_dim,ntotal,itype(1:ntotal),etype(1:etotal),dataOutputPath)

!---------------------------------------------------- BIL2 Formulation  discrete gamma------------------------------------------------   
    write(*,*) "********************************************************************************************************* "
    write(*,*) "BIL2 Formulation gamma dsicrete BIL CASE =14....   "
    
    ! Values of functions are used as analytical values
    fncn=0.D0
    fncn=fncn0
    
    !Values of gradients are analytical values at boundary and for particles close to boundary
    dfncn=0.D0
    dfncn=dfncn0

    ! Determine function's laplacian  𝜵∙ 𝜵(f_a) using BIL2 formulation 
    lap_fncn=0.D0

    ! Boundary vales of fucntions are its analytical values at the boundary
    do s= 1, etotal
       dfncn_s(:,s)=dfncn(:,nedge_rel_edge(s))
    enddo
  
    kk=1.D0
        
    call LaplacianSPHTraditional(lap_fncn, fncn, kk, x(1:SPH_dim,1:ntotal), &
                 & gamma_discrt, SPH_dim, niac, pair_i, pair_j, dwdx, mass, rho, itype, ntotal)
    
    allocate(bdryVal(SPH_dim,etotal))
    bdryVal=0.D0
    call BIL2(lap_fncn, dfncn_s, gamma_discrt, SPH_dim, eniac, epair_a, epair_s, &
                        & del_gamma_as, etype, etotal, ntotal, etype_periodic, &
                        & etype_SolidWall2,etype_SolidWall1,bdryVal,etotal, 2.D0)
    deallocate(bdryval)
    !relative error
    call relErrorFormulation(lap_fncn0,lap_fncn,RE(:,14),REB(:,14),REI(:,14),relErrorLogical)
    
    ! Output of the BIL2 formulation
    call outputFncnApprx(int(14),x, fncn0, fncn, dfncn0, dfncn, lap_fncn0, lap_fncn, nedge_rel_edge(1:etotal),edge(SPH_dim,1:etotal), etotal,SPH_dim,ntotal,itype(1:ntotal),etype(1:etotal),dataOutputPath)



!---------------------------------------------------- BIL1 Corrected Formulation ------------------------------------------------   
    write(*,*) "********************************************************************************************************* "
    write(*,*) "BIL1 Corrected Formulation BIL CASE =15....  "
    
    ! Values of functions are used as analytical values
    fncn=0.D0
    fncn=fncn0
    
    !Values of gradients are analytical values at boundary and for particles close to boundary
    dfncn=0.D0
    dfncn=dfncn0

    ! Determine function's laplacian  𝜵∙ 𝜵(f_a) using BIL1 formulation where β is assumed to be negligible in taylor series to obtan lap
    lap_fncn=0.D0
    
    ! Determine function's laplacian  𝜵∙ 𝜵(f_a) using BIL1 formulation
    call LaplacianSPHTraditionalCorrected(lap_fncn, fncn, G, x(1:SPH_dim,1:ntotal), &
        & gamma_mat_inv, SPH_dim, niac, pair_i, pair_j, dwdx, mass, rho, itype, ntotal)
    call BIL1Corrected(lap_fncn, dfncn, gamma_mat_inv, SPH_dim, eniac, &
        & epair_a, epair_s, del_gamma_as, etype, etotal, ntotal, 2.D0)


    !relative error
    call relErrorFormulation(lap_fncn0,lap_fncn,RE(:,15),REB(:,15),REI(:,15),relErrorLogical)
    
    ! Output of the BIL1 formulation
    call outputFncnApprx(int(15),x, fncn0, fncn, dfncn0, dfncn, lap_fncn0, lap_fncn, nedge_rel_edge(1:etotal),edge(SPH_dim,1:etotal), etotal,SPH_dim,ntotal,itype(1:ntotal),etype(1:etotal),dataOutputPath)

!---------------------------------------------------- BIL2 Corrected Formulation ------------------------------------------------   
    write(*,*) "********************************************************************************************************* "
    write(*,*) "BIL2 Corrected Formulation BIL CASE =16....   "
    
    ! Values of functions are used as analytical values
    fncn=0.D0
    fncn=fncn0
    
    !Values of gradients are analytical values at boundary and for particles close to boundary
    dfncn=0.D0
    dfncn=dfncn0

    ! Determine function's laplacian  𝜵∙ 𝜵(f_a) using BIL2 formulation 
    lap_fncn=0.D0

    ! Determine function's laplacian  𝜵∙ 𝜵(f_a) using BIL2 formulation
    call LaplacianSPHTraditionalCorrected(lap_fncn, fncn, G, x(1:SPH_dim,1:ntotal), &
        & gamma_mat_inv, SPH_dim, niac, pair_i, pair_j, dwdx, mass, rho, itype, ntotal)
    ! Boundary values of fucntions are its analytical values at the boundary
    do s= 1, etotal
        dfncn_s(:,s)=dfncn(:,nedge_rel_edge(s))
    enddo
        
    call BIL2Corrected(lap_fncn, dfncn_s, gamma_mat_inv,&
        & SPH_dim, eniac, epair_a, epair_s, del_gamma_as, etype, etotal, ntotal, 2.D0)   

    !relative error
    call relErrorFormulation(lap_fncn0,lap_fncn,RE(:,16),REB(:,16),REI(:,16),relErrorLogical)
    
    ! Output of the BIL2 formulation
    call outputFncnApprx(int(16),x, fncn0, fncn, dfncn0, dfncn, lap_fncn0, lap_fncn, nedge_rel_edge(1:etotal),edge(SPH_dim,1:etotal), etotal,SPH_dim,ntotal,itype(1:ntotal),etype(1:etotal),dataOutputPath)


!---------------------------------------------------- BIL1 Corrected Formulation inside only------------------------------------------------   
    write(*,*) "********************************************************************************************************* "
    write(*,*) "BIL1 Corrected Formulation BIL CASE =17....  "
    
    ! Values of functions are used as analytical values
    fncn=0.D0
    fncn=fncn0
    
    !Values of gradients are analytical values at boundary and for particles close to boundary
    dfncn=0.D0
    dfncn=dfncn0

    ! Determine function's laplacian  𝜵∙ 𝜵(f_a) using BIL1 formulation where β is assumed to be negligible in taylor series to obtan lap
    lap_fncn=0.D0
  
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
    call LaplacianSPHTraditionalCorrected(lap_fncn, fncn, G, x(1:SPH_dim,1:ntotal), &
        & gamma_mat_approx_inv, SPH_dim, niac, pair_i, pair_j, dwdx, mass, rho, itype, ntotal)
    call BIL1Corrected(lap_fncn, dfncn, gamma_mat_approx_inv, SPH_dim, eniac, &
        & epair_a, epair_s, del_gamma_as, etype, etotal, ntotal, 2.D0)

    deallocate(gamma_mat_approx_inv)    
    

    !relative error
    call relErrorFormulation(lap_fncn0,lap_fncn,RE(:,17),REB(:,17),REI(:,17),relErrorLogical)
    
    ! Output of the BIL1 formulation
    call outputFncnApprx(int(17),x, fncn0, fncn, dfncn0, dfncn, lap_fncn0, lap_fncn, nedge_rel_edge(1:etotal),edge(SPH_dim,1:etotal), etotal,SPH_dim,ntotal,itype(1:ntotal),etype(1:etotal),dataOutputPath)

!---------------------------------------------------- BIL2 Corrected Formulation inside only ------------------------------------------------   
    write(*,*) "********************************************************************************************************* "
    write(*,*) "BIL2 Corrected Formulation BIL CASE =18....   "
    
    ! Values of functions are used as analytical values
    fncn=0.D0
    fncn=fncn0
    
    !Values of gradients are analytical values at boundary and for particles close to boundary
    dfncn=0.D0
    dfncn=dfncn0

    ! Determine function's laplacian  𝜵∙ 𝜵(f_a) using BIL2 formulation 
    lap_fncn=0.D0

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
        
    ! Determine function's laplacian  𝜵∙ 𝜵(f_a) using BIL2 formulation
    call LaplacianSPHTraditionalCorrected(lap_fncn, fncn, G, x(1:SPH_dim,1:ntotal), &
        & gamma_mat_approx_inv, SPH_dim, niac, pair_i, pair_j, dwdx, mass, rho, itype, ntotal)
    ! Boundary values of fucntions are its analytical values at the boundary
    do s= 1, etotal
        dfncn_s(:,s)=dfncn(:,nedge_rel_edge(s))
    enddo
        
    call BIL2Corrected(lap_fncn, dfncn_s, gamma_mat_approx_inv,&
        & SPH_dim, eniac, epair_a, epair_s, del_gamma_as, etype, etotal, ntotal, 2.D0)   
    
    deallocate(gamma_mat_approx_inv)
        
    !relative error
    call relErrorFormulation(lap_fncn0,lap_fncn,RE(:,18),REB(:,18),REI(:,18),relErrorLogical)
    
    ! Output of the BIL2 formulation
    call outputFncnApprx(int(18),x, fncn0, fncn, dfncn0, dfncn, lap_fncn0, lap_fncn, nedge_rel_edge(1:etotal),edge(SPH_dim,1:etotal), etotal,SPH_dim,ntotal,itype(1:ntotal),etype(1:etotal),dataOutputPath)
 
!---------------------------------------------------- BIL1 Corrected Formulation, but gamma_cont for boundary integral------------------------------------------------   
    write(*,*) "********************************************************************************************************* "
    write(*,*) "BIL1 Corrected Formulation BIL CASE =19....  "
    
    ! Values of functions are used as analytical values
    fncn=0.D0
    fncn=fncn0
    
    !Values of gradients are analytical values at boundary and for particles close to boundary
    dfncn=0.D0
    dfncn=dfncn0

    ! Determine function's laplacian  𝜵∙ 𝜵(f_a) using BIL1 formulation where β is assumed to be negligible in taylor series to obtan lap
    lap_fncn=0.D0
  
    ! Determine function's laplacian  𝜵∙ 𝜵(f_a) using BIL1 formulation
    call LaplacianSPHTraditionalCorrected(lap_fncn, fncn, G, x(1:SPH_dim,1:ntotal), &
        & gamma_mat_inv, SPH_dim, niac, pair_i, pair_j, dwdx, mass, rho, itype, ntotal)
    
    call BIL1(lap_fncn, dfncn, gamma_cont, SPH_dim, eniac, epair_a, epair_s, &
                        & del_gamma_as, etype, etotal, ntotal, etype_periodic, &
                        & etype_SolidWall2,etype_SolidWall1,bdryVal, 2.D0)
    !relative error
    call relErrorFormulation(lap_fncn0,lap_fncn,RE(:,19),REB(:,19),REI(:,19),relErrorLogical)
    
    ! Output of the BIL1 formulation
    call outputFncnApprx(int(19),x, fncn0, fncn, dfncn0, dfncn, lap_fncn0, lap_fncn, nedge_rel_edge(1:etotal),edge(SPH_dim,1:etotal), etotal,SPH_dim,ntotal,itype(1:ntotal),etype(1:etotal),dataOutputPath)

!---------------------------------------------------- BIL2 Corrected, but gamma_cont for boundary integral------------------------------------------------   
    write(*,*) "********************************************************************************************************* "
    write(*,*) "BIL2 Corrected Formulation BIL CASE =20....   "
    
    ! Values of functions are used as analytical values
    fncn=0.D0
    fncn=fncn0
    
    !Values of gradients are analytical values at boundary and for particles close to boundary
    dfncn=0.D0
    dfncn=dfncn0

    ! Determine function's laplacian  𝜵∙ 𝜵(f_a) using BIL2 formulation 
    lap_fncn=0.D0

    ! Determine function's laplacian  𝜵∙ 𝜵(f_a) using BIL2 formulation
    call LaplacianSPHTraditionalCorrected(lap_fncn, fncn, G, x(1:SPH_dim,1:ntotal), &
        & gamma_mat_inv, SPH_dim, niac, pair_i, pair_j, dwdx, mass, rho, itype, ntotal)
    ! Boundary values of fucntions are its analytical values at the boundary
    do s= 1, etotal
        dfncn_s(:,s)=dfncn(:,nedge_rel_edge(s))
    enddo
        
    allocate(bdryVal(SPH_dim,etotal))
    bdryVal=0.D0
    call BIL2(lap_fncn, dfncn_s, gamma_cont, SPH_dim, eniac, epair_a, epair_s, &
                        & del_gamma_as, etype, etotal, ntotal, etype_periodic, &
                        & etype_SolidWall2,etype_SolidWall1,bdryVal,etotal, 2.D0)
    deallocate(bdryval)
    
    !relative error
    call relErrorFormulation(lap_fncn0,lap_fncn,RE(:,20),REB(:,20),REI(:,20),relErrorLogical)
    
    ! Output of the BIL2 formulation
    call outputFncnApprx(int(20),x, fncn0, fncn, dfncn0, dfncn, lap_fncn0, lap_fncn, nedge_rel_edge(1:etotal),edge(SPH_dim,1:etotal), etotal,SPH_dim,ntotal,itype(1:ntotal),etype(1:etotal),dataOutputPath)
  
    !---------------------------------------------------- BIL1 Corrected Formulation with Xi------------------------------------------------   
    write(*,*) "********************************************************************************************************* "
    write(*,*) "BIL1 Corrected Formulation BIL CASE =21....  "
    
    ! Values of functions are used as analytical values
    fncn=0.D0
    fncn=fncn0
    
    !Values of gradients are analytical values at boundary and for particles close to boundary
    dfncn=0.D0
    dfncn=dfncn0

    ! Determine function's laplacian  𝜵∙ 𝜵(f_a) using BIL1 formulation where β is assumed to be negligible in taylor series to obtan lap
    lap_fncn=0.D0
  
    ! Determine function's laplacian  𝜵∙ 𝜵(f_a) using BIL1 formulation
    call LaplacianSPHTraditionalCorrected(lap_fncn, fncn, G, x(1:SPH_dim,1:ntotal), &
        & xi1_mat_inv, SPH_dim, niac, pair_i, pair_j, dwdx, mass, rho, itype, ntotal)
    
    call BIL1(lap_fncn, dfncn, gamma_cont, SPH_dim, eniac, epair_a, epair_s, &
                        & del_gamma_as, etype, etotal, ntotal, etype_periodic, &
                        & etype_SolidWall2,etype_SolidWall1,bdryVal, 2.D0)
    !relative error
    call relErrorFormulation(lap_fncn0,lap_fncn,RE(:,21),REB(:,21),REI(:,21),relErrorLogical)
    
    ! Output of the BIL1 formulation
    call outputFncnApprx(int(21),x, fncn0, fncn, dfncn0, dfncn, lap_fncn0, lap_fncn, nedge_rel_edge(1:etotal),edge(SPH_dim,1:etotal), etotal,SPH_dim,ntotal,itype(1:ntotal),etype(1:etotal),dataOutputPath)

!---------------------------------------------------- BIL2 Corrected Formulation with Xi ------------------------------------------------   
    write(*,*) "********************************************************************************************************* "
    write(*,*) "BIL2 Corrected Formulation BIL CASE =22....   "
    
    ! Values of functions are used as analytical values
    fncn=0.D0
    fncn=fncn0
    
    !Values of gradients are analytical values at boundary and for particles close to boundary
    dfncn=0.D0
    dfncn=dfncn0

    ! Determine function's laplacian  𝜵∙ 𝜵(f_a) using BIL2 formulation 
    lap_fncn=0.D0

    ! Determine function's laplacian  𝜵∙ 𝜵(f_a) using BIL2 formulation
    call LaplacianSPHTraditionalCorrected(lap_fncn, fncn, G, x(1:SPH_dim,1:ntotal), &
        & xi1_mat_inv, SPH_dim, niac, pair_i, pair_j, dwdx, mass, rho, itype, ntotal)
    ! Boundary values of fucntions are its analytical values at the boundary
    do s= 1, etotal
        dfncn_s(:,s)=dfncn(:,nedge_rel_edge(s))
    enddo
        
    allocate(bdryVal(SPH_dim,etotal))
    bdryVal=0.D0
    call BIL2(lap_fncn, dfncn_s, gamma_cont, SPH_dim, eniac, epair_a, epair_s, &
                        & del_gamma_as, etype, etotal, ntotal, etype_periodic, &
                        & etype_SolidWall2,etype_SolidWall1,bdryVal,etotal, 2.D0)
    deallocate(bdryval)
    
    !relative error
    call relErrorFormulation(lap_fncn0,lap_fncn,RE(:,22),REB(:,22),REI(:,22),relErrorLogical)
    
    ! Output of the BIL2 formulation
    call outputFncnApprx(int(22),x, fncn0, fncn, dfncn0, dfncn, lap_fncn0, lap_fncn, nedge_rel_edge(1:etotal),edge(SPH_dim,1:etotal), etotal,SPH_dim,ntotal,itype(1:ntotal),etype(1:etotal),dataOutputPath)
 
    
!---------------------------------------------------- BIL1 Corrected Formulation Xi With γ_a------------------------------------------------   
    write(*,*) "********************************************************************************************************* "
    write(*,*) "BIL1 Corrected Formulation BIL CASE =23....  "
    
    ! Values of functions are used as analytical values
    fncn=0.D0
    fncn=fncn0
    
    !Values of gradients are analytical values at boundary and for particles close to boundary
    dfncn=0.D0
    dfncn=dfncn0

    ! Determine function's laplacian  𝜵∙ 𝜵(f_a) using BIL1 formulation where β is assumed to be negligible in taylor series to obtan lap
    lap_fncn=0.D0
    
    ALLOCATE( gamma_mat_approx_inv(SPH_dim,SPH_dim,ntotal)) 
    gamma_mat_approx_inv=0.D0
    do a=1,ntotal
        gamma_mat_approx_inv(:,:,a)=xi1_mat_inv(:,:,a)/gamma_cont(a)
    enddo
  
    ! Determine function's laplacian  𝜵∙ 𝜵(f_a) using BIL1 formulation
    call LaplacianSPHTraditionalCorrected(lap_fncn, fncn, G, x(1:SPH_dim,1:ntotal), &
        & gamma_mat_approx_inv, SPH_dim, niac, pair_i, pair_j, dwdx, mass, rho, itype, ntotal)
    
    call BIL1(lap_fncn, dfncn, gamma_cont, SPH_dim, eniac, epair_a, epair_s, &
                        & del_gamma_as, etype, etotal, ntotal, etype_periodic, &
                        & etype_SolidWall2,etype_SolidWall1,bdryVal, 2.D0)
    
    deallocate(gamma_mat_approx_inv)
    
    !relative error
    call relErrorFormulation(lap_fncn0,lap_fncn,RE(:,23),REB(:,23),REI(:,23),relErrorLogical)
    
    ! Output of the BIL1 formulation
    call outputFncnApprx(int(23),x, fncn0, fncn, dfncn0, dfncn, lap_fncn0, lap_fncn, nedge_rel_edge(1:etotal),edge(SPH_dim,1:etotal), etotal,SPH_dim,ntotal,itype(1:ntotal),etype(1:etotal),dataOutputPath)

!---------------------------------------------------- BIL2 Corrected Formulation Xi With γ_a ------------------------------------------------   
    write(*,*) "********************************************************************************************************* "
    write(*,*) "BIL2 Corrected Formulation BIL CASE =24....   "
    
    ! Values of functions are used as analytical values
    fncn=0.D0
    fncn=fncn0
    
    !Values of gradients are analytical values at boundary and for particles close to boundary
    dfncn=0.D0
    dfncn=dfncn0

    ! Determine function's laplacian  𝜵∙ 𝜵(f_a) using BIL2 formulation 
    lap_fncn=0.D0
    ALLOCATE( gamma_mat_approx_inv(SPH_dim,SPH_dim,ntotal)) 
    gamma_mat_approx_inv=0.D0
    do a=1,ntotal
        gamma_mat_approx_inv(:,:,a)=xi1_mat_inv(:,:,a)/gamma_cont(a)
    enddo

    ! Determine function's laplacian  𝜵∙ 𝜵(f_a) using BIL2 formulation
    call LaplacianSPHTraditionalCorrected(lap_fncn, fncn, G, x(1:SPH_dim,1:ntotal), &
        & gamma_mat_approx_inv, SPH_dim, niac, pair_i, pair_j, dwdx, mass, rho, itype, ntotal)
    ! Boundary values of fucntions are its analytical values at the boundary
    do s= 1, etotal
        dfncn_s(:,s)=dfncn(:,nedge_rel_edge(s))
    enddo
        
    allocate(bdryVal(SPH_dim,etotal))
    bdryVal=0.D0
    call BIL2(lap_fncn, dfncn_s, gamma_cont, SPH_dim, eniac, epair_a, epair_s, &
                        & del_gamma_as, etype, etotal, ntotal, etype_periodic, &
                        & etype_SolidWall2,etype_SolidWall1,bdryVal,etotal, 2.D0)
    deallocate(bdryval,gamma_mat_approx_inv)
    
    !relative error
    call relErrorFormulation(lap_fncn0,lap_fncn,RE(:,24),REB(:,24),REI(:,24),relErrorLogical)
    
    ! Output of the BIL2 formulation
    call outputFncnApprx(int(24),x, fncn0, fncn, dfncn0, dfncn, lap_fncn0, lap_fncn, nedge_rel_edge(1:etotal),edge(SPH_dim,1:etotal), etotal,SPH_dim,ntotal,itype(1:ntotal),etype(1:etotal),dataOutputPath)
 
 !---------------------------------------------------- BIL-USAW Formulation corrected------------------------------------------------   
    write(*,*) "********************************************************************************************************* "
    write(*,*) "BIL-USAW Formulation BIL CASE = 25 ....  "
    ! Values of functions are used as analytical values
    fncn=0.D0
    fncn=fncn0
    
    !Values of gradients are analytical values at boundary and for particles close to boundary
    dfncn=0.D0
    dfncn=dfncn0

    ! Determine function's laplacian  𝜵∙kk 𝜵(f_a) using USAW formulation
    lap_fncn=0.D0
    
    ! Boundary vales of fucntions are its analytical values at the boundary
    do s= 1, etotal
       dfncn_s(:,s)=dfncn(:,nedge_rel_edge(s))
    enddo
    
    ! Determine function's laplacian  𝜵∙kk 𝜵(f_a) using USAW formulation    
    ! we set kk as 1 in 𝜵∙kk 𝜵(f_a)
    kk=1.D0
    
    call LaplacianSPHTraditionalCorrected(lap_fncn, fncn, kk, x(1:SPH_dim,1:ntotal),&
            & gamma_mat_inv, SPH_dim, niac, pair_i, pair_j, dwdx, mass, rho, itype, ntotal)
    
    call BIL1Corrected(lap_fncn, dfncn, gamma_mat_inv, SPH_dim, eniac,&
            & epair_a, epair_s, del_gamma_as, etype, etotal, ntotal, 1.D0)
    call BIL2Corrected(lap_fncn, dfncn_s, gamma_mat_inv, SPH_dim, eniac,&
            & epair_a, epair_s, del_gamma_as, etype, etotal, ntotal, 1.D0)
    
    !relative error
    call relErrorFormulation(lap_fncn0,lap_fncn,RE(:,25),REB(:,25),REI(:,25),relErrorLogical)
    
    ! Output of the KCBI-CSPM formulation
    call outputFncnApprx(int(25),x, fncn0, fncn, dfncn0, dfncn, lap_fncn0, lap_fncn, nedge_rel_edge(1:etotal),edge(SPH_dim,1:etotal), etotal,SPH_dim,ntotal,itype(1:ntotal),etype(1:etotal),dataOutputPath)
   
 !---------------------------------------------------- BIL-USAW Formulation corrected,  inside only------------------------------------------------   
    write(*,*) "********************************************************************************************************* "
    write(*,*) "BIL-USAW Formulation BIL CASE = 26 ....  "
    ! Values of functions are used as analytical values
    fncn=0.D0
    fncn=fncn0
    
    !Values of gradients are analytical values at boundary and for particles close to boundary
    dfncn=0.D0
    dfncn=dfncn0

    ! Determine function's laplacian  𝜵∙kk 𝜵(f_a) using USAW formulation
    lap_fncn=0.D0
    
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
    
    ! Boundary vales of fucntions are its analytical values at the boundary
    do s= 1, etotal
       dfncn_s(:,s)=dfncn(:,nedge_rel_edge(s))
    enddo
    
    ! Determine function's laplacian  𝜵∙kk 𝜵(f_a) using USAW formulation    
    ! we set kk as 1 in 𝜵∙kk 𝜵(f_a)
    kk=1.D0
    
    call LaplacianSPHTraditionalCorrected(lap_fncn, fncn, kk, x(1:SPH_dim,1:ntotal),&
            & gamma_mat_approx_inv, SPH_dim, niac, pair_i, pair_j, dwdx, mass, rho, itype, ntotal)
    
    allocate(bdryVal(SPH_dim,etotal))
    bdryVal=0.D0
    call BIL1(lap_fncn, dfncn, gamma_cont, SPH_dim, eniac, epair_a, epair_s, &
                        & del_gamma_as, etype, etotal, ntotal, etype_periodic, &
                        & etype_SolidWall2,etype_SolidWall1,bdryVal, 1.D0)
    call BIL2(lap_fncn, dfncn_s, gamma_cont, SPH_dim, eniac, epair_a, epair_s, &
                        & del_gamma_as, etype, etotal, ntotal, etype_periodic, &
                        & etype_SolidWall2,etype_SolidWall1,bdryVal,etotal, 1.D0)
    
    deallocate(bdryval,gamma_mat_approx_inv)
    
    !relative error
    call relErrorFormulation(lap_fncn0,lap_fncn,RE(:,26),REB(:,26),REI(:,26),relErrorLogical)
    
    ! Output of the KCBI-CSPM formulation
    call outputFncnApprx(int(26),x, fncn0, fncn, dfncn0, dfncn, lap_fncn0, lap_fncn, nedge_rel_edge(1:etotal),edge(SPH_dim,1:etotal), etotal,SPH_dim,ntotal,itype(1:ntotal),etype(1:etotal),dataOutputPath)

!---------------------------------------------------- BIL-USAW Formulation corrected, but gamma_cont for boundary integral------------------------------------------------   
    write(*,*) "********************************************************************************************************* "
    write(*,*) "BIL-USAW Formulation BIL CASE = 27 ....  "
    ! Values of functions are used as analytical values
    fncn=0.D0
    fncn=fncn0
    
    !Values of gradients are analytical values at boundary and for particles close to boundary
    dfncn=0.D0
    dfncn=dfncn0

    ! Determine function's laplacian  𝜵∙kk 𝜵(f_a) using USAW formulation
    lap_fncn=0.D0
    
    ! Boundary vales of fucntions are its analytical values at the boundary
    do s= 1, etotal
       dfncn_s(:,s)=dfncn(:,nedge_rel_edge(s))
    enddo
    
    ! Determine function's laplacian  𝜵∙kk 𝜵(f_a) using USAW formulation    
    ! we set kk as 1 in 𝜵∙kk 𝜵(f_a)
    kk=1.D0
    
    call LaplacianSPHTraditionalCorrected(lap_fncn, fncn, kk, x(1:SPH_dim,1:ntotal),&
            & gamma_mat_inv, SPH_dim, niac, pair_i, pair_j, dwdx, mass, rho, itype, ntotal)
    
    allocate(bdryVal(SPH_dim,etotal))
    bdryVal=0.D0
    call BIL1(lap_fncn, dfncn, gamma_cont, SPH_dim, eniac, epair_a, epair_s, &
                        & del_gamma_as, etype, etotal, ntotal, etype_periodic, &
                        & etype_SolidWall2,etype_SolidWall1,bdryVal, 1.D0)
    call BIL2(lap_fncn, dfncn_s, gamma_cont, SPH_dim, eniac, epair_a, epair_s, &
                        & del_gamma_as, etype, etotal, ntotal, etype_periodic, &
                        & etype_SolidWall2,etype_SolidWall1,bdryVal,etotal, 1.D0)
    deallocate(bdryval)
    
    !relative error
    call relErrorFormulation(lap_fncn0,lap_fncn,RE(:,27),REB(:,27),REI(:,27),relErrorLogical)
    
    ! Output of the KCBI-CSPM formulation
    call outputFncnApprx(int(27),x, fncn0, fncn, dfncn0, dfncn, lap_fncn0, lap_fncn, nedge_rel_edge(1:etotal),edge(SPH_dim,1:etotal), etotal,SPH_dim,ntotal,itype(1:ntotal),etype(1:etotal),dataOutputPath)

  !---------------------------------------------------- BIL-USAW Formulation corrected with gamma_discrete------------------------------------------------   
    write(*,*) "********************************************************************************************************* "
    write(*,*) "BIL-USAW Formulation BIL CASE = 28 ....  "
    ! Values of functions are used as analytical values
    fncn=0.D0
    fncn=fncn0
    
    !Values of gradients are analytical values at boundary and for particles close to boundary
    dfncn=0.D0
    dfncn=dfncn0

    ! Determine function's laplacian  𝜵∙kk 𝜵(f_a) using USAW formulation
    lap_fncn=0.D0
    
    ! Boundary vales of fucntions are its analytical values at the boundary
    do s= 1, etotal
       dfncn_s(:,s)=dfncn(:,nedge_rel_edge(s))
    enddo
    
    ! Determine function's laplacian  𝜵∙kk 𝜵(f_a) using USAW formulation    
    ! we set kk as 1 in 𝜵∙kk 𝜵(f_a)
    kk=1.D0
    
    call LaplacianSPHTraditional(lap_fncn, fncn, kk, x(1:SPH_dim,1:ntotal), &
                 & gamma_discrt, SPH_dim, niac, pair_i, pair_j, dwdx, mass, rho, itype, ntotal)

    
    allocate(bdryVal(SPH_dim,etotal))
    bdryVal=0.D0
    call BIL1(lap_fncn, dfncn, gamma_discrt, SPH_dim, eniac, epair_a, epair_s, &
                        & del_gamma_as, etype, etotal, ntotal, etype_periodic, &
                        & etype_SolidWall2,etype_SolidWall1,bdryVal, 1.D0)
    call BIL2(lap_fncn, dfncn_s, gamma_discrt, SPH_dim, eniac, epair_a, epair_s, &
                        & del_gamma_as, etype, etotal, ntotal, etype_periodic, &
                        & etype_SolidWall2,etype_SolidWall1,bdryVal,etotal, 1.D0)
    deallocate(bdryval)
    
    !relative error
    call relErrorFormulation(lap_fncn0,lap_fncn,RE(:,28),REB(:,28),REI(:,28),relErrorLogical)
    
    ! Output of the KCBI-CSPM formulation
    call outputFncnApprx(int(28),x, fncn0, fncn, dfncn0, dfncn, lap_fncn0, lap_fncn, nedge_rel_edge(1:etotal),edge(SPH_dim,1:etotal), etotal,SPH_dim,ntotal,itype(1:ntotal),etype(1:etotal),dataOutputPath)
   
 !---------------------------------------------------- BIL-USAW Formulation corrected with gamma_discrete,  inside only------------------------------------------------   
    write(*,*) "********************************************************************************************************* "
    write(*,*) "BIL-USAW Formulation BIL CASE = 29 ....  "
    ! Values of functions are used as analytical values
    fncn=0.D0
    fncn=fncn0
    
    !Values of gradients are analytical values at boundary and for particles close to boundary
    dfncn=0.D0
    dfncn=dfncn0

    ! Determine function's laplacian  𝜵∙kk 𝜵(f_a) using USAW formulation
    lap_fncn=0.D0
    
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
    
    ! Boundary vales of fucntions are its analytical values at the boundary
    do s= 1, etotal
       dfncn_s(:,s)=dfncn(:,nedge_rel_edge(s))
    enddo
    
    ! Determine function's laplacian  𝜵∙kk 𝜵(f_a) using USAW formulation    
    ! we set kk as 1 in 𝜵∙kk 𝜵(f_a)
    kk=1.D0
    
    call LaplacianSPHTraditional(lap_fncn, fncn, kk, x(1:SPH_dim,1:ntotal), &
                 & gamma_approx, SPH_dim, niac, pair_i, pair_j, dwdx, mass, rho, itype, ntotal)
    
    allocate(bdryVal(SPH_dim,etotal))
    bdryVal=0.D0
    call BIL1(lap_fncn, dfncn, gamma_cont, SPH_dim, eniac, epair_a, epair_s, &
                        & del_gamma_as, etype, etotal, ntotal, etype_periodic, &
                        & etype_SolidWall2,etype_SolidWall1,bdryVal, 1.D0)
    call BIL2(lap_fncn, dfncn_s, gamma_cont, SPH_dim, eniac, epair_a, epair_s, &
                        & del_gamma_as, etype, etotal, ntotal, etype_periodic, &
                        & etype_SolidWall2,etype_SolidWall1,bdryVal,etotal, 1.D0)
    
    deallocate(bdryval,gamma_approx)
    
    !relative error
    call relErrorFormulation(lap_fncn0,lap_fncn,RE(:,29),REB(:,29),REI(:,29),relErrorLogical)
    
    ! Output of the KCBI-CSPM formulation
    call outputFncnApprx(int(29),x, fncn0, fncn, dfncn0, dfncn, lap_fncn0, lap_fncn, nedge_rel_edge(1:etotal),edge(SPH_dim,1:etotal), etotal,SPH_dim,ntotal,itype(1:ntotal),etype(1:etotal),dataOutputPath)

!---------------------------------------------------- BIL-USAW Formulation corrected with gamma_discrete, but gamma_cont for boundary integral------------------------------------------------   
    write(*,*) "********************************************************************************************************* "
    write(*,*) "BIL-USAW Formulation BIL CASE = 30 ....  "
    ! Values of functions are used as analytical values
    fncn=0.D0
    fncn=fncn0
    
    !Values of gradients are analytical values at boundary and for particles close to boundary
    dfncn=0.D0
    dfncn=dfncn0

    ! Determine function's laplacian  𝜵∙kk 𝜵(f_a) using USAW formulation
    lap_fncn=0.D0
    
    ! Boundary vales of fucntions are its analytical values at the boundary
    do s= 1, etotal
       dfncn_s(:,s)=dfncn(:,nedge_rel_edge(s))
    enddo
    
    ! Determine function's laplacian  𝜵∙kk 𝜵(f_a) using USAW formulation    
    ! we set kk as 1 in 𝜵∙kk 𝜵(f_a)
    kk=1.D0
    
    call LaplacianSPHTraditional(lap_fncn, fncn, kk, x(1:SPH_dim,1:ntotal), &
                 & gamma_discrt, SPH_dim, niac, pair_i, pair_j, dwdx, mass, rho, itype, ntotal)
    
    allocate(bdryVal(SPH_dim,etotal))
    bdryVal=0.D0
    call BIL1(lap_fncn, dfncn, gamma_cont, SPH_dim, eniac, epair_a, epair_s, &
                        & del_gamma_as, etype, etotal, ntotal, etype_periodic, &
                        & etype_SolidWall2,etype_SolidWall1,bdryVal, 1.D0)
    call BIL2(lap_fncn, dfncn_s, gamma_cont, SPH_dim, eniac, epair_a, epair_s, &
                        & del_gamma_as, etype, etotal, ntotal, etype_periodic, &
                        & etype_SolidWall2,etype_SolidWall1,bdryVal,etotal, 1.D0)
    deallocate(bdryval)
    
    !relative error
    call relErrorFormulation(lap_fncn0,lap_fncn,RE(:,30),REB(:,30),REI(:,30),relErrorLogical)
    
    ! Output of the KCBI-CSPM formulation
    call outputFncnApprx(int(30),x, fncn0, fncn, dfncn0, dfncn, lap_fncn0, lap_fncn, nedge_rel_edge(1:etotal),edge(SPH_dim,1:etotal), etotal,SPH_dim,ntotal,itype(1:ntotal),etype(1:etotal),dataOutputPath)
   
   

 !---------------------------------------------------------------------------------------------------------------------------
    
 ! Output Relative Errors (L2 norm) and Max error (L inf)
    call outputError(int(30),RE,REB,REI,dlog10(dx_r),dataOutputPath,nreal,hsml_factor_name,hsml_factor_code)
    
    
deallocate(fncn0,fncn,fncn_s0)
deallocate(dfncn0,dfncn,dfncn_s)
deallocate(lap_fncn0,lap_fncn)
    
  end       
