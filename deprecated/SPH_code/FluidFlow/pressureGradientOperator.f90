!****************************************************************************
!
!  SUBROUTINE: pressureDensityOperattor
!
!  PURPOSE:  This subroutine determines pressure gradient   
!
!   CREATED:        08/21/2022       by  PARIKSHIT BOREGOWDA
!   Last Modified:  01/13/2023       by  PARIKSHIT BOREGOWDA 
!****************************************************************************
subroutine pressureGradientOperator(dstress,fncn,rho,oprtrTyype)
    use config_parameter, only: SPH_dim
    use particle_data, only: pair_i,pair_j,niac,ntotal,itype, dwdx, &
        & epair_a, epair_s, eniac, etotal, etype, del_gamma_as,nedge_rel_edge, &
        & gamma_cont,gamma_discrt, gamma_mat, gamma_mat_inv,del_gamma_as,xi1_mat_inv,xi_cont_mat_inv, &
        & mass

    implicit none
    
    real(8) :: fncn(ntotal), dstress(SPH_dim,ntotal), rho(ntotal)
    integer(4) :: oprtrTyype
    integer(4) a, k, d, s
    real(8),DIMENSION(:),ALLOCATABLE:: fncn_s, rho_s
    real(8),DIMENSION(:,:),ALLOCATABLE:: dfncn
    real(8),DIMENSION(:,:,:),ALLOCATABLE:: gamma_mat_approx_inv
    
    ALLOCATE(fncn_s(etotal),dfncn(SPH_dim,ntotal))
    fncn_s=0.D0
    dfncn=0.D0
    
    ! This evaluates function value at the boundary. This can be further generalized
    do s= 1, etotal
       fncn_s(s)=fncn(nedge_rel_edge(s))
    enddo
    
    if(oprtrTyype .eq. 1) then    
        !---------------------------------------------------- KCBI Formulation ------------------------------------------------   
        ! Determine function's first order derivative approximation using KCBI formulation
        call GradientKCBI(dfncn, fncn, fncn_s, gamma_cont, SPH_dim, niac, pair_i, pair_j, dwdx, eniac, epair_a, epair_s, del_gamma_as, mass, rho, itype, ntotal, etotal)
    
    elseif(oprtrTyype .eq. 2) then
        !---------------------------------------------------- PCBI Formulation ------------------------------------------------   
        ! Determine function's first order derivative approximation using PCBI formulation
        call GradientPCBI( dfncn, fncn, fncn_s, gamma_mat_inv, SPH_dim, niac, pair_i, pair_j, dwdx, eniac, epair_a, epair_s, del_gamma_as, mass, rho, itype, ntotal, etotal)

     
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
    
        ! Determine function's first order derivative approximation using WCBI1 formulation
        call GradientPCBI( dfncn, fncn, fncn_s, gamma_mat_approx_inv, SPH_dim, niac, pair_i, pair_j, dwdx, eniac, epair_a, epair_s, del_gamma_as, mass, rho, itype, ntotal, etotal)

        DEALLOCATE( gamma_mat_approx_inv )
    elseif(oprtrTyype .eq. 6) then
    !---------------------------------------------------- CSPM Formulation ------------------------------------------------   

        ! Determine function's first order derivative approximation using CSPM formulation
        call GradientCSPM( dfncn, fncn, xi1_mat_inv, SPH_dim, niac, pair_i, pair_j, dwdx, mass, rho, itype, ntotal)

    
    elseif(oprtrTyype .eq. 7) then
     !---------------------------------------------------- KCBI-CSPM Formulation ------------------------------------------------   
    
        ! Determine function's first order derivative approximation using KCBI-CSPM formulation
        call GradientKCBI_CSPM(dfncn, fncn, xi_cont_mat_inv, SPH_dim, niac, pair_i, pair_j, dwdx, eniac, epair_a, del_gamma_as, mass, rho, itype, ntotal, etotal)
    
    elseif(oprtrTyype .eq. 8) then
        !---------------------------------------------------- Ferrannd Gradient Formulation ------------------------------------------------   
        ALLOCATE(rho_s(etotal))
        rho_s=0.D0
         do s= 1, etotal
            rho_s(s)=rho(nedge_rel_edge(s))
        enddo
        ! Determine function's first order derivative approximation using PCBI formulation
        call GradientFerrand( dfncn, fncn, fncn_s, gamma_cont, SPH_dim, niac, pair_i, pair_j, dwdx, eniac, epair_a, epair_s, del_gamma_as, mass, rho,rho_s, itype, ntotal, etotal)

        DEALLOCATE(rho_s)
        
    elseif(oprtrTyype .eq. 9) then
        !---------------------------------------------------- symmetric-Conservative Inside, with PCBI bdry------------------------------------------------   
        ALLOCATE(rho_s(etotal))
        rho_s=0.D0
         do s= 1, etotal
            rho_s(s)=rho(nedge_rel_edge(s))
        enddo
        ! Determine function's first order derivative approximation using PCBI formulation
        call GradientPCBIbdrySymIns( dfncn, fncn, fncn_s, gamma_cont,gamma_mat_inv, SPH_dim, niac, pair_i, pair_j, dwdx, eniac, epair_a, epair_s, del_gamma_as, mass, rho,rho_s, itype, ntotal, etotal)
        DEALLOCATE(rho_s)
    
    elseif(oprtrTyype .eq. 10) then
        !---------------------------------------------------- asymmertic-Conistent Inside, with PCBI bdry ------------------------------------------------   
        
        ALLOCATE(rho_s(etotal))
        rho_s=0.D0
         do s= 1, etotal
            rho_s(s)=rho(nedge_rel_edge(s))
        enddo
        ! Use PCBI correction matrix, for Ferrand form instead of continuous_gamma
        call GradientPCBIbdryAssymIns( dfncn, fncn, fncn_s, gamma_cont,gamma_mat_inv, SPH_dim, niac, pair_i, pair_j, dwdx, eniac, epair_a, epair_s, del_gamma_as, mass, rho,rho_s, itype, ntotal, etotal)
        DEALLOCATE(rho_s)
        
    elseif(oprtrTyype .eq. 11) then
        !---------------------------------------------------- Kulasegaram Pressure Gradient Formulation ------------------------------------------------   
        
        ALLOCATE(rho_s(etotal))
        rho_s=0.D0
         do s= 1, etotal
            rho_s(s)=rho(nedge_rel_edge(s))
        enddo
        ! Determine function's first order derivative approximation using PCBI formulation
        call GradientKulasegaram( dfncn, fncn, fncn_s, gamma_cont, SPH_dim, niac, pair_i, pair_j, dwdx, eniac, epair_a, epair_s, del_gamma_as, mass, rho,rho_s, itype, ntotal, etotal)

        DEALLOCATE(rho_s)
        
        
     elseif(oprtrTyype .eq. 12) then
        !---------------------------------------------------- Kulasegaram Pressure Gradient Formulation ------------------------------------------------   
        
        ALLOCATE(rho_s(etotal))
        rho_s=0.D0
         do s= 1, etotal
            rho_s(s)=rho(nedge_rel_edge(s))
        enddo
        ! Determine function's first order derivative approximation using PCBI formulation
        call GradientKCKulasegaram( dfncn, fncn, fncn_s, gamma_cont, SPH_dim, niac, pair_i, pair_j, dwdx, eniac, epair_a, epair_s, del_gamma_as, mass, rho,rho_s, itype, ntotal, etotal)

        DEALLOCATE(rho_s)
        
        
     elseif(oprtrTyype .eq. 13) then
        !---------------------------------------------------- asymmertic-Conservative Inside, with WCBI1 bdry ------------------------------------------------   
        
        ALLOCATE( gamma_mat_approx_inv(SPH_dim,SPH_dim,ntotal)) 
        gamma_mat_approx_inv=0.D0
        do a=1,ntotal
            do d=1,SPH_dim
                    gamma_mat_approx_inv(d,d,a)= 1.D0/gamma_cont(a)
            enddo
        enddo

        ALLOCATE(rho_s(etotal))
        rho_s=0.D0
         do s= 1, etotal
            rho_s(s)=rho(nedge_rel_edge(s))
        enddo
        ! Use PCBI correction matrix, for Ferrand form instead of continuous_gamma
        call GradientPCBIbdrySymIns( dfncn, fncn, fncn_s, gamma_cont,gamma_mat_approx_inv, SPH_dim, niac, pair_i, pair_j, dwdx, eniac, epair_a, epair_s, del_gamma_as, mass, rho,rho_s, itype, ntotal, etotal)
        
        DEALLOCATE(rho_s,gamma_mat_approx_inv)
        
     elseif(oprtrTyype .eq. 14) then
        !---------------------------------------------------- Ferrannd Gradient Formulation ------------------------------------------------   
        ALLOCATE(rho_s(etotal))
        rho_s=0.D0
         do s= 1, etotal
            rho_s(s)=rho(nedge_rel_edge(s))
         enddo
         
        ! Set boundary pressure as zero to be used in Ferrand gradient, so boundary value can be included         
         fncn_s=0.D0
         
        ! Determine function's first order derivative approximation using PCBI formulation
        call GradientFerrand( dfncn, fncn, fncn_s, gamma_cont, SPH_dim, niac, pair_i, pair_j, dwdx, eniac, epair_a, epair_s, del_gamma_as, mass, rho,rho_s, itype, ntotal, etotal)

        call bdryPressureGradient( dfncn, fncn, gamma_cont)
        
        DEALLOCATE(rho_s)
        
    elseif(oprtrTyype .eq. 15) then
        !---------------------------------------------------- chiron Gradient Formulation ------------------------------------------------   
        ALLOCATE(rho_s(etotal))
        rho_s=0.D0
         do s= 1, etotal
            rho_s(s)=rho(nedge_rel_edge(s))
         enddo
         
        ! Set boundary pressure as zero to be used in Ferrand gradient, so boundary value can be included         
         fncn_s=0.D0
         
        ! Determine function's first order derivative approximation using PCBI formulation
        call GradientChiron( dfncn, fncn, fncn_s, gamma_cont, SPH_dim, niac, pair_i, pair_j, dwdx, eniac, epair_a, epair_s, del_gamma_as, mass, rho,rho_s, itype, ntotal, etotal)

        call bdryPressureGradient( dfncn, fncn, gamma_cont)
        
        DEALLOCATE(rho_s)
    endif
    
!---------------------------------------------------------------------------------------------------------------------------
    
    ! Now store  −𝜕_𝑖 𝑝 
    do d=1,SPH_dim
        dstress(d,1:ntotal)= dstress(d,1:ntotal) - dfncn(d,1:ntotal)
    enddo
    
    
deallocate(fncn_s,dfncn)
    
end    
