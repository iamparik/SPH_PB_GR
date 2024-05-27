!****************************************************************************
!
!  SUBROUTINE: continuityDensityOperattor
!
!  PURPOSE:  This subroutine determines density using continuity equation   
!
!   CREATED:        08/21/2022       by  PARIKSHIT BOREGOWDA
!   Last Modified:  01/13/2023       by  PARIKSHIT BOREGOWDA 
!****************************************************************************
subroutine continuityDensityOperator(divfncn,fncn,rho_in,oprtrTyype)
    use config_parameter, only: SPH_dim
    use particle_data, only: pair_i,pair_j,niac,ntotal,itype, dwdx, &
        & epair_a, epair_s, eniac, etotal, etype, del_gamma_as, nedge_rel_edge, &
        & gamma_cont,gamma_discrt, gamma_mat, gamma_mat_inv,&
        & del_gamma_as,xi1_mat_inv,xi_cont_mat_inv, &
        & mass, FreeSurfaceVar, gamma_density_cont

    implicit none
    
    real(8) :: fncn(SPH_dim,ntotal), divfncn(ntotal), rho_in(ntotal)
    integer(4) :: oprtrTyype
    integer(4) a, k, d, s
    real(8),DIMENSION(:,:),ALLOCATABLE:: fncn_s
    real(8),DIMENSION(:,:,:),ALLOCATABLE:: gamma_mat_approx_inv
    
    ALLOCATE(fncn_s(SPH_dim,etotal))
    fncn_s=0.D0
    
    ! This evaluates function value at the boundary. This can be further generalized
    do s= 1, etotal
       fncn_s(:,s)=fncn(:,nedge_rel_edge(s))
    enddo
    
    if(oprtrTyype .eq. 1) then    
        !---------------------------------------------------- KCBI Formulation ------------------------------------------------   
        ! Determine function's first order derivative approximation using KCBI formulation
        call DivergenceKCBI( divfncn, fncn, fncn_s, gamma_cont, SPH_dim, niac, pair_i, pair_j, dwdx, eniac, epair_a, epair_s, del_gamma_as, mass, rho_in, itype, ntotal, etotal)
            
    elseif(oprtrTyype .eq. 2) then
        !---------------------------------------------------- PCBI Formulation ------------------------------------------------   
        ! Determine function's first order derivative approximation using PCBI formulation
        call DivergencePCBI( divfncn, fncn, fncn_s, gamma_mat_inv, SPH_dim, niac, pair_i, pair_j, dwdx, eniac, epair_a, epair_s, del_gamma_as, mass, rho_in, itype, ntotal, etotal)

     
    elseif(oprtrTyype .eq. 3) then
        !---------------------------------------------------- PCBI based weakly consistent Formulation I ------------------------------------------------   
       ! This is a weakly consistent formulation where [¯𝛤_𝑎^(−1) ]  ≈1/𝛾_𝑎   [𝐼], where 𝛾_𝑎 is continuous
        ALLOCATE( gamma_mat_approx_inv(SPH_dim,SPH_dim,ntotal)) 
        gamma_mat_approx_inv=0.D0
        do a=1,ntotal
            do d=1,SPH_dim
                    gamma_mat_approx_inv(d,d,a)= 1.D0/gamma_cont(a)
            enddo
        enddo
        
        ! Determine function's first order derivative approximation using WCBI1 formulation
        call DivergencePCBI( divfncn, fncn, fncn_s, gamma_mat_approx_inv, SPH_dim, niac, pair_i, pair_j, dwdx, eniac, epair_a, epair_s, del_gamma_as, mass, rho_in, itype, ntotal, etotal)

        DEALLOCATE( gamma_mat_approx_inv )
    
    elseif(oprtrTyype .eq. 4) then
    !---------------------------------------------------- PCBI based weakly consistent Formulation 2------------------------------------------------   
    ! This is a weakly consistent formulation where [¯𝛤_𝑎^(−1) ]  ≈1/𝛾_𝑎   [𝐼], where 𝛾_𝑎 is discrete
        ALLOCATE( gamma_mat_approx_inv(SPH_dim,SPH_dim,ntotal)) 
        gamma_mat_approx_inv=0.D0
        do a=1,ntotal
            do d=1,SPH_dim
                    gamma_mat_approx_inv(d,d,a)= 1.D0/gamma_discrt(a)
            enddo
        enddo
    
        ! Determine function's first order derivative approximation using WCBI2 formulation
        call DivergencePCBI( divfncn, fncn, fncn_s, gamma_mat_approx_inv, SPH_dim, niac, pair_i, pair_j, dwdx, eniac, epair_a, epair_s, del_gamma_as, mass, rho_in, itype, ntotal, etotal)

        DEALLOCATE( gamma_mat_approx_inv )
    
    
    elseif(oprtrTyype .eq. 5) then
    !---------------------------------------------------- PCBI based weakly consistent Formulation 3------------------------------------------------   
    ! This is a weakly consistent formulation where only diagonal terms of 𝛤_𝑎 are retained for simplifying the invere matrix    
        ALLOCATE( gamma_mat_approx_inv(SPH_dim,SPH_dim,ntotal)) 
        gamma_mat_approx_inv=0.D0
        do a=1,ntotal
            do d=1,SPH_dim
                    gamma_mat_approx_inv(d,d,a)= 1.D0/gamma_mat(d,d,a)
            enddo
        enddo
    
        ! Determine function's first order derivative approximation using WCBI1 formulation
        call DivergencePCBI( divfncn, fncn, fncn_s, gamma_mat_approx_inv, SPH_dim, niac, pair_i, pair_j, dwdx, eniac, epair_a, epair_s, del_gamma_as, mass, rho_in, itype, ntotal, etotal)

        DEALLOCATE( gamma_mat_approx_inv )
    elseif(oprtrTyype .eq. 6) then
    !---------------------------------------------------- CSPM Formulation ------------------------------------------------   

        ! Determine function's first order derivative approximation using CSPM formulation
        call DivergenceCSPM(divfncn, fncn, xi1_mat_inv, SPH_dim, niac, pair_i, pair_j, dwdx, mass, rho_in, itype, ntotal)

    
    elseif(oprtrTyype .eq. 7) then
     !---------------------------------------------------- KCBI-CSPM Formulation ------------------------------------------------   
    
        ! Determine function's first order derivative approximation using KCBI-CSPM formulation
        call DivergenceKCBI_CSPM(divfncn, fncn, xi_cont_mat_inv, SPH_dim, niac, pair_i, pair_j, dwdx, eniac, epair_a, del_gamma_as, mass, rho_in, itype, ntotal, etotal)
    
    
    elseif(oprtrTyype .eq. 8) then
        !---------------------------------------------------- PCBI Formulation at wall boundary------------------------------------------------           
        ! Determine function's first order derivative approximation using PCBI formulation at the boundary but WCBI1 inside
        call DivergencePCBIBdry( divfncn, fncn, fncn_s, gamma_cont,gamma_mat_inv, SPH_dim, niac, pair_i, pair_j, dwdx, eniac, epair_a, epair_s, del_gamma_as, mass, rho_in, itype, ntotal, etotal)
   
    elseif(oprtrTyype .eq. 9) then
        !---------------------------------------------------- Summation density at wall boundary and continuity formulation everywhere else------------------------------------------------           
        ! Use summation fensity at wall boundary and WCBI1 inside
        call DivergenceSummationBdry( divfncn, fncn, rho_in)

    elseif(oprtrTyype .eq. 10) then
        
        !---------------------------------------------------- PCBI based weakly consistent Formulation I ------------------------------------------------   
       ! This is a weakly consistent formulation where [¯𝛤_𝑎^(−1) ]  ≈1/𝛾_𝑎   [𝐼], where 𝛾_𝑎 is continuous
        ALLOCATE( gamma_mat_approx_inv(SPH_dim,SPH_dim,ntotal)) 
        gamma_mat_approx_inv=0.D0
        do a=1,ntotal
            do d=1,SPH_dim
                    gamma_mat_approx_inv(d,d,a)= 1.D0/gamma_cont(a)
            enddo
        enddo
        
        ! Determine function's first order derivative approximation using WCBI1 formulation
        call DivergencePCBI( divfncn, fncn, fncn_s, gamma_mat_approx_inv, SPH_dim, niac, pair_i, pair_j, dwdx, eniac, epair_a, epair_s, del_gamma_as, mass, rho_in, itype, ntotal, etotal)

        DEALLOCATE( gamma_mat_approx_inv )    
    endif
    
!---------------------------------------------------------------------------------------------------------------------------

deallocate(fncn_s)


do a=1,ntotal
    divfncn(a)=-rho_in(a)*divfncn(a)
enddo

    
end    
