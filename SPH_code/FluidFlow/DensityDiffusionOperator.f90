!****************************************************************************
!
!  SUBROUTINE: DensityDiffusionOperator
!
!  PURPOSE:  This subroutine determines density using continuity equation   
!
!   CREATED:        08/01/2023       by  PARIKSHIT BOREGOWDA
!   Last Modified:  08/01/2023       by  PARIKSHIT BOREGOWDA 
!****************************************************************************
subroutine DensityDiffusionOperator(drho, rho, densDiffType, delta_SPH)
    use config_parameter, only: SPH_dim, itype_real_min, &
        & itype_real_max, itype_virtual, hsml_const, c_sound
    use particle_data, only: pair_i,pair_j,niac,ntotal,itype, dwdx, &
        & epair_a, epair_s, eniac, etotal, etype, del_gamma_as,nedge_rel_edge, &
        & gamma_cont,gamma_discrt, gamma_mat, gamma_mat_inv,del_gamma_as,xi1_mat_inv,xi_cont_mat_inv, &
        & mass, x, pBC_edges


    implicit none
    
    real(8), intent(in) :: rho(ntotal), delta_SPH
    real(8), intent(inout) :: drho(ntotal)
    integer(4), intent(in):: densDiffType
    integer(4) a,b, k, d, s
    real(8) diff_term(SPH_dim)
    real(8),DIMENSION(:,:),ALLOCATABLE:: gradrho
    real(8),DIMENSION(:),ALLOCATABLE:: rho_s
    real(8),DIMENSION(:,:,:),ALLOCATABLE:: gamma_mat_approx_inv
    
    ALLOCATE(gradrho(SPH_dim,ntotal), rho_s(etotal))
    gradrho=0.D0
    rho_s=0.D0
    
    ! This evaluates function value at the boundary. This can be further generalized
    do s= 1, etotal
       rho_s(s)=rho(nedge_rel_edge(s))
    enddo
    
    
    
    if(densDiffType .eq. 1) then    
    !---------------------------------------------------- KCBI Formulation ------------------------------------------------   
    ! Determine function's first order derivative approximation using KCBI formulation
        call GradientKCBI(gradrho, rho, rho_s, gamma_cont, SPH_dim, niac, pair_i, pair_j, dwdx, eniac, epair_a, epair_s, del_gamma_as, mass, rho, itype, ntotal, etotal)
    
    elseif(densDiffType .eq. 2) then
    !---------------------------------------------------- PCBI Formulation ------------------------------------------------   
    ! Determine function's first order derivative approximation using PCBI formulation
        call GradientPCBI( gradrho, rho, rho_s, gamma_mat_inv, SPH_dim, niac, pair_i, pair_j, dwdx, eniac, epair_a, epair_s, del_gamma_as, mass, rho, itype, ntotal, etotal)

     
    elseif(densDiffType .eq. 3) then
    !---------------------------------------------------- PCBI based weakly consistent Formulation I ------------------------------------------------   
       
        ALLOCATE( gamma_mat_approx_inv(SPH_dim,SPH_dim,ntotal)) 
        gamma_mat_approx_inv=0.D0
        do a=1,ntotal
            do d=1,SPH_dim
                    gamma_mat_approx_inv(d,d,a)= 1.D0/gamma_cont(a)
            enddo
        enddo
        
        ! Determine function's first order derivative approximation using WCBI1 formulation
        call GradientPCBI( gradrho, rho, rho_s, gamma_mat_approx_inv, SPH_dim, niac, pair_i, pair_j, dwdx, eniac, epair_a, epair_s, del_gamma_as, mass, rho, itype, ntotal, etotal)

        DEALLOCATE( gamma_mat_approx_inv )
    
    elseif(densDiffType .eq. 4) then
    !---------------------------------------------------- PCBI based weakly consistent Formulation 2------------------------------------------------   
    
        ALLOCATE( gamma_mat_approx_inv(SPH_dim,SPH_dim,ntotal)) 
        gamma_mat_approx_inv=0.D0
        do a=1,ntotal
            do d=1,SPH_dim
                    gamma_mat_approx_inv(d,d,a)= 1.D0/gamma_discrt(a)
            enddo
        enddo
    
        ! Determine function's first order derivative approximation using WCBI2 formulation
        call GradientPCBI( gradrho, rho, rho_s, gamma_mat_approx_inv, SPH_dim, niac, pair_i, pair_j, dwdx, eniac, epair_a, epair_s, del_gamma_as, mass, rho, itype, ntotal, etotal)

        DEALLOCATE( gamma_mat_approx_inv )
    
    
    elseif(densDiffType .eq. 5) then
    !---------------------------------------------------- PCBI based weakly consistent Formulation 3------------------------------------------------   
        
        ALLOCATE( gamma_mat_approx_inv(SPH_dim,SPH_dim,ntotal)) 
        gamma_mat_approx_inv=0.D0
        do a=1,ntotal
            do d=1,SPH_dim
                    gamma_mat_approx_inv(d,d,a)= 1.D0/gamma_mat(d,d,a)
            enddo
        enddo
    
        ! Determine function's first order derivative approximation using WCBI1 formulation
        call GradientPCBI( gradrho, rho, rho_s, gamma_mat_approx_inv, SPH_dim, niac, pair_i, pair_j, dwdx, eniac, epair_a, epair_s, del_gamma_as, mass, rho, itype, ntotal, etotal)

        DEALLOCATE( gamma_mat_approx_inv )
        
    elseif(densDiffType .eq. 6) then
    !---------------------------------------------------- CSPM Formulation ------------------------------------------------   

    ! Determine function's first order derivative approximation using CSPM formulation
        call GradientCSPM( gradrho, rho, xi1_mat_inv, SPH_dim, niac, pair_i, pair_j, dwdx, mass, rho, itype, ntotal)

    
    elseif(densDiffType .eq. 7) then
    !---------------------------------------------------- KCBI-CSPM Formulation ------------------------------------------------   
    
        ! Determine function's first order derivative approximation using KCBI-CSPM formulation
        call GradientKCBI_CSPM(gradrho, rho, xi_cont_mat_inv, SPH_dim, niac, pair_i, pair_j, dwdx, eniac, epair_a, del_gamma_as, mass, rho, itype, ntotal, etotal)
    
        elseif(densDiffType .eq. 8) then
        !---------------------------------------------------- Ferrannd Gradient Formulation ------------------------------------------------   
        ALLOCATE(rho_s(etotal))
        rho_s=0.D0
            do s= 1, etotal
            rho_s(s)=rho(nedge_rel_edge(s))
        enddo
        ! Determine function's first order derivative approximation using PCBI formulation
        call GradientFerrand( gradrho, rho, rho_s, gamma_cont, SPH_dim, niac, pair_i, pair_j, dwdx, eniac, epair_a, epair_s, del_gamma_as, mass, rho,rho_s, itype, ntotal, etotal)

        DEALLOCATE(rho_s)
        
    elseif(densDiffType .eq. 9) then
    !---------------------------------------------------- symmetric-Conservative Inside, with PCBI bdry------------------------------------------------   
        ALLOCATE(rho_s(etotal))
        rho_s=0.D0
            do s= 1, etotal
            rho_s(s)=rho(nedge_rel_edge(s))
        enddo
        ! Determine function's first order derivative approximation using PCBI formulation
        call GradientPCBIbdrySymIns( gradrho, rho, rho_s, gamma_cont,gamma_mat_inv, SPH_dim, niac, pair_i, pair_j, dwdx, eniac, epair_a, epair_s, del_gamma_as, mass, rho,rho_s, itype, ntotal, etotal)

    
    elseif(densDiffType .eq. 10) then
    !---------------------------------------------------- asymmertic-Conistent Inside, with PCBI bdry ------------------------------------------------   
        
        ALLOCATE(rho_s(etotal))
        rho_s=0.D0
            do s= 1, etotal
            rho_s(s)=rho(nedge_rel_edge(s))
        enddo
        ! Use PCBI correction matrix, for Ferrand form instead of continuous_gamma
        call GradientPCBIbdryAssymIns( gradrho, rho, rho_s, gamma_cont,gamma_mat_inv, SPH_dim, niac, pair_i, pair_j, dwdx, eniac, epair_a, epair_s, del_gamma_as, mass, rho,rho_s, itype, ntotal, etotal)

        
    elseif(densDiffType .eq. 11) then
    !---------------------------------------------------- Kulasegaram Pressure Gradient Formulation ------------------------------------------------   
        
        ALLOCATE(rho_s(etotal))
        rho_s=0.D0
            do s= 1, etotal
            rho_s(s)=rho(nedge_rel_edge(s))
        enddo
        ! Determine function's first order derivative approximation using PCBI formulation
        call GradientKulasegaram( gradrho, rho, rho_s, gamma_cont, SPH_dim, niac, pair_i, pair_j, dwdx, eniac, epair_a, epair_s, del_gamma_as, mass, rho,rho_s, itype, ntotal, etotal)

        DEALLOCATE(rho_s)
        
        
    elseif(densDiffType .eq. 12) then
    !---------------------------------------------------- Kulasegaram Pressure Gradient Formulation ------------------------------------------------   
        
        ALLOCATE(rho_s(etotal))
        rho_s=0.D0
            do s= 1, etotal
            rho_s(s)=rho(nedge_rel_edge(s))
        enddo
        ! Determine function's first order derivative approximation using PCBI formulation
        call GradientKCKulasegaram( gradrho, rho, rho_s, gamma_cont, SPH_dim, niac, pair_i, pair_j, dwdx, eniac, epair_a, epair_s, del_gamma_as, mass, rho,rho_s, itype, ntotal, etotal)

        
        
    elseif(densDiffType .eq. 13) then
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
        call GradientPCBIbdrySymIns( gradrho, rho, rho_s, gamma_cont,gamma_mat_approx_inv, SPH_dim, niac, pair_i, pair_j, dwdx, eniac, epair_a, epair_s, del_gamma_as, mass, rho,rho_s, itype, ntotal, etotal)
        
        DEALLOCATE(gamma_mat_approx_inv)
    
     
    endif
     
    DEALLOCATE(rho_s)
!---------------------------------------------------------------------------------------------------------------------------

     ! if there was periodic Bcapplied, for periodic particles and points update variables used 
    if (Allocated(pBC_edges)) then
        call PeriodicParameter(gradrho(1,:))
        call PeriodicParameter(gradrho(2,:))
    endif
        

    do k=1,niac
        
        a=pair_i(k)
        b=pair_j(k)        

        if((mod(itype(a),itype_virtual) .le. itype_real_max) .and. (mod(itype(a),itype_virtual) .ge. itype_real_min)) then
            diff_term=2.D0*(rho(b)-rho(a))*(x(:,b)-x(:,a))/(norm2(x(:,b)-x(:,a))**2)- (gradrho(:,a)+gradrho(:,b))
            if(mass(b) .gt. 0.D0) then
                if(isNAN(norm2(diff_term))) then
                    write(*,*) "NAN encountered for a=",a," and b=",b
                else
                    drho(a)= drho(a) + delta_SPH*hsml_const*c_sound*dot_product(diff_term, dwdx(:,k))*mass(b)/rho(b)
                endif
                
            endif
            
        endif
        
        if((mod(itype(b),itype_virtual) .le. itype_real_max) .and. (mod(itype(b),itype_virtual) .ge. itype_real_min)) then
            diff_term=2.D0*(rho(a)-rho(b))*(x(:,a)-x(:,b))/(norm2(x(:,a)-x(:,b))**2)- (gradrho(:,b)+gradrho(:,a))
            if(mass(a) .gt. 0.D0) then
                if(isNAN(norm2(diff_term))) then
                    write(*,*) "NAN encountered for a=",a," and b=",b
                else
                    drho(b)= drho(b) + delta_SPH*hsml_const*c_sound*dot_product(diff_term, -dwdx(:,k))*mass(a)/rho(a)
                endif
                
            endif
            
        endif
        
    enddo
        
    
!---------------------------------------------------------------------------------------------------------------------------

    
    deallocate(gradrho)
end    
