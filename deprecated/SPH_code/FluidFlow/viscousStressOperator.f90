!****************************************************************************
!
!  SUBROUTINE: viscousStressOperattor
!
!  PURPOSE:  This subroutine determines vscous Stress    
!
!   CREATED:        08/21/2022       by  PARIKSHIT BOREGOWDA
!   Last Modified:  01/28/2023       by  PARIKSHIT BOREGOWDA 
!****************************************************************************
subroutine ViscousStressOperator(dstress,fncn,rho,oprtrTyype)
    use config_parameter, only: SPH_dim, etype_periodic, etype_SolidWall1, & 
                & etype_SolidWall2, WallBoundaryLayer
    use particle_data, only: pair_i,pair_j,niac,ntotal,itype, dwdx, &
        & epair_a, epair_s, eniac, etotal, etype, del_gamma_as, nedge_rel_edge, surf_norm, &
        & gamma_cont,gamma_discrt, gamma_mat, gamma_mat_inv,del_gamma_as,xi1_mat_inv,xi_cont_mat_inv, &
        & mass, mu,x, pBC_edges !,bdryval_Vel

    implicit none
    
    real(8) :: dstress(SPH_dim,ntotal), fncn(SPH_dim,ntotal), rho(ntotal)
    integer(4) :: oprtrTyype
    integer(4) a, k, d, s, i,j
    real(8) G_s(etotal), G(ntotal) , bdryVal_vel(SPH_dim,etotal)
    real(8),DIMENSION(:),ALLOCATABLE:: gamma_approx, temp1D_1, temp1D_2, temp1D_3, temp1D_4
    real(8),DIMENSION(:,:,:),ALLOCATABLE:: gamma_mat_approx_inv
    real(8),DIMENSION(:,:),ALLOCATABLE:: fncn_s,dfncn_s,dfncn, dstrain, x_s, bdryValWallVel, temp2D, temp2D_2
    real(8),DIMENSION(:,:,:),ALLOCATABLE:: strain
    
    ALLOCATE(fncn_s(SPH_dim,etotal),strain(SPH_dim,SPH_dim,ntotal),dstrain(SPH_dim,ntotal),dfncn_s(SPH_dim,ntotal))
    fncn_s=0.D0
    strain=0.D0
    dstrain=0.D0
    dfncn_s=0.D0
    G=1.D0
    G_s= 1.D0
    ! This evaluates function value at the boundary. This can be further generalized
    do s= 1, etotal
       fncn_s(:,s)=fncn(:,nedge_rel_edge(s))
    enddo
    
    if(oprtrTyype .eq. 1) then    
        !---------------------------------------------------- KCBI Formulation ------------------------------------------------   
        ! Determine function's first order derivative approximation using KCBI formulation
         do i= 1,SPH_dim
            !Calculate 𝜀_𝑖𝑗=⟨⟨𝜕_𝑗 𝑣_𝑖 ⟩⟩
            call GradientKCBI( strain(i,:,:), fncn(i,:), fncn_s(i,:), gamma_cont, SPH_dim, niac, pair_i, pair_j, dwdx, eniac, epair_a, epair_s, del_gamma_as, mass, rho, itype, ntotal, etotal)

        enddo
        ! For periodic condition, update dfncn
        if (Allocated(pBC_edges)) then
            do i=1,SPH_dim
                do j=1,SPH_dim
                    call PeriodicParameter(strain(i,j,:))
                enddo
            enddo
        endif
        
        ! Determine the laplacian as ∇.∇v_i 
        do i=1,SPH_dim            
            ! This evaluates function first order derivative value at the boundary. This can be further generalized for numerical integration along edge...
            do s= 1, etotal
                dfncn_s(:,s)=strain(i,:,nedge_rel_edge(s))
            enddo
            
            call DivergenceKCBI( dstrain(i,:), strain(i,:,:), dfncn_s, gamma_cont, SPH_dim, niac, pair_i, pair_j, dwdx, eniac, epair_a, epair_s, del_gamma_as, mass, rho, itype, ntotal, etotal)
        
        enddo
        
        ! 𝜇 𝜕_𝑗 𝜕_𝑗 𝑣_𝑖
        do a=1,ntotal
            dstress(:,a)=dstress(:,a) + mu(a)*dstrain(:,a)    
        enddo

    elseif(oprtrTyype .eq. 2) then
        !---------------------------------------------------- PCBI Formulation ------------------------------------------------   
        ! Determine function's first order derivative approximation using PCBI formulation
        do i= 1,SPH_dim
            !Calculate 𝜀_𝑖𝑗=⟨⟨𝜕_𝑗 𝑣_𝑖 ⟩⟩
            call GradientPCBI( strain(i,:,:), fncn(i,:), fncn_s(i,:), gamma_mat_inv, SPH_dim, niac, pair_i, pair_j, dwdx, eniac, epair_a, epair_s, del_gamma_as, mass, rho, itype, ntotal, etotal)
            
        enddo
        ! For periodic condition, update dfncn
        if (Allocated(pBC_edges)) then
            do i=1,SPH_dim
                do j=1,SPH_dim
                    call PeriodicParameter(strain(i,j,:))
                enddo
            enddo
        endif
        
        ! Determine the laplacian as ∇.∇v_i 
        do i=1,SPH_dim            
            ! This evaluates function first order derivative value at the boundary. This can be further generalized for numerical integration along edge...
            do s= 1, etotal
                dfncn_s(:,s)=strain(i,:,nedge_rel_edge(s))
            enddo
            
            call DivergencePCBI( dstrain(i,:), strain(i,:,:), dfncn_s, gamma_mat_inv, SPH_dim, niac, pair_i, pair_j, dwdx, eniac, epair_a, epair_s, del_gamma_as, mass, rho, itype, ntotal, etotal)
        
        enddo
        
        ! 𝜇 𝜕_𝑗 𝜕_𝑗 𝑣_𝑖
        do a=1,ntotal
            dstress(:,a)=dstress(:,a) + mu(a)*dstrain(:,a)    
        enddo

     
    elseif(oprtrTyype .eq. 3) then
        !---------------------------------------------------- PCBI based weakly consistent Formulation I ------------------------------------------------   
       
        ALLOCATE( gamma_mat_approx_inv(SPH_dim,SPH_dim,ntotal)) 
        gamma_mat_approx_inv=0.D0
        do a=1,ntotal
            do d=1,SPH_dim
                    gamma_mat_approx_inv(d,d,a)= 1.D0/gamma_cont(a)
            enddo
        enddo
        
        ! Determine function's first order derivative approximation using PCBI formulation
        do i= 1,SPH_dim
            !Calculate 𝜀_𝑖𝑗=⟨⟨𝜕_𝑗 𝑣_𝑖 ⟩⟩
            call GradientPCBI( strain(i,:,:), fncn(i,:), fncn_s(i,:), gamma_mat_approx_inv, SPH_dim, niac, pair_i, pair_j, dwdx, eniac, epair_a, epair_s, del_gamma_as, mass, rho, itype, ntotal, etotal)
            
        enddo
        ! For periodic condition, update dfncn
        if (Allocated(pBC_edges)) then
            do i=1,SPH_dim
                do j=1,SPH_dim
                    call PeriodicParameter(strain(i,j,:))
                enddo
            enddo
        endif
        
        ! Determine the laplacian as ∇.∇v_i 
        do i=1,SPH_dim            
            ! This evaluates function first order derivative value at the boundary. This can be further generalized for numerical integration along edge...
            do s= 1, etotal
                dfncn_s(:,s)=strain(i,:,nedge_rel_edge(s))
            enddo
            
            call DivergencePCBI( dstrain(i,:), strain(i,:,:), dfncn_s, gamma_mat_approx_inv, SPH_dim, niac, pair_i, pair_j, dwdx, eniac, epair_a, epair_s, del_gamma_as, mass, rho, itype, ntotal, etotal)
        
        enddo
        
        ! 𝜇 𝜕_𝑗 𝜕_𝑗 𝑣_𝑖
        do a=1,ntotal
            dstress(:,a)=dstress(:,a) + mu(a)*dstrain(:,a)    
        enddo
        
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
    
        ! Determine function's first order derivative approximation using PCBI formulation
        do i= 1,SPH_dim
            !Calculate 𝜀_𝑖𝑗=⟨⟨𝜕_𝑗 𝑣_𝑖 ⟩⟩
            call GradientPCBI( strain(i,:,:), fncn(i,:), fncn_s(i,:), gamma_mat_approx_inv, SPH_dim, niac, pair_i, pair_j, dwdx, eniac, epair_a, epair_s, del_gamma_as, mass, rho, itype, ntotal, etotal)
            
        enddo
        ! For periodic condition, update dfncn
        if (Allocated(pBC_edges)) then
            do i=1,SPH_dim
                do j=1,SPH_dim
                    call PeriodicParameter(strain(i,j,:))
                enddo
            enddo
        endif
        
        ! Determine the laplacian as ∇.∇v_i 
        do i=1,SPH_dim            
            ! This evaluates function first order derivative value at the boundary. This can be further generalized for numerical integration along edge...
            do s= 1, etotal
                dfncn_s(:,s)=strain(i,:,nedge_rel_edge(s))
            enddo
            
            call DivergencePCBI( dstrain(i,:), strain(i,:,:), dfncn_s, gamma_mat_approx_inv, SPH_dim, niac, pair_i, pair_j, dwdx, eniac, epair_a, epair_s, del_gamma_as, mass, rho, itype, ntotal, etotal)
        enddo
        
        ! 𝜇 𝜕_𝑗 𝜕_𝑗 𝑣_𝑖
        do a=1,ntotal
            dstress(:,a)=dstress(:,a) + mu(a)*dstrain(:,a)    
        enddo
        
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
    
        ! Determine function's first order derivative approximation using PCBI formulation
        do i= 1,SPH_dim
            !Calculate 𝜀_𝑖𝑗=⟨⟨𝜕_𝑗 𝑣_𝑖 ⟩⟩
            call GradientPCBI( strain(i,:,:), fncn(i,:), fncn_s(i,:), gamma_mat_approx_inv, SPH_dim, niac, pair_i, pair_j, dwdx, eniac, epair_a, epair_s, del_gamma_as, mass, rho, itype, ntotal, etotal)
            
        enddo
        ! For periodic condition, update dfncn
        if (Allocated(pBC_edges)) then
            do i=1,SPH_dim
                do j=1,SPH_dim
                    call PeriodicParameter(strain(i,j,:))
                enddo
            enddo
        endif
        
        ! Determine the laplacian as ∇.∇v_i 
        do i=1,SPH_dim            
            ! This evaluates function first order derivative value at the boundary. This can be further generalized for numerical integration along edge...
            do s= 1, etotal
                dfncn_s(:,s)=strain(i,:,nedge_rel_edge(s))
            enddo
            
            call DivergencePCBI( dstrain(i,:), strain(i,:,:), dfncn_s, gamma_mat_approx_inv, SPH_dim, niac, pair_i, pair_j, dwdx, eniac, epair_a, epair_s, del_gamma_as, mass, rho, itype, ntotal, etotal)
        enddo
        
        ! 𝜇 𝜕_𝑗 𝜕_𝑗 𝑣_𝑖
        do a=1,ntotal
            dstress(:,a)=dstress(:,a) + mu(a)*dstrain(:,a)    
        enddo
        
        DEALLOCATE( gamma_mat_approx_inv )
    elseif(oprtrTyype .eq. 6) then
    !---------------------------------------------------- CSPM Formulation ----------------------------------------------
        write(*,*) "operation not included in the subroutine"
        pause
    
    elseif(oprtrTyype .eq. 7) then
     !---------------------------------------------------- KCBI-CSPM Formulation ------------------------------------------------   
    
        write(*,*) "operation not included in the subroutine"
        pause
        
    elseif(oprtrTyype .eq. 8) then
        !---------------------------------------------------- BIL-PCG Formulation (using BIL1+BIL2)------------------------------------------------   
    
        ! Determine function's first order derivative approximation using PCBI formulation
        do i= 1,SPH_dim
            
            !Calculate 𝜀_𝑖𝑗=⟨⟨𝜕_𝑗 𝑣_𝑖 ⟩⟩
            call GradientPCBI( strain(i,:,:), fncn(i,:), fncn_s(i,:), gamma_mat_inv, SPH_dim, niac,pair_i, pair_j, dwdx, eniac, epair_a, epair_s, del_gamma_as, mass, rho, itype, ntotal, etotal)
            
        enddo
        ! For periodic condition, update dfncn
        if (Allocated(pBC_edges)) then
            do i=1,SPH_dim
                do j=1,SPH_dim
                    call PeriodicParameter(strain(i,j,:))
                enddo
            enddo
        endif
        
        allocate(dfncn(SPH_dim,ntotal),bdryValWallVel(SPH_dim,eniac))
        dfncn=0.D0
        bdryValWallVel=0.D0
        
        if(WallBoundaryLayer) call BdryLayerVel(SPH_dim,bdryValWallVel)

        ! Determine function's laplacian  𝜵∙k 𝜵(v_a) using USAW formulation        
        do i=1,SPH_dim                 

			if(WallBoundaryLayer) then
                ! This evaluates function first order derivative value at the boundary. This can be further generalized for numerical integration along edge...
                do s= 1, etotal
                    dfncn_s(:,s)=mu(nedge_rel_edge(s))* strain(i,:,nedge_rel_edge(s))
                enddo 
                
                do a=1,ntotal
                    dfncn(:,a)= mu(a)*strain(i,:,a)
                enddo
                
                call LaplacianSPHTraditional(dstrain(i,:), fncn(i,:), mu, x(1:SPH_dim,1:ntotal), &
                    & gamma_cont, SPH_dim, niac, pair_i, pair_j, dwdx, mass, rho, itype, ntotal)
                
                call BIL1(dstrain(i,:), dfncn, gamma_cont, SPH_dim, eniac, epair_a, epair_s, del_gamma_as,&
                    & etype, etotal, ntotal, etype_periodic, etype_SolidWall2, etype_SolidWall1, bdryVal_vel(i,:), 1.D0)
                
                call BIL2(dstrain(i,:), dfncn_s, gamma_cont, SPH_dim, eniac, epair_a, epair_s,del_gamma_as, &
                    & etype, etotal, ntotal, etype_periodic, etype_SolidWall2,etype_SolidWall1, bdryValWallVel(i,:),eniac, 1.D0)
            else
                do s= 1, etotal
                    dfncn_s(:,s)=strain(i,:,nedge_rel_edge(s))
                enddo 
                
                do a=1,ntotal
                    dfncn(:,a)= strain(i,:,a)
                enddo
                
                call LaplacianSPHTraditional(dstrain(i,:), fncn(i,:), G, x(1:SPH_dim,1:ntotal), &
                 & gamma_cont, SPH_dim, niac, pair_i, pair_j, dwdx, mass, rho, itype, ntotal)
                
                call BIL1(dstrain(i,:), dfncn, gamma_cont, SPH_dim, eniac, epair_a, epair_s, del_gamma_as,&
                & etype, etotal, ntotal, etype_periodic, etype_SolidWall2, etype_SolidWall1, bdryVal_vel(i,:), 1.D0)
                
                call BIL2(dstrain(i,:), dfncn_s, gamma_cont, SPH_dim, eniac, epair_a, epair_s,del_gamma_as, &
                    & etype, etotal, ntotal, etype_periodic, etype_SolidWall2,etype_SolidWall1, bdryVal_vel(i,:),etotal, 1.D0)
            
            endif
        enddo

        
        deallocate(dfncn,bdryValWallVel)
        
        if(WallBoundaryLayer) then
        ! 𝜇 𝜕_𝑗 𝜕_𝑗 𝑣_𝑖
            do a=1,ntotal
                dstress(:,a)=dstress(:,a) + dstrain(:,a)    
            enddo
        else
            do a=1,ntotal
                dstress(:,a)=dstress(:,a) + mu(a)*dstrain(:,a)    
            enddo
        endif
        
     elseif(oprtrTyype .eq. 9) then    
    !---------------------------------------------------- BIL-PCG formulation Formulation (using 2*BIL1)------------------------------------------------   

         !Determine function's first order derivative approximation using PCBI formulation
        do i= 1,SPH_dim
            !Calculate 𝜀_𝑖𝑗=⟨⟨𝜕_𝑗 𝑣_𝑖 ⟩⟩
            call GradientPCBI( strain(i,:,:), fncn(i,:), fncn_s(i,:), gamma_mat_inv, SPH_dim, &
                & niac, pair_i, pair_j, dwdx, eniac, epair_a, epair_s, del_gamma_as, mass, rho, itype, ntotal, etotal)
        enddo
        
        ! For periodic condition, update dfncn
        if (Allocated(pBC_edges)) then
            do i=1,SPH_dim
                do j=1,SPH_dim
                    call PeriodicParameter(strain(i,j,:))
                enddo
            enddo
        endif

        ! Determine function's laplacian  𝜵∙ 𝜵(v_a) using BIL1 formulation
        do i=1,SPH_dim    
            
            call LaplacianSPHTraditional(dstrain(i,:), fncn(i,:),G, x(1:SPH_dim,1:ntotal),&
                & gamma_cont, SPH_dim, niac, pair_i, pair_j, dwdx, mass, rho, itype, ntotal )
       
            call BIL1(dstrain(i,:), strain(i,:,:), gamma_cont, SPH_dim, eniac, epair_a, epair_s, del_gamma_as, etype, etotal, ntotal,&
                & etype_periodic, etype_SolidWall2, etype_SolidWall1, bdryVal_vel(i,:), 2.D0)

        enddo
        
        ! 𝜇 𝜕_𝑗 𝜕_𝑗 𝑣_𝑖
        do a=1,ntotal
            dstress(:,a)=dstress(:,a) + mu(a)*dstrain(:,a)    
        enddo
    
     elseif(oprtrTyype .eq. 10) then 
    !---------------------------------------------------- BIL-PCG/BIL-BLT formulation Formulation (using 2*BIL2 = BIL-PCG, wall boundary layer = BIL-BLT)------------------------------------------------   
        !Determine function's first order derivative approximation using PCBI formulation      
        do i= 1,SPH_dim
            !Calculate 𝜀_𝑖𝑗=⟨⟨𝜕_𝑗 𝑣_𝑖 ⟩⟩
            allocate(temp2D(SPH_dim,ntotal),temp1D_1(ntotal),temp1D_2(etotal))
            temp2D(:,:)=strain(i,:,:)
            temp1D_1=fncn(i,:)
            temp1D_2=fncn_s(i,:)
            call GradientPCBI( temp2D, temp1D_1, temp1D_2, gamma_mat_inv, SPH_dim, niac,pair_i, pair_j, dwdx, eniac, epair_a, epair_s, del_gamma_as, mass, rho, itype, ntotal, etotal)
            strain(i,:,:)=temp2D(:,:)
            deallocate(temp2D,temp1D_1,temp1D_2)
        enddo
        ! For periodic condition, update dfncn
        if (Allocated(pBC_edges)) then
            do i=1,SPH_dim
                do j=1,SPH_dim
                    call PeriodicParameter(strain(i,j,:))
                enddo
            enddo
        endif
        
        
        allocate(dfncn(SPH_dim,ntotal),bdryValWallVel(SPH_dim,eniac))
        dfncn=0.D0
        bdryValWallVel=0.D0
        
        if(WallBoundaryLayer) call BdryLayerVel(SPH_dim,bdryValWallVel)

        ! Determine function's laplacian  𝜵∙k 𝜵(v_a) using USAW formulation        
        do i=1,SPH_dim                 

			if(WallBoundaryLayer) then
                ! This evaluates function first order derivative value at the boundary. This can be further generalized for numerical integration along edge...
                do s= 1, etotal
                    dfncn_s(:,s)=mu(nedge_rel_edge(s))* strain(i,:,nedge_rel_edge(s))
                enddo 
                
                allocate(temp1D_1(ntotal), temp2D(SPH_dim,ntotal))
                temp1D_1=fncn(i,:)
                temp2D=x(:,1:ntotal)
                call LaplacianSPHTraditional(dstrain(i,:), temp1D_1, mu, temp2D, &
                 & gamma_cont, SPH_dim, niac, pair_i, pair_j, dwdx, mass, rho, itype, ntotal)       
                deallocate(temp1D_1,temp2D)
                
                allocate(temp1D_1(eniac))
                temp1D_1=bdryValWallVel(i,:)
                call BIL2(dstrain(i,:), dfncn_s, gamma_cont, SPH_dim, eniac, epair_a, epair_s,del_gamma_as, &
                    & etype, etotal, ntotal, etype_periodic, etype_SolidWall2,etype_SolidWall1, temp1D_1,eniac, 2.D0)                
                deallocate(temp1D_1)
                    
            else
                do s= 1, etotal
                    dfncn_s(:,s)=strain(i,:,nedge_rel_edge(s))
                enddo 
                
                !allocate(temp1D_1(ntotal),temp1D_2(ntotal),temp2D(SPH_dim,ntotal))
                allocate(temp1D_1(ntotal), temp2D(SPH_dim,ntotal))
                temp1D_1=fncn(i,:)
                temp2D=x(:,1:ntotal)
                call LaplacianSPHTraditional(dstrain(i,:), temp1D_1, G, temp2D, &
                 & gamma_cont, SPH_dim, niac, pair_i, pair_j, dwdx, mass, rho, itype, ntotal)
                deallocate(temp1D_1,temp2D)
                
                allocate(temp1D_1(eniac))
                temp1D_1=bdryValWallVel(i,:)
				call BIL2(dstrain(i,:), dfncn_s, gamma_cont, SPH_dim, eniac, epair_a, epair_s,del_gamma_as, &
                    & etype, etotal, ntotal, etype_periodic, etype_SolidWall2,etype_SolidWall1, temp1D_1,etotal, 2.D0)
                deallocate(temp1D_1)
            endif
        enddo

        
        deallocate(dfncn,bdryValWallVel)
        
        if(WallBoundaryLayer) then
        ! 𝜇 𝜕_𝑗 𝜕_𝑗 𝑣_𝑖
            do a=1,ntotal
                dstress(:,a)=dstress(:,a) + dstrain(:,a)    
            enddo
        else
            do a=1,ntotal
                dstress(:,a)=dstress(:,a) + mu(a)*dstrain(:,a)    
                !dstress(2,a)=dstress(2,a) + mu(a)*dstrain(2,a) 
            enddo
        endif
        
      
    
     elseif(oprtrTyype .eq. 11) then 
    
    !---------------------------------------------------- BIL-Macia Formulation------------------------------------------------   
        
        ! Determine function's laplacian  𝜵∙ 𝜵(v_a) using BIL1 formulation
         
        ALLOCATE(x_s(SPH_dim,etotal))       
        x_s=0.D0
        
        do s= 1, etotal
            x_s(:,s)=x(:,nedge_rel_edge(s))
        enddo
        
        do i=1,SPH_dim    
            
            call LaplacianSPHTraditional(dstrain(i,:), fncn(i,:), G, x(1:SPH_dim,1:ntotal), &
                 & gamma_cont, SPH_dim, niac, pair_i, pair_j, dwdx, mass, rho, itype, ntotal)


            ! Boundary vales of fucntions are its analytical values at the boundary
            do s= 1, etotal
               G_s(s)=1.D0
            enddo
    
            call BILMacia(dstrain(i,:), fncn(i,:), fncn_s(i,:), x(1:SPH_dim,1:ntotal), x_s, gamma_cont, &
                & SPH_dim, eniac, epair_a, epair_s, del_gamma_as, etype, etotal, ntotal, G_s)
  
        enddo
  
        ! 𝜇 𝜕_𝑗 𝜕_𝑗 𝑣_𝑖
        do a=1,ntotal
            dstress(:,a)=dstress(:,a) + mu(a)*dstrain(:,a)    
        enddo
        
        deallocate(x_s)
        
   
     elseif(oprtrTyype .eq. 31) then
        !---------------------------------------------------- BIL-NTG Formulation (using continuous gamma)------------------------------------------------   
    
        ! Determine function's first order derivative approximation using PCBI formulation
        do i= 1,SPH_dim
            !Calculate 𝜀_𝑖𝑗=⟨⟨𝜕_𝑗 𝑣_𝑖 ⟩⟩
            allocate(temp2D(SPH_dim,ntotal),temp1D_1(ntotal),temp1D_2(etotal))
            temp2D(:,:)=strain(i,:,:)
            temp1D_1=fncn(i,:)
            temp1D_2=fncn_s(i,:)
            call GradientPCBI( temp2D, temp1D_1, temp1D_2, gamma_mat_inv, SPH_dim, niac,pair_i, pair_j, dwdx, eniac, epair_a, epair_s, del_gamma_as, mass, rho, itype, ntotal, etotal)
            strain(i,:,:)=temp2D(:,:)
            deallocate(temp2D,temp1D_1,temp1D_2)
        enddo
        ! For periodic condition, update dfncn
        if (Allocated(pBC_edges)) then
            do i=1,SPH_dim
                do j=1,SPH_dim
                    call PeriodicParameter(strain(i,j,:))
                enddo
            enddo
        endif
        
        allocate(dfncn(SPH_dim,ntotal),bdryValWallVel(SPH_dim,eniac), x_s(SPH_dim,etotal))
        
        x_s=0.D0
        dfncn=0.D0
        bdryValWallVel=0.D0
        
        do s= 1, etotal
            x_s(:,s)=x(:,nedge_rel_edge(s))
        enddo
        
        if(WallBoundaryLayer) call BdryLayerVel(SPH_dim,bdryValWallVel)

        ! Determine function's laplacian  𝜵∙k 𝜵(v_a) using USAW formulation        
        do i=1,SPH_dim                 

			if(WallBoundaryLayer) then

                do s= 1, etotal
                    G_s(s)=1.D0
                enddo 
                
                
                call LaplacianSPHTraditional(dstrain(i,:), fncn(i,:), mu, x(1:SPH_dim,1:ntotal), &
                 & gamma_cont, SPH_dim, niac, pair_i, pair_j, dwdx, mass, rho, itype, ntotal)

                call BILMEA(dstrain(i,:), fncn(i,:), fncn_s(i,:), x(1:SPH_dim,1:ntotal), x_s, gamma_cont, SPH_dim, eniac, &
                & epair_a, epair_s, del_gamma_as, etype, etotal, ntotal,  etype_periodic, etype_SolidWall2,etype_SolidWall1, &
                & bdryValWallVel(i,:), etotal, surf_norm, G, G_s)
                
            else
                
                do s= 1, etotal
                    G_s(s)=mu(nedge_rel_edge(s))
                enddo 
                
                allocate(temp1D_1(ntotal), temp2D(SPH_dim,ntotal))
                temp1D_1=fncn(i,:)
                temp2D=x(:,1:ntotal)                
                call LaplacianSPHTraditional(dstrain(i,:), temp1D_1, mu, temp2D, &
                 & gamma_cont, SPH_dim, niac, pair_i, pair_j, dwdx, mass, rho, itype, ntotal)
                deallocate(temp1D_1,temp2D)
                
                allocate(temp1D_1(ntotal), temp2D(SPH_dim,ntotal),temp1D_2(etotal),temp2D_2(SPH_dim,etotal),temp1D_4(eniac))
                temp1D_1=fncn(i,:)
                temp1D_2=fncn_s(i,:)
                temp2D=x(:,1:ntotal)   
                temp2D_2=surf_norm(:,1:etotal)
                temp1D_4=bdryVal_vel(i,:)
                call BILMEA(dstrain(i,:), temp1D_1, temp1D_2,temp2D, x_s, gamma_cont, SPH_dim, eniac, &
                & epair_a, epair_s, del_gamma_as, etype, etotal, ntotal,  etype_periodic, etype_SolidWall2,etype_SolidWall1, &
                & temp1D_4, etotal, temp2D_2, mu, G_s)
                deallocate(temp1D_1,temp2D, temp1D_2,temp2D_2, temp1D_4)

            endif
        enddo

        
        deallocate(dfncn,bdryValWallVel, x_s)
        
        do a=1,ntotal
            dstress(:,a)=dstress(:,a) + dstrain(:,a)    
        enddo
        
        
    endif   
    
       
    
!---------------------------------------------------------------------------------------------------------------------------

deallocate(dfncn_s,dstrain,strain, fncn_s)
    
end    
