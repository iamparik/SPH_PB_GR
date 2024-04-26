!****************************************************************************
!
!  SUBROUTINE: EulerIntegration_fluid
!
!  PURPOSE:  Subroutine to march in time for fluid flow problem
!
!   CREATED:        08/20/2022       by  PARIKSHIT BOREGOWDA
!   Last Modified:  07/11/2023       by  PARIKSHIT BOREGOWDA 
!****************************************************************************

subroutine EulerIntegration_fluid(itimestep,dt)

    use config_parameter, only: SPH_dim, itype_real_max, itype_real_min, &
        & ConDivtype,PSTtype, SumDenstype, summationDensity, MLS_density_bound,&
        & MLS_step, densDiffType, delta_SPH, rho_init, HG_density_correction, save_step, g_const
    use particle_data, only: ntotal, etotal, itype,&
        & p,rho, rho_prev,drho, temp, dtemp, temp_prev,  &
        & x,vx, dstress, &
        &  rho, gamma_mat_inv, pBC_edges, b_MLS, &
        &  gamma_discrt, gamma_cont, niac, pair_i, pair_j, w, w_aa, mass,&
        &  max_vel, KE,PE, TE, delC, delCAvg, delCMax, del_gamma, delCL2, &
        & ga_Max,ga_Avg,ga_L2, g_a_min, g_a_max, FreeSurfaceVar, gamma_density_cont

    
    implicit none
    
    real(8), intent(in):: dt
    integer(4), intent(in)::itimestep
    integer(4) :: a,b, i, j , energyStepT
    logical:: AnalyticalSimulation= .false.
    real(8),DIMENSION(:,:),ALLOCATABLE :: vxCorrected
    real(8), dimension(:), ALLOCATABLE :: rho_temp
    real(8) dgammaMAX
    
     
     
    !Either run the Analytical Simulation/SPHSimulation
    if(AnalyticalSimulation) then
        !if runnign Analytical Simulation, make sure the solution for the given case exists and 
        ! is included in the file 
        
        ! The last input defines the analytical case icnluded in the code
        call AnalyticalFluidFlow(itimestep,dt,vx(1,1:ntotal),AnalyticalSimulation)
        
    else
        
        ! if there was periodic BC applied, for periodic particles and points update variables used 
        if (Allocated(pBC_edges)) then
            call PeriodicParameter(vx(1,:))
            call PeriodicParameter(vx(2,:))
            call PeriodicParameter(rho)
            call PeriodicParameter(p)
        endif
        
!------------------------DENSITY UPDATE -----------------------------!
        
        if(summationDensity) then

            call summationDensityOperator(rho,SumDenstype)
            
                    
            if (MLS_density_bound .and. (mod(itimestep,MLS_step) .eq. 0) ) then
            
                allocate(rho_temp(ntotal))
            
                call MLSfunctionApprox(rho_temp, rho, b_MLS,x, niac, pair_i, pair_j, w, w_aa, mass, rho, itype, gamma_cont, ntotal)
 
                rho(1:ntotal)=rho_temp(1:ntotal)

                deallocate(rho_temp)
            endif

        else         
                       
            allocate(drho(ntotal),rho_prev(ntotal))
            drho=0.D0
            
            rho_prev(1:ntotal)=rho(1:ntotal)

            !Update Density of particles using Continuity Equation
            call continuityDensityOperator(drho,vx,rho_prev, ConDivtype)
            
            !Use density diffusion term (δ-SPH)
            if (densDiffType .ne. 0) call DensityDiffusionOperator(drho, rho_prev, densDiffType, delta_SPH)

            
            if (itimestep .eq. 1) then
                do a=1,ntotal
                    if((itype(a) .le. itype_real_max) .and. (itype(a) .gt. itype_real_min)) then
                        rho(a)=rho_prev(a)
                    endif      
                enddo
            else
                call FreeSurfaceDetection
                
                ! The below is used for continuity desity updates
                do a=1,ntotal
                    !drho(a)= - drho(a)*rho_prev(a)
                    if((itype(a) .le. itype_real_max) .and. (itype(a) .gt. itype_real_min)) then                        
                        rho(a)=rho_prev(a) + dt*drho(a)
                        
                        !Below to have all free surf particles except those next to wall have 0 pressure boundary
                        !if(mod(itimestep,10) .eq. 0) then
                        !    if((FreeSurfaceVar(a) .le. 1.5D0) .and. (gamma_cont(a) .eq. 1.D0)) rho(a) =rho_init
                        !endif
                        
                        !below is supposidly hughes and grammes condition
                        if (HG_density_correction) then
                            if(rho(a) .le. rho_init) rho(a)=rho_init
                        endif
                    endif      
                enddo
                
              deallocate(FreeSurfaceVar)  
            endif
            
                    
            if (MLS_density_bound .and. (mod(itimestep,MLS_step) .eq. 0) ) then
            
                allocate(rho_temp(ntotal))
            
                call MLSfunctionApprox(rho_temp, rho, b_MLS,x, niac, pair_i, pair_j, w, w_aa, mass, rho, itype, gamma_cont, ntotal)

                !call MLSfunctionApprox(rho_temp, rho, b_MLS,x, niac, pair_i, pair_j, w, w_aa, mass, rho_prev, itype, gamma_cont, ntotal)

                rho(1:ntotal)=rho_temp(1:ntotal)

                deallocate(rho_temp)
            endif
            
            deallocate(drho,rho_prev)
            
            call gamma_density_continuous_leroy            
            !do a=1,ntotal
            !    if((itype(a) .le. itype_real_max) .and. (itype(a) .gt. itype_real_min)) then
            !        rho(a)=rho(a)/gamma_density_cont(a)
            !    endif      
            !enddo
            
        endif
        

        
        !Update desnity of periodic particles, if problem contians periodicBC
        if (Allocated(pBC_edges)) then
            call PeriodicParameter(rho)
        endif
        
        !   Calcualte denisty value at the boundary, to use for pressure boundary calculation.
        !   Also update Boundary density for periodicBC within the subroutine
        call boundaryDensityUpdate
        
        !Update desnity of periodic particles, if problem contians periodicBC
        if (Allocated(pBC_edges)) then
            call PeriodicParameter(rho)
        endif
        
        

!------------------------Acceleration UPDATE -----------------------------!
    
        allocate(dstress(SPH_dim,ntotal) )
        dstress=0.D0
        
        call FluidMomentum(dstress, x, vx, rho, p, itype, ntotal)
        
        !--------------Velocity update-----------------------
        ! calculate velocity from the momentum equation
        do a=1,ntotal
            dstress(:,a)=dstress(:,a)/rho(a)
            if((itype(a) .le. itype_real_max) .and. (itype(a) .gt. itype_real_min)) then
                do i =1,SPH_dim
                    vx(i,a)=vx(i,a) + dt*(dstress(i,a))                      
                enddo
            endif
        enddo
        
        deallocate(dstress)
 !------------------------POSITION UPDATE -----------------------------!   
        ! Update position from Velocity
        do a=1,ntotal
            if((itype(a) .le. itype_real_max) .and. (itype(a) .gt. itype_real_min)) then
               x(:,a)= x(:,a) + dt*vx(:,a)
            endif
        enddo       
        
    endif        
!---------------------- free surface detection and PST algorithm -------------------------------------!
    call FreeSurfaceDetection
    if(PSTtype .gt. 0) then
        call ParticleShiftingTechnique
    endif
    
!------------ Update parameter evolution of a system ----------------------!
    if( .not. allocated(delC)) call ConcGradient  
    !    
    !call evolvingVariables(itimestep,save_step)

    !if ((mod(itimestep,save_step).eq.0) .or. (itimestep.eq.1)) call output_time_flow(itimestep,dt)    
    !deallocate(delC)
    
    !if(mod(itimestep,save_step).eq. 0) deallocate(KE,PE,TE,max_vel,delcAvg,delCMax,delCL2,ga_Max,ga_Avg,ga_L2, g_a_min, g_a_max)
    
    if ((mod(itimestep,save_step).eq.0) .or. (itimestep.eq.1)) call output_flow_simplified(itimestep,dt)   
    deallocate(delC)
    deallocate(FreeSurfaceVar)
end
     

    