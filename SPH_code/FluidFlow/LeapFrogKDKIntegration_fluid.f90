!****************************************************************************
!
!  SUBROUTINE: LeapFrogKDKIntegration_fluid
!
!  PURPOSE:  Subroutine to march in time for fluid flow problem
!
!   CREATED:        07/11/2023       by  PARIKSHIT BOREGOWDA
!   Last Modified:  07/11/2023       by  PARIKSHIT BOREGOWDA 
!****************************************************************************

subroutine LeapFrogKDKIntegration_fluid(itimestep,dt)

    use config_parameter, only: SPH_dim, itype_real_max, itype_real_min
    use config_geometry, only:ConDivtype,SumDenstype, summationDensity, save_step
    use particle_data, only: ntotal, etotal, itype, maxn, &
        & p,rho, rho_prev,drho, temp, dtemp, temp_prev,  &
        & x,vx, dstress, &
        &  rho, gamma_mat_inv, pBC_edges, gamma_cont, &
        &  gamma_discrt, niac, pair_i, pair_j, w, w_aa, mass
    
    
    implicit none
    
    
    real(8), intent(in):: dt
    integer(4), intent(in)::itimestep
    integer(4) :: a, i, j 
    logical:: AnalyticalSimulation= .false.
    real(8),DIMENSION(:,:),ALLOCATABLE :: vxCorrected
    real(8),DIMENSION(:),ALLOCATABLE :: rho_temp
    
    ! if there was periodic Bcapplied, for periodic particles and points update variables used 
    if (Allocated(pBC_edges)) then
        do i =1,SPH_dim
            call PeriodicParameter(vx(i,:))
        enddo
        call PeriodicParameter(rho)
        call PeriodicParameter(p)
    endif  
   
    
!------------------------HALF Step update of DENSITY and VELOCITY -----------------------------!
    
    if (itimestep .ne. 1) then       
        
        !update HALF step DENSITY
        if(summationDensity) then
            call summationDensityOperator(rho,SumDenstype)
        else
            if(ConDivtype .eq. 9) then
                allocate(rho_temp(maxn))
                call summationDensityOperator(rho_temp,SumDenstype)   
                do a=1,ntotal 
                    if((itype(a) .le. itype_real_max) .and. (itype(a) .gt. itype_real_min)) then
                        if((gamma_cont(a) .lt. 1.D0) .and. (gamma_cont(a) .gt. 0.D0)) then
                            rho_prev(a)=rho_temp(a)
                            drho(a)=0.D0
                        endif
                    endif
                enddo
                
                deallocate(rho_temp)
            endif
            
            !The below is used for continuity desity updates
            do a=1,ntotal
                if((itype(a) .le. itype_real_max) .and. (itype(a) .gt. itype_real_min)) then
                    rho(a)=rho_prev(a) + (dt/2.D0)*drho(a)
                endif      
            enddo 
            !deallocate(drho,rho_prev)
        endif
        
        !update HALF step VELOCITY
        do a=1,ntotal
            if((itype(a) .le. itype_real_max) .and. (itype(a) .gt. itype_real_min)) then
                do i =1,SPH_dim
                    vx(i,a)=vx(i,a) + (dt/2.D0)*(dstress(i,a)/rho(a))                      
                enddo
            endif
        enddo    
        
    endif   
        
    !Update desnity,velocity and pressure of periodic particles, if problem contians periodicBC
    if (Allocated(pBC_edges)) then
        call PeriodicParameter(rho)
        do i =1,SPH_dim
            call PeriodicParameter(vx(i,:))
        enddo
        call PeriodicParameter(p)
    endif
        
    !   Calcualte denisty value at the boundary, to use for pressure boundary calculation.
    !   Also update Boundary density for periodicBC within the subroutine
    call boundaryDensityUpdate
        
    !Update density of periodic particles, if problem contians periodicBC
    if (Allocated(pBC_edges)) then
        call PeriodicParameter(rho)
    endif

!------------------------Acceleration UPDATE -----------------------------!
    
    if(.not.(allocated(dstress))) allocate(dstress(SPH_dim,ntotal))
    dstress=0.D0
        
    call FluidMomentum(dstress, x, vx, rho, p, itype, ntotal)
    
    if(itimestep .eq. 1) then
        
        !--------------Velocity update for half a timestep-----------------------
        ! calculate velocity from the momentum equation
        do a=1,ntotal
            dstress(:,a)=dstress(:,a)/rho(a)
            if((itype(a) .le. itype_real_max) .and. (itype(a) .gt. itype_real_min)) then
                do i =1,SPH_dim
                    vx(i,a)=vx(i,a) + (dt/2.D0)*dstress(i,a)               
                enddo
            endif
        enddo
        
        !deallocate(dstress)
        
        !--------------Density update for half a timestep-----------------------
        if(.not. summationDensity) then
            
            allocate(drho(ntotal),rho_prev(ntotal))
            drho=0.D0
    
            !Update Density of particles using Continuity Equation
            call continuityDensityOperator(drho,vx,rho, ConDivtype)
            
        
            !The below is used for continuity desity updates
            do a=1,ntotal
                drho(a)= -drho(a)*rho(a)
                if((itype(a) .le. itype_real_max) .and. (itype(a) .gt. itype_real_min)) then                    
                    rho(a)=rho(a) + (dt/2.D0)*drho(a)
                endif      
                rho_prev(a)=rho(a)  
            enddo
            
            !deallocate(drho,rho_prev)
        endif
        
    else
        !--------------Velocity update for m-1/2 time step from m-3/2 time step-----------------------
        ! calculate velocity from the momentum equation
        do a=1,ntotal
            dstress(:,a)=dstress(:,a)/rho(a)
            if((itype(a) .le. itype_real_max) .and. (itype(a) .gt. itype_real_min)) then
                do i =1,SPH_dim
                    vx(i,a)=vx(i,a) + dt*dstress(i,a)                    
                enddo
            endif
        enddo     
        
        !--------------Density update for m-1/2 time step from m-3/2 time step-----------------------
        if(.not. summationDensity) then
            
            drho=0.D0
    
            !Update Density of particles using Continuity Equation
            call continuityDensityOperator(drho,vx,rho, ConDivtype)
            
            if(ConDivtype .eq. 9) then                
                do a=1,ntotal 
                    if((itype(a) .le. itype_real_max) .and. (itype(a) .gt. itype_real_min)) then
                        if((gamma_cont(a) .lt. 1.D0) .and. (gamma_cont(a) .gt. 0.D0)) then
                            rho_prev(a)=rho(a)
                        endif
                    endif
                enddo
            endif
        
            !The below is used for continuity desity updates
            do a=1,ntotal
                drho(a)= -drho(a)*rho(a)
                if((itype(a) .le. itype_real_max) .and. (itype(a) .gt. itype_real_min)) then                    
                    rho(a)=rho_prev(a) + dt*drho(a)
                endif      
                rho_prev(a)=rho(a)                
            enddo
            
            !deallocate(drho,rho_prev)
        endif
        
    endif
    
!------------------------POSITION UPDATE -----------------------------! 
    !𝒙_𝑎^𝑚=𝒙_𝑎^(𝑚−1)+𝛿𝑡 𝒙 ̇_𝑎^(𝑚−1/2)
    ! Update position from Velocity
    do a=1,ntotal
        if((itype(a) .le. itype_real_max) .and. (itype(a) .gt. itype_real_min)) then
            x(:,a)= x(:,a) + dt*vx(:,a)
        endif
    enddo         


    if ((mod(itimestep,save_step).eq.0) .or. (itimestep.eq.1)) call output_time_flow(itimestep,dt)    

end
     

    