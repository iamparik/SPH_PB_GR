!  Case2PeriodicCylinder.f90 
!
!  FUNCTIONS:
!  Case2PeriodicCylinder - Entry point of console application.
!

!****************************************************************************
!
!  PROGRAM: Case2PeriodicCylinder
!
!  PURPOSE:  Entry point for the console application.
!
!****************************************************************************

    program Case2PeriodicCylinder
    
    use particle_data ,   only: rho,p,mu, x,vx,mass, temp, &
        & hsml,itype,surf_norm,edge,nedge_rel_edge,itimestep, etype, &
        & dgrho_prev, drho , rho_prev, xStart, ntotal
    use config_geometry

    implicit none
    
    integer(4) maxtimestep, current_ts
    real(8) dt

    ! Read all variables from the config Files for the given case geometry
    call read_config_file
    
    ! Read the minimum step required to run the above simulation
    call WCSPHminTimeStep(dt)
    
    ! Maximum number of timesteps, and ith time step variables are initiatilized as 0 
    maxtimestep=0
    
    ! Input particle and boundary condiguration
    if(ExtInputMeshType .gt. 0) then
        call inputExt(maxTimeStep)
        current_ts=0
    endif
    
    
    do itimestep = current_ts+1, current_ts+maxtimestep
        
        !   [itimestep mod print step] like 6mod4=2 is measured to check if
        !   itimestep is a multiple of printstep, if yes then value of modulo
        !   will be zero and the below if loop is executed    
        if (mod(itimestep,print_step).eq.0) then 
            write(*,*)'______________________________________________'
            write(*,*)'  current    number     of     time     step =',     itimestep
            write(*,*)'______________________________________________'
        endif 
        
        ! Subroutine to recaclulate kernels, kernel derivative and 
        ! all correction factors that dpend on position of particle
        call PositionDependentFactors

            
        ! if there was periodic BC applied, for periodic particles and points update variables used 
        if (Allocated(pBC_edges)) then
            call PeriodicParameter(vx(1,:))
            call PeriodicParameter(vx(2,:))
            call PeriodicParameter(rho)
            call PeriodicParameter(p)
        endif
        
        !------------------------DENSITY UPDATE -----------------------------!
        
        allocate(drho(ntotal),rho_prev(ntotal))
            drho=0.D0            
            rho_prev(1:ntotal)=rho(1:ntotal)
            
            !Update Density of particles using Continuity Equation
            call continuityDensityOperator(drho,vx,rho_prev, ConDivtype)
   
            if (itimestep .eq. 1) then
                do a=1,ntotal
                    if((itype(a) .le. itype_real_max) .and. (itype(a) .gt. itype_real_min)) then
                        rho(a)=rho_prev(a)
                    endif      
                enddo
            else
                !call FreeSurfaceDetection                
                ! The below is used for continuity density updates
                do a=1,ntotal
                    !drho(a)= - drho(a)*rho_prev(a)
                    if((itype(a) .le. itype_real_max) .and. (itype(a) .gt. itype_real_min)) then                        
                        rho(a)=rho_prev(a) + dt*drho(a)
                        
                        !Below to have all free surf particles except those next to wall have 0 pressure boundary
                        !if((FreeSurfaceVar(a) .le. 1.5D0) .and. (gamma_cont(a) .eq. 1.D0)) rho(a) =rho_init
                        
                        !below is supposidly hughes and grammes condition
                        if (HG_density_correction) then
                            if(rho(a) .le. rho_init) rho(a)=rho_init
                        endif
                    endif      
                enddo                
              !deallocate(FreeSurfaceVar)  
            endif
            
           deallocate(drho,rho_prev) 
            
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
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
        
        ! Periodic boundary resets particles
        if (Allocated(pBC_edges)) then
            call PeriodicBCreset
        endif
        
        ! the below needs to be deallocated here because periodic bc
        ! can sometimes increase the total number of particles in next loop
        deallocate(w_aa, gamma_discrt,del_gamma_as, del_gamma, xi1_mat, beta_mat,gamma_mat, &
                   & gamma_mat_inv,xi1_mat_inv,b_MLS )
        
    enddo
      
    deallocate(gamma_cont, epair_a, epair_s, pair_i, pair_j, w,dwdx)    
    
    if (Allocated(pBC_edges)) then
        deallocate(pBC_epair_a, pBC_epair_s)
    endif
        
        
    end program Case2PeriodicCylinder

