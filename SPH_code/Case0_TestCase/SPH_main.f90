program SPH3D
 
!-------------------------------------------------------
! Program for running Boundary Integral SPH
! Author- Parikshit B            
!
!   CREATED:        04/13/2022       by  PARIKSHIT BOREGOWDA
!   Last Modified:  04/17/2022       by  PARIKSHIT BOREGOWDA 
!-------------------------------------------------------

use particle_data 
use config_parameter

implicit none

!--------------------------------------
!input_timeStep : An integer Paramter that decides the number of time step iterations 
!                to be run for the simulation
!maztimestep : An integer Paramter that logs the maximum time step iterations
!yesorno :   An integer paramter,that implies no for 0
!            and yes for 1
!ic1,ic2 :   Counters to note clock time 
!crate1, cmax1  :   Itneger kind 8 coutns time in Microseconds
!                    and is used for calling function SYSTEM_CLOCK
!--------------------------------------

real(8) dt
integer(4) input_timeStep,max_timeSteps, yesorno, input_file_type, current_ts
integer(8) :: ic1, crate1, cmax1, ic2
logical ::  runSPH = .true.
integer(4) :: k, a, b , s, Scalar0Matrix1
real(8) :: scalar_factor, Pr_s
real(8), DIMENSION(:), allocatable  :: F_a, F_b, Cdwdx_a, Cdwdx_b, Cdgmas
real(8), DIMENSION(:,:), allocatable :: matrix_factor, delP
real(8), DIMENSION(:), allocatable :: div_vel


   
!start system clock to evealuate time taken to run
! System_clock returns number of seconds from 00:00 CUT on 1 Jan 1970
    call system_clock(count=ic1, count_rate=crate1, count_max=cmax1)

! input particle configuration data
    call inputSPHConfig
    
    allocate(F_a(SPH_dim), F_b(SPH_dim), Cdwdx_a(SPH_dim), Cdwdx_b(SPH_dim), Cdgmas(SPH_dim), &
    & matrix_factor(SPH_dim,SPH_dim))

! Add boundary conditions to the boundary:
    allocate( bdryVal_vel(SPH_dim,maxedge),bdryVal_prs(maxedge), bdryVal_rho(maxedge), bdryVal_temp(maxedge))
    do s =1, etotal
        call BCinputValue(etype(s),bdryVal_vel(:,s), bdryVal_prs(s), bdryVal_rho(s), bdryVal_temp(s))
    enddo

! theoretical maximum for time step.
    call maxTimeStep
    dt = time_step


! Maximum number of timesteps, and ith time step variables are initiatilized as 0 
    max_timeSteps=0
    itimestep=0

!The below requires the user to input on the terminal the maximum time step to run the simulation
    write(*,*)'  ***************************************************'
    write(*,*)'          Please input the additional number of time steps to run simulations'
    write(*,*)'  ***************************************************'
    read(*,*) max_timeSteps

!
    max_timeSteps=max_timeSteps+itimestep
    Allocate(xStart(SPH_dim,ntotal))
    xStart=x
    
    current_ts= itimestep
    
    do itimestep = current_ts+1, max_timesteps
        
        
        call printTimeStep(itimestep,print_step)
        
        call PositionDependentFactors
        
        if (itimestep .eq. 1) call output_initial_config
        
        ! if there was periodic BC applied, for periodic particles and points update variables used 
        if (Allocated(pBC_edges)) then
            call PeriodicParameter(vx(1,:))
            call PeriodicParameter(vx(2,:))
            call PeriodicParameter(rho)
            call PeriodicParameter(p)
        endif
        
        ! Allocate variables necessary
        ! Density updated related paramters
        !allocate(drho(ntotal),rho_prev(ntotal))
        !drho=0.D0
        
        allocate(div_vel(ntotal), delP(SPH_dim, ntotal)) ! this can be reduced by accoutnign for nreal and nedge correctly
        div_vel=0.D0
        delP=0.D0
    
        
        ! Use all particle-particle interaction to find non boundary terms
        do k= 1,niac
            a= pair_i(k)
            b= pair_j(k)
            
            !------------------- Find divergence of velocity (to be used in continuity equation) -------------------------!
            call CorrectionFactorParsing(1,Scalar0Matrix1,scalar_factor,matrix_factor, &
                & gamma_cont(a), gamma_discrt(a), gamma_mat(:,:,a), gamma_mat_inv(:,:,a), xi1_mat_inv(:,:,a), SPH_dim)
            
            Cdwdx_a=dwdx(:,k)
            call CorrectedKernelGradient(Cdwdx_a, scalar_factor, matrix_factor, Scalar0Matrix1, SPH_dim)
            Cdwdx_b=-dwdx(:,k)
            call CorrectedKernelGradient(Cdwdx_b, scalar_factor, matrix_factor, Scalar0Matrix1, SPH_dim)    
            F_a = vx(:,a)
            F_b = vx(:,b)
            call VectorDivergencePtoP(div_vel(a),div_vel(b),F_a,F_b,Cdwdx_a, Cdwdx_b, mass(a), mass(b), rho(a), rho(b), SPH_dim, 1)
            ! -------------------------------------------------------------------------------------------------------------!
            
            !-------------- Find Pressure Gradient term (to be used in momentum equation) --------------!
            
            call CorrectionFactorParsing(1,Scalar0Matrix1,scalar_factor,matrix_factor, &
                & gamma_cont(a), gamma_discrt(a), gamma_mat(:,:,a), gamma_mat_inv(:,:,a), xi1_mat_inv(:,:,a), SPH_dim)
            
            Cdwdx_a=dwdx(:,k)
            call CorrectedKernelGradient(Cdwdx_a, scalar_factor, matrix_factor, Scalar0Matrix1, SPH_dim)
            Cdwdx_b=-dwdx(:,k)
            call CorrectedKernelGradient(Cdwdx_b, scalar_factor, matrix_factor, Scalar0Matrix1, SPH_dim)    
            call ScalarGradientPtoP(delP(:,a),div_vel(b),P(a),P(b),Cdwdx_a, Cdwdx_b, mass(a), mass(b), rho(a), rho(b), SPH_dim, 2)
            !-------------------------------------------------------------------------------------------------------------!
            
            
            
            
        enddo
        
        ! Use all particle-edge interactions to find boundary terms
        do k= 1,eniac
            a= epair_a(k)
            s= epair_s(k)
            b= nedge_rel_edge(s)

            
            !------ Find divergence of velocity for calculating density -------------!
            call CorrectionFactorParsing(1,Scalar0Matrix1,scalar_factor,matrix_factor, &
                & gamma_cont(a), gamma_discrt(a), gamma_mat(:,:,a), gamma_mat_inv(:,:,a), xi1_mat_inv(:,:,a), SPH_dim)
            
            Cdgmas=del_gamma_as(:,k)
            call CorrectedKernelGradient(Cdgmas, scalar_factor, matrix_factor, Scalar0Matrix1, SPH_dim)  
            F_a = vx(:,a)
            F_b = vx(:,b)
            call VectorDivergencePtoB(div_vel(a),vx(:,a),vx(:,b),Cdgmas,SPH_dim, 1)
            ! -----------------------------------------------------------------------!
            
            !------ Find Pressure Gradient term (to be used in momentum equation) -------------!
            call CorrectionFactorParsing(1,Scalar0Matrix1,scalar_factor,matrix_factor, &
                & gamma_cont(a), gamma_discrt(a), gamma_mat(:,:,a), gamma_mat_inv(:,:,a), xi1_mat_inv(:,:,a), SPH_dim)
            
            Cdgmas=del_gamma_as(:,k)
            call CorrectedKernelGradient(Cdgmas, scalar_factor, matrix_factor, Scalar0Matrix1, SPH_dim)  
            Pr_s = P(a) +rho(a)*c_sound*dot_product(vx(:,a)-vx(:,b), surf_norm(:,s))
            call ScalarGradientPtoB(delP(:,a),P(a),Pr_s,Cdgmas,SPH_dim, 2)
            ! -----------------------------------------------------------------------!
            
        enddo
        
        
        ! Update variables for the next time step
        do a =1, ntotal
            if((itype(a) .le. itype_real_max) .and. (itype(a) .gt. itype_real_min)) then
                ! Calcualate density as (𝐷𝜌_𝑎)/𝐷𝑡=− 𝜌_𝑎  ∇∙𝑣_𝑎
                rho(a) = rho(a) - dt* rho(a) * div_vel(a)
                ! Use Hughes density correction if necessary
                if ((HG_density_correction) .and. (rho(a) .le. rho_init)) rho(a)=rho_init
            
            
                ! Update Pressure as it depends on density for WCSPH
                call ParticlePressureEOS(p(a), rho(a), itype(a), itype_virtual) 
                
                
                !Update acceleration terms
                
                
                !Update Velocity
                vx(:,a) = vx(:,a) + dt* (delP(:,a)/rho(a))
            
            
                !Update position
                x(:,a) = x(:,a) + dt* vx(:,a)
            
            endif
        enddo
        
         deallocate(div_vel,delP)
        
        
        !---------------------- free surface detection and PST algorithm -------------------------------------!
        call FreeSurfaceDetection
        if(PSTtype .gt. 0) then
            call ParticleShiftingTechnique
        endif
        
        if( .not. allocated(delC)) call ConcGradient  
        
        
        if ((mod(itimestep,save_step).eq.0) .or. (itimestep.eq.1)) call output_flow_simplified(itimestep,dt)   
        deallocate(delC)
        deallocate(FreeSurfaceVar)
        
       !call EulerIntegration_fluid(itimestep,dt)
        
        if (Allocated(pBC_edges)) call PeriodicBCreset
        
        ! the below needs to be deallocated here because periodic bc can sometimes increase the total number of particles in next loop
        deallocate(w_aa, gamma_discrt,del_gamma_as, del_gamma, xi1_mat, beta_mat,gamma_mat,xi_cont_mat, &
                   & gamma_mat_inv,xi1_mat_inv,xi_cont_mat_inv,b_MLS )
        
        
        ! if (mod(itimestep,backup_step).eq.0)  call backupOutput(itimestep)
        
    enddo
        




!Now net time is calculated
call system_clock(count=ic2)
write (*,*)'        Elapsed time = ', (ic2-ic1)/real(crate1), 'sec'

!All allocated variables need to be deallocated so that the memory is freed

if (Allocated(pBC_edges)) deallocate(pBC_epair_a, pBC_epair_s)

DEALLOCATE(x,vx,mass, rho, p, hsml, itype, mu, temp, xStart) 
Deallocate(surf_norm, edge,nedge_rel_edge, etype )

if(Allocated(dgrho_prev)) DEALLOCATE(dgrho_prev)
if(Allocated(drho)) DEALLOCATE(drho)
if(Allocated(rho_prev)) DEALLOCATE(rho_prev)

!Pause is used so that the terminal is closed only after hitting enter/return on keyboard
pause 
end program