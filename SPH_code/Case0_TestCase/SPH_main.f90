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
integer(8) :: ic1, crate1, cmax1, ic2,s 
logical ::  runSPH = .true.


   
!start system clock to evealuate time taken to run
! System_clock returns number of seconds from 00:00 CUT on 1 Jan 1970
    call system_clock(count=ic1, count_rate=crate1, count_max=cmax1)

! input particle configuration data
    call inputSPHConfig

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
        
        !------------------------DENSITY UPDATE -----------------------------!
        allocate(drho(ntotal),rho_prev(ntotal))
        drho=0.D0
        
        call EulerIntegration_fluid(itimestep,dt)
        
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