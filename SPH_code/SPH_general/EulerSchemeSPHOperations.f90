!****************************************************************************
!
!  SUBROUTINE: EulerSchemeSPHOperations
!
!  PURPOSE:  Subroutine to test different SPH operations for Euler Integration Scheme
!
!   CREATED:        04/13/2022       by  PARIKSHIT BOREGOWDA
!   Last Modified:  07/11/2023        by  PARIKSHIT BOREGOWDA 
!****************************************************************************

subroutine EulerSchemeSPHOperations(maxtimestep)

use config_parameter, only: SPH_dim, itype_real_max, itype_real_min, &
        & time_step ,NumericalSimCase, print_step, save_step,backup_step         
use particle_data, only:itimestep, w_aa, w, dwdx, &
        & gamma_discrt, gamma_cont, del_gamma_as, del_gamma, &
        & xi1_mat, beta_mat,gamma_mat,xi_cont_mat, &
        & gamma_mat_inv,xi1_mat_inv,xi_cont_mat_inv, &
        & pBC_edges, pBC_epair_a, pBC_epair_s, &
        & epair_a,epair_s,  pair_i, pair_j, &
        & dynamicProblem, b_MLS, &
        & max_vel, KE,PE, TE, delC, delCAvg, delCMax, vx, ntotal, itype, mass

implicit none

integer(4) maxtimestep, current_ts, yesorno,a,b
real(8) dt


! The configaration of the input geometry and features are verified by
! printing the input data to be graphically displayed

! The below can be updated to take input time step
yesorno=0
if (yesorno .eq. 0) dt=time_step

! Output timestep with parameters at 0th timestep
!call output_time_thermal(itimestep, dt)

! store the current timestep seperately, before starting the time integration scheme
current_ts= itimestep



do itimestep= current_ts+1, current_ts+maxtimestep
    

!   [itimestep mod print step] like 6mod4=2 is measured to check if
!   itimestep is a multiple of printstep, if yes then value of modulo
!   will be zero and the below if loop is executed    
    if (mod(itimestep,print_step).eq.0) then 
        write(*,*)'______________________________________________'
        write(*,*)'  current    number     of     time     step =',     itimestep
        write(*,*)'______________________________________________'
    endif    
    
    
  
    ! Ensure neighbor search and correction factors are caclulated only once for static problem
    if( (itimestep .eq.  current_ts+1) .or. dynamicProblem) then
    
    ! Subroutine to recaclulate kernels, kernel derivative and 
    ! all correction factors that dpend on position of particle
        call PositionDependentFactors
    
    ! Account for periodic boundary, by settign up periodic bc and copying particles   
        if (itimestep .eq. 1) then
            call output_initial_config
            !call  output_time_thermal(0, dt)
        endif
    
    endif
    
!Call a time integration scehme/ operator test function to make various calulations
    if((NumericalSimCase .eq. 1) .or. (NumericalSimCase .eq. 2) ) then
        call EulerIntegration_thermal(itimestep,dt)
        call AnalyticalError(itimestep,maxtimestep,dt)
    elseif(NumericalSimCase .eq. 3 ) then
        call EulerIntegration_fluid(itimestep,dt)
        call AnalyticalError(itimestep,maxtimestep,dt)
    else
        call EulerIntegration_fluid(itimestep,dt)
    endif

    
    ! Periodic boundary resets particles only for dynamic flows
    if( dynamicProblem) then   
        if (Allocated(pBC_edges)) then
            call PeriodicBCreset
        endif
   
        ! the below needs to be deallocated here because periodic bc can sometimes increase the total number of particles in next loop
        deallocate(w_aa, gamma_discrt,del_gamma_as, del_gamma, xi1_mat, beta_mat,gamma_mat,xi_cont_mat, &
                   & gamma_mat_inv,xi1_mat_inv,xi_cont_mat_inv,b_MLS )
    endif

        !gamma_cont and gamma_cont_prev was initialized with maxn

    if (mod(itimestep,backup_step).eq.0)  call backupOutput(itimestep)
enddo

deallocate(gamma_cont, epair_a, epair_s, pair_i, pair_j, w,dwdx)

!deallocate all correction factors for static problems, reset periodic BC
if( .not. dynamicProblem) then
    if (Allocated(pBC_edges)) then
            call PeriodicBCreset
    endif
    
    deallocate(w_aa, gamma_discrt,del_gamma_as, del_gamma, xi1_mat, beta_mat,gamma_mat,xi_cont_mat, &
                   & gamma_mat_inv,xi1_mat_inv,xi_cont_mat_inv,b_MLS )
endif
     

if (Allocated(pBC_edges)) then
    deallocate(pBC_epair_a, pBC_epair_s)
endif

! epairs, pairs, w, dw/dx, correction factors, all need to be deallocated.


end