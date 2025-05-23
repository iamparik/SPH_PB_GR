program BIPI_SPH
 
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
real(8) ::temp_scalar, driac, w_temp


   
!start system clock to evealuate time taken to run
! System_clock returns number of seconds from 00:00 CUT on 1 Jan 1970
    call system_clock(count=ic1, count_rate=crate1, count_max=cmax1)
    
! Read parameters to configure the SPH simulation from config_data.dat
    call read_config_file    
    
! run packing algorithm according to user input
    call inputExt   
    

    

!Now net time is calculated
call system_clock(count=ic2)
write (*,*)'        Elapsed time = ', (ic2-ic1)/real(crate1), 'sec'

!All allocated variables need to be deallocated so that the memory is freeda

!Pause is used so that the terminal is closed only after hitting enter/return on keyboard
write(*,*) 'The code has finished executing, press return to exit'
pause 
end program