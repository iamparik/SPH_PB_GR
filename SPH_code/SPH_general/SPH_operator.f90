!****************************************************************************
!
!  SUBROUTINE: SPH_operator
!
!  PURPOSE:  Subroutine to test different SPH operations
!
!   CREATED:        04/13/2022       by  PARIKSHIT BOREGOWDA
!   Last Modified:  07/11/2023         by  PARIKSHIT BOREGOWDA 
!****************************************************************************

subroutine SPH_operator(maxtimestep)

use config_geometry, only: timeIntegrationScheme             


implicit none

integer(4) maxtimestep, current_ts, yesorno 
real(8) dt

!Determine if the problem is time indepedent scheme (like steady state or for testing SPHOperators)

if(timeIntegrationScheme .eq. 0) then
    call TimeIndependentSPHOperations
elseif(timeIntegrationScheme .eq. 1) then
    call EulerSchemeSPHOperations(maxtimestep)
elseif(timeIntegrationScheme .eq. 2) then
    call LeapFrogKDKSchemeSPHOperations(maxtimestep)
endif



end