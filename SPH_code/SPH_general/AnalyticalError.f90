!****************************************************************************
!
!  SUBROUTINE:AnalyticalError
!
!  PURPOSE:  Subroutine to store L2 norms of error, Rel Error, analytical solution
!           for each time step and then output a datafile with the same
!
!   CREATED:        12/23/2022       by  PARIKSHIT BOREGOWDA
!   Last Modified:  01/24/2023        by  PARIKSHIT BOREGOWDA 
!****************************************************************************

subroutine AnalyticalError(itimestep,maxtimestep,dt)
use config_parameter, only:dataOutputErrors
use config_geometry, only: BILtype, NumericalSimCase, save_step
use particle_data, only: ErrorL2norm, ntotal, nreal, &
    & vx, temp, mass,rho, gamma_cont, f0max

integer(4), intent(in) :: itimestep,maxtimestep
real(8), intent(in):: dt
integer(4) a
real(8), dimension(:),allocatable::fncn,fncn0,fncnDenom,Error,RelError,Error2,RelError2
real(8), dimension(:),allocatable:: fncnAnalytical, fncnNumerical
real(8), dimension(:,:),allocatable:: L2normfncn
!real(8):: largeNumber
logical:: bdryval

bdryval=.false.
!largeNumber=huge(1.D0)

! allocate Error2Lnorm to include as many  paramters.
!ErrorL2(maxtimestep,3), allows 3 variables to be stored (L2(Error),L2(RelError),L2(fncn0))
if ( .not.allocated(ErrorL2norm) ) then
    allocate(ErrorL2norm(maxtimestep,5))
    ErrorL2norm=0.D0
endif


allocate(fncnAnalytical(ntotal),fncnNumerical(ntotal), fncnDenom(ntotal))
fncnAnalytical=0.D0
fncnNumerical=0.D0
fncnDenom=0.D0

if((NumericalSimCase .eq. 1) .or. (NumericalSimCase .eq. 2) ) then
     call AnalyticalThermalDiffusion(itimestep,dt,fncnAnalytical,bdryval)
     fncnNumerical(1:ntotal)=temp(1:ntotal)
     fncnDenom(1:ntotal)=fncnAnalytical(1:ntotal)
elseif((NumericalSimCase .eq. 3) .or. (NumericalSimCase .eq. 4) ) then
    call AnalyticalFluidFlow(itimestep,dt,fncnAnalytical,bdryval)
    fncnNumerical(1:ntotal)=vx(1,1:ntotal)
    fncnDenom(1:ntotal)=f0max
endif


allocate(fncn(nreal),fncn0(nreal),Error(nreal),RelError(nreal),Error2(nreal),RelError2(nreal))
fncn=0.D0
fncn0=0.D0
Error=0.D0
RelError=0.D0
Error2=0.D0
RelError2=0.D0

 do a = 1, nreal
     fncn(a)=fncnNumerical(a)*dsqrt(mass(a)/rho(a))
     fncn0(a)=fncnAnalytical(a)*dsqrt(mass(a)/rho(a))
     Error(a)=fncn(a)-fncn0(a)
     RelError(a)=Error(a)/fncnDenom(a)
     if(gamma_cont(a) .lt. 1.D0) then
         Error2(a)=fncn(a)-fncn0(a)
         RelError2(a)=Error(a)/fncnDenom(a)
     endif
 enddo
 

! ErrorL2=(∑_i(u_i-u_i^h )^2  V_i )^(1/2)
 ErrorL2norm(itimestep,1)=norm2(Error)
 ErrorL2norm(itimestep,2)=norm2(Error2)
! RelErrorL2=(∑_i((u_i-u_i^h)/u_i )^2 V_i )^(1/2)
 ErrorL2norm(itimestep,3)=dlog10(norm2(RelError))
 ErrorL2norm(itimestep,4)=dlog10(norm2(RelError2))
! trueL2=(∑_i(u_i)^2  V_i )^(1/2)
 ErrorL2norm(itimestep,5)= norm2(fncn0)



if((mod(itimestep,save_step).eq.0)) then    
    call outputErrorAllTime(maxtimestep,itimestep,ErrorL2norm,int(5),dataOutputErrors,BILtype,dt) 
    if(itimestep .eq. maxtimestep) deallocate(ErrorL2norm)
endif


deallocate(fncn,fncn0,Error,RelError,Error2,RelError2,fncnAnalytical,fncnNumerical)


end