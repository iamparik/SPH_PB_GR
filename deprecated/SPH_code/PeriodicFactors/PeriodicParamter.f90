!****************************************************************************
!
!  SUBROUTINE: PeriodicParamters
!
!  PURPOSE: For simulations having periodic boundary condition
!           the duplicate particles from periodicity are updated with
!           their original particle's parameters
!
!   CREATED:        09/24/2022       by  PARIKSHIT BOREGOWDA
!   Last Modified:  09/24/2022       by  PARIKSHIT BOREGOWDA 
!****************************************************************************
    
    
subroutine PeriodicParameter(fx)
    use particle_data, only: nperiodic,pBC_duplicate_pair, ntotal
    implicit none
    integer(4) a
    real(8) ::fx(ntotal)
    
    
    do  a= 1,nperiodic
        fx(pBC_duplicate_pair(2,a))=fx(pBC_duplicate_pair(1,a))
    enddo
    
end

