!****************************************************************************
!
!  SUBROUTINE: ExternalForceAcceleration
!
!  PURPOSE:  This subroutine determines vscous Stress    
!
!   CREATED:        01/13/2022       by  PARIKSHIT BOREGOWDA
!   Last Modified:  01/13/2023       by  PARIKSHIT BOREGOWDA 
!****************************************************************************
subroutine ExternalForceAcceleration(dstress,rho, itype, ntotal)

    use config_parameter, only: SPH_dim, itype_real_max, itype_real_min
    use config_geometry, only: F_ext, g_const
    
    
    implicit none
    
    integer(4), intent(in) :: ntotal
    integer(2), intent(in) :: itype(ntotal)
    real(8), intent(in) :: rho(ntotal)
    real(8) :: dstress(SPH_dim,ntotal)
    integer(4) a, i
    real(8) g
    
    g= -abs(g_const)
    
    
    do a=1,ntotal
        if((itype(a) .le. itype_real_max) .and. (itype(a) .gt. itype_real_min)) then
            do i=1,SPH_dim
                dstress(i,a) = dstress(i,a) + F_ext(i)*rho(a)
            enddo
            dstress(2,a) = dstress(2,a) + g*rho(a)
        endif
    enddo
    
    
    
end
    