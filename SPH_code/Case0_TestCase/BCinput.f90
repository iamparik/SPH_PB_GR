!****************************************************************************
!
!  SUBROUTINE: BCinput
!
!  PURPOSE:  Subroutine that allows imposing different types of BC for 
!               for different boundary groups
!
!****************************************************************************
subroutine BCinput(etype, tempType)
    implicit none
    integer etype, tempType

    if(tempType .eq. 1) etype = 200
        
    if(tempType .eq. 2) etype = 200
    
    if(tempType .eq. 3) etype = 2
        
    if(tempType .eq. 4) etype = 2  


end subroutine

    
    
subroutine BCinputValue(etype,bdryVal_vel, bdryVal_prs, bdryVal_rho, bdryVal_temp)
    use config_parameter , only : rho_init
    implicit none
    integer(4), intent(in) :: etype
    real(8), dimension(2), intent(out) :: bdryVal_vel
    real(8), intent(out) :: bdryVal_prs
    real(8), intent(out) :: bdryVal_rho
    real(8), intent(out) :: bdryVal_temp
    
    if(etype .eq. 2) then
        
        bdryVal_vel(1)= 0.D0
        bdryVal_vel(2)= 0.D0

        bdryVal_prs  = 0.D0 
        bdryVal_rho  = rho_init 
        bdryVal_temp =300 
    endif
    
end subroutine



    
    