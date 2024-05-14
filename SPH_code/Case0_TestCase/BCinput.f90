!****************************************************************************
!
!  SUBROUTINE: BCinput
!
!  PURPOSE:  Subroutine that allows imposing different types of BC  
!               for different boundary groups
!
!****************************************************************************
subroutine BCinput(etype, tempType)
    implicit none
    integer etype, tempType

    if(tempType .eq. 1) etype = 200
        
    if(tempType .eq. 2) etype = 200
    
    if(tempType .eq. 3) etype = 2
        
    if(tempType .eq. 4) etype = 4  


end subroutine

    
    
subroutine BCinputValue(etype,bdryVal_var,num_bdry_var, current_time)
    use config_parameter , only : rho_init
    implicit none
    integer(4), intent(in) :: etype, num_bdry_var
    real(8), dimension(num_bdry_var), intent(out) :: bdryVal_var
    real(8), intent(in) :: current_time
    
    if(etype .eq. 2) then
        
        bdryVal_var(1)= 0.D0
        bdryVal_var(2)= 0.D0
        bdryVal_var(3)= rho_init 
    endif
    
    if(etype .eq. 4) then
        
        bdryVal_var(1)= 0.D0
        bdryVal_var(2)= 1.D-2*1.D2*cos(current_time/(1.D-2))
        bdryVal_var(3)= rho_init
    endif
    
end subroutine



    
    