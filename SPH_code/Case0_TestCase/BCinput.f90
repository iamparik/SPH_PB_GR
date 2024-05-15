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

    
    
subroutine BCinputValue(bdryVal_seg,num_bdry_seg, bdryVal_ve, num_bdry_ve, etype, dim, current_time)
    use config_parameter , only : rho_init
    implicit none
    integer(4), intent(in) :: etype, num_bdry_seg, num_bdry_ve, dim
    real(8), dimension(num_bdry_seg), intent(inout) :: bdryVal_seg
    real(8), dimension(num_bdry_ve,dim), intent(inout) :: bdryVal_ve
    real(8), intent(in) :: current_time
    integer(4) :: d
    
    
    if(etype .eq. 2) then
        
        bdryVal_seg(1)= 0.D0
        bdryVal_seg(2)= 0.D0
        bdryVal_seg(3)= rho_init 
    endif
    
    if(etype .eq. 4) then
        
        do d=1,dim
            bdryVal_ve(4,d)= 1.D-2*1.D2*cos(current_time/(1.D-2))
        enddo
        
        
        bdryVal_seg(1)= 0.D0
        bdryVal_seg(2)= (bdryVal_ve(4,1)+bdryVal_ve(4,2))/2.D0
        !bdryVal_seg(2)= 1.D-2*1.D2*cos(current_time/(1.D-2))
        bdryVal_seg(3)= rho_init
    endif
    
end subroutine



    
    