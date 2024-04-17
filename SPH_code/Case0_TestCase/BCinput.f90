!****************************************************************************
!
!  SUBROUTINE: BCinput
!
!  PURPOSE:  Subroutine that allows imposing different types of BC for 
!               for different boundary groups
!
!****************************************************************************
subroutine BCinput(etype, tempType)

integer etype, tempType

    if(tempType .eq. 1) etype = 3
        
    if(tempType .eq. 2) etype = 3
    
    if(tempType .eq. 3) etype = 3
        
    if(tempType .eq. 4) etype = 3  


end subroutine