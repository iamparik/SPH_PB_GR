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

    if(tempType .eq. 1) etype = 200
        
    if(tempType .eq. 2) etype = 200
    
    if(tempType .eq. 3) etype = 2
        
    if(tempType .eq. 4) etype = 2  


end subroutine