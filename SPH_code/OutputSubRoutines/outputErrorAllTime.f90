!****************************************************************************
!
!  SUBROUTINE: OutputErrors
!
!  PURPOSE: This subroutine will output the errors
!
!   CREATED:        12/17/2022       by  PARIKSHIT BOREGOWDA
!   Last Modified:  12/23/2022       by  PARIKSHIT BOREGOWDA 
!****************************************************************************      

subroutine outputErrorAllTime(maxtimestep, currentstep, fx,isize,dataOutputPath,runtype,dt)
    implicit none
!----------------------------------------------------------------------           
!                   i:      A variable used as a counter in loop
!                   d:      A variable used to iterate dimensions
!                   xname:  Path of the file which needs to be created
!------------------------------------------------------------------------
    
    integer(4), intent(in):: maxtimestep,currentstep, isize,runtype
    real(8), intent(in):: fx(maxtimestep,isize)
    real(8), intent(in):: dt
    character(LEN=*) :: dataOutputPath
    character (40) :: xname
    logical :: text_print, NANInfDetect
    integer(4):: i,d,ii    
    real(8):: largeNumber

    largeNumber=huge(1.D0)
    
    
!     output result in tecplot format
    !write(xname, '(A,A,I1,A)') dataOutputPath,'/',runtype,'error.dat'
    
    if(runtype.eq.0) then
            write(xname, '(A,A, I1, A)') dataOutputPath,'/000000', runtype,'error.dat'
    elseif(runtype>0.and.runtype<10) then
            write(xname, '(A,A, I1, A)') dataOutputPath,'/000000', runtype,'error.dat'
    elseif(runtype>9.and.runtype<100) then
            write(xname, '(A,A, I2, A)') dataOutputPath,'/00000', runtype,'error.dat'      
    endif

    
!   Open the file with the path xname, call this as 1
    open (1, file = xname)

!  Write the title on the open file
    write (1, '(A)') 'title="Error across time"'
    
 
!   Specify variables that need to be read by tecplot to correspond to the data fed
    write(1,'(A)') 'variables = time, L2Norm'
    do ii=1, isize
        write(1,'(A,I2,A,I2,A)')'ZONE T="L2norm',runtype,'_',ii,'",F=Point,C=Red'
        do d =1 ,currentstep
            if(isNAN(fx(d,ii))) then
                write(1,1001) dt*d, largeNumber 
                
            elseif( abs(fx(d,ii)) .gt. largeNumber) then
                write(1,1001) dt*d, largeNumber 
                
            else
                write(1,1001) dt*d, fx(d,ii)
            endif
            
        enddo
        
    enddo
    
    
    
        
        
      
!Here is a readable format in which the data is stored
1001  format(12(e22.10,2x))    ! i am using Maximum, it was 9 before
! 2X means skip 2 spaces
      
      close(1)

end
      