!****************************************************************************
!
!  SUBROUTINE: OutputErrors
!
!  PURPOSE: This subroutine will output the errors
!
!   CREATED:        12/17/2022       by  PARIKSHIT BOREGOWDA
!   Last Modified:  02/06/2023       by  PARIKSHIT BOREGOWDA 
!****************************************************************************      
subroutine  outputError(n,RE,REB,REI,s,dataOutputPath,nreal, hsml_factor_name,hfc)
    implicit none
!----------------------------------------------------------------------           
!                   i:      A variable used as a counter in loop
!                   d:      A variable used to iterate dimensions
!                   xname:  Path of the file which needs to be created
!------------------------------------------------------------------------
    
    integer(4), intent(in):: n, hfc
    real(8), intent(in), dimension(2,n):: RE, REB, REI
    real(8), intent(in):: s
    character(9), intent(in) :: hsml_factor_name
    character(LEN=*) :: dataOutputPath
    character (80) :: xname
    logical :: text_print
    integer(4):: nreal,i,d
    
!     output result in tecplot format
      if(nreal.eq.0) then
            write(xname, '(A,A,I1, A, I1, A, A)') dataOutputPath,'/',hfc,'_nreal=000000', nreal,hsml_factor_name,'error.dat'
      else if(nreal>0.and.nreal<10) then
            write(xname, '(A,A,I1,A, I1, A, A)') dataOutputPath,'/',hfc,'_nreal=000000', nreal,hsml_factor_name,'error.dat'
      else if(nreal>9.and.nreal<100) then
            write(xname, '(A,A,I1,A, I2, A, A)') dataOutputPath,'/',hfc,'_nreal=00000', nreal,hsml_factor_name,'error.dat'
      else if(nreal>99.and.nreal<1000) then
            write(xname, '(A,A,I1, A,I3,A, A)') dataOutputPath,'/',hfc,'_nreal=0000', nreal,hsml_factor_name,'error.dat'
      else if(nreal>999.and.nreal<10000) then                                   
            write(xname, '(A,A,I1, A,I4,A, A)') dataOutputPath,'/',hfc,'_nreal=000', nreal,hsml_factor_name,'error.dat'
      else if(nreal>9999.and.nreal<100000) then
            write(xname, '(A,A,I1, A,I5,A, A)') dataOutputPath,'/',hfc,'_nreal=00', nreal,hsml_factor_name,'error.dat'
      else if(nreal>99999.and.nreal<1000000) then
            write(xname, '(A,A,I1, A,I6,A, A)') dataOutputPath,'/',hfc,'_nreal=0', nreal,hsml_factor_name,'error.dat'
      else if(nreal>999999.and.nreal<10000000) then
            write(xname, '(A,A,I1, A,I7,A, A)') dataOutputPath,'/',hfc,'_nreal=', nreal,hsml_factor_name,'error.dat'
      endif

!   Open the file with the path xname, call this as 1
    open (1, file = xname)

!  Write the title on the open file
    write (1, '(A)') 'title="Error and Relative"'
    
 
!   Specify variables that need to be read by tecplot to correspond to the data fed
    write(1,'(A)') 'variables = s, L2Error, LInfError, L2ErrorBdry, LInfErrorBdry, L2ErrorIns, LInfErrorIns'

    
    do i=1,n
        write(1,*)'ZONE T="Op',i,'",F=Point,C=Red'
        write(1,1001) s, (RE(d, i), d=1,2),(REB(d, i), d=1,2), (REI(d, i), d=1,2)
         
    enddo
    
    !do i=1,n
    !    write(1,1001) s, (RE(d, i), d=1,2),(REB(d, i), d=1,2), (REI(d, i), d=1,2)
    !enddo
    

      
!Here is a readable format in which the data is stored
1001  format(12(e22.10,2x))    ! i am using Maximum, it was 9 before
! 2X means skip 2 spaces
      
      close(1)

end
      