!****************************************************************************
!
!  SUBROUTINE: caseBasedOutput 
!
!****************************************************************************      
subroutine  caseBasedOutput(iterStep)
    use config_parameter, only:dataOutputPath, SPH_dim
    use particle_data ,   only: x, vx,mass, nreal, xstart, delC
    
    implicit none
!----------------------------------------------------------------------           
!                   i:      A variable used as a counter in loop
!                   d:      A variable used to iterate dimensions
!                   xname:  Path of the file which needs to be created
!------------------------------------------------------------------------
    integer(4):: a,d,s, iterStep,f_dim,saveStep
    character (40) :: xname, x2name
    real(8) :: delC_avg, TPD,KE
    logical :: text_print


    if (iterStep .eq. 1) call createOutputFolder(dataOutputPath)

    
    TPD=0.D0
    delC_avg=0.D0
    KE=0.D0
    do a=1,nreal
        TPD = TPD + norm2(xStart(:,a)-x(:,a))/nreal
        delC_avg= delC_avg + norm2(delC(:,a))/nreal
        KE= KE+ 0.5D0*mass(a)*norm2(vx(:,a))**2
    enddo
    
    
!     output result as a continous plot against time   
    if (iterStep .eq. 1) then
        write(x2name, '(A,A)') dataOutputPath,'/energy.dat'
        open (1, file = x2name)
        write(1,'(A)') 'variables = Iterations , TPD , delC_avg , KE'
        write(1,*) iterStep, TPD, delC_avg, KE
        close(1)
    else
        write(x2name, '(A,A)') dataOutputPath,'/energy.dat'
        open (1, file = x2name,  position='append')
        write(1,*) iterStep, TPD, delC_avg, KE
        close(1)
    endif
    


end
      