!****************************************************************************
!
!  SUBROUTINE: caseBasedOutput 
!
!****************************************************************************      
subroutine  caseBasedOutput(iterStep,dt)
    use config_parameter, only:dataOutputPath, SPH_dim
    use particle_data ,   only: x, vx,mass, nreal, xstart, delC, p
    
    implicit none
!----------------------------------------------------------------------           
!                   i:      A variable used as a counter in loop
!                   d:      A variable used to iterate dimensions
!                   xname:  Path of the file which needs to be created
!------------------------------------------------------------------------
    integer(4), intent(in) :: iterStep
    real(8), intent(in) :: dt    
    integer(4):: a,d
    character (40) ::x2name
    real(8) :: p_error, TPD,KE, p_true

    if (iterStep .eq. 1) call createOutputFolder(dataOutputPath)

    
    TPD=0.D0
    p_error=0.D0
    KE=0.D0
    do a=1,nreal
        TPD = TPD + norm2(xStart(:,a)-x(:,a))/nreal
        p_true=9.81D0*(0.5D0-x(2,a))* 1000.D0
        p_error= p_error + abs(p_true-p(a))/nreal       
        KE= KE+ 0.5D0*mass(a)*norm2(vx(:,a))**2
    enddo
    
    
!     output result as a continous plot against time   
    if (iterStep .eq. 1) then
        write(x2name, '(A,A)') dataOutputPath,'/energy.dat'
        open (1, file = x2name)
        write(1,'(A)') 'variables = time , TPD , P_error , KE'
        write(1,*) iterStep*dt, TPD, p_error, KE
        close(1)
    else
        write(x2name, '(A,A)') dataOutputPath,'/energy.dat'
        open (1, file = x2name,  position='append')
        write(1,*) iterStep*dt, TPD, p_error, KE
        close(1)
    endif
    


end
      