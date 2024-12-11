!****************************************************************************
!
!  SUBROUTINE: caseBasedOutput 
!
!****************************************************************************      
subroutine  caseBasedOutput(iterStep,dt)
    use config_parameter, only:dataOutputPath, SPH_dim, save_step
    use particle_data ,   only: x, vx,p, mass, nreal, xstart, delC
    
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
    real(8) :: delC_avg, TPD,KE

    if (iterStep .eq. 1) call createOutputFolder(dataOutputPath)
    
    if ((mod(iterStep,save_step).eq.0) .or. (iterStep.eq.1)) then

        !     output result as a continous plot against time   
        !write(x2name, '(A,A,I,A)') dataOutputPath,'/p',iterStep,'.dat'
        write(x2name, '(A,A, "p", I0, ".dat")') dataOutputPath,'/', iterStep
        open (1, file = x2name)
        write(1,'(A)') 'variables = x,y,vx,vy,p'
        do a=1,nreal
            write(1,'(e22.10, A, e22.10, A, e22.10, A, e22.10, A, e22.10)') &
                        & x(1,a), ',', x(2,a), ',', vx(1,a), ',', vx(2,a), ',', p(a)
        enddo
        close(1)
    endif
    

    


end
      