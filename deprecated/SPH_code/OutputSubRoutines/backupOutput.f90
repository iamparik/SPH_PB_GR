!****************************************************************************
!
!  SUBROUTINE: backupOutput 
!
!  PURPOSE: This subroutine will output the input geomtrical data,
!           and can be used to test if the input file is being
!           read correctly
!
!   CREATED:        9/13/2017       by  PARIKSHIT BOREGOWDA
!   Last Modified:  04/18/2022        by  PARIKSHIT BOREGOWDA 
!****************************************************************************      
subroutine  backupOutput(itimestep)
    use config_parameter, only:DataConfigPath, SPH_dim, itype_real_max, itype_real_min
    use particle_data ,   only: itype, etype, x,vx, mass, rho,p, hsml, mu, temp, &
        & nreal, ntotal,etotal, edge
    
    implicit none
!----------------------------------------------------------------------           
!                   i:      A variable used as a counter in loop
!                   d:      A variable used to iterate dimensions
!                   xname:  Path of the file which needs to be created
!------------------------------------------------------------------------
    integer(4), intent(inout) :: itimestep
    integer(4) i,d,s, a,b, dd
    real(8) dgammaMAX
    character (40) :: xname, x2name
    logical :: text_print
    
    

    write(*,*)'  **************************************************'
    write(*,*)'       Write backup of the SPH domain configuration for timestep = ', itimestep
    
    write(xname, '(A,A)') DataConfigPath,'/backup_SimSize.dat'
    open(1,file = xname)
    write(1,1000) nreal, etotal, itimestep
    close(1)

    if(.not. allocated(temp)) then
        allocate(temp(nreal))
        temp=0.D0
    endif
    
    write(*,*)'  **************************************************'
    write(*,*)'       Write backup of the SPH particles...   '
    
    write(xname, '(A,A)') DataConfigPath,'/backup_Particles.dat'
    open(1,file = xname)

! Store all real particle data
    do i=1,nreal
       if((itype(i) .le. itype_real_max) .and. (itype(i) .gt. itype_real_min))  then
           write(1,1001) itype(i), (x(d, i), d=1,SPH_dim), (vx(d, i), d = 1, SPH_dim),mass(i),rho(i), p(i), temp(i) 
       endif
    enddo
    close(1)
    
   
    write(*,*)'  **************************************************'
    write(*,*)'       Write backup of the SPH edges...   '
    
    write(xname, '(A,A)') DataConfigPath,'/backup_Edges.dat'
    open(1,file = xname)
 
! Store all real particle data
    do s=1,etotal
        write(1,1002) etype(s), ((x(d,edge(dd,s)), d=1,SPH_dim), dd=1,SPH_dim) 
    enddo
       
    close(1)


!Here is a readable format in which the data is stored
1000  FORMAT(3(I10,2X))    ! i am using Maximum, it was 9 before
1001  format(I10, 2X, 9(e22.10,2x))    ! i am using Maximum, it was 9 before
1002  format(I10, 2X, 4(e22.10,2x))    ! i am using Maximum, it was 9 before

      

end
      