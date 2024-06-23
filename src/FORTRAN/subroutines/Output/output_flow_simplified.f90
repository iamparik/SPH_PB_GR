!****************************************************************************
!
!  SUBROUTINE: output_flow_simplified 
!
!  PURPOSE: This subroutine will output simple flow variables
!
!   CREATED:        09/13/2017       by  PARIKSHIT BOREGOWDA
!   Last Modified:  03/25/2024        by  PARIKSHIT BOREGOWDA 
!****************************************************************************      
subroutine  output_flow_simplified(itimestep,dt)
    use config_parameter, only:dataOutputPath, SPH_dim, itype_real_max, itype_real_min, &
        & etype_real_max, etype_real_min, etype_virtual, itype_virtual, itype_periodic, &
        & save_step, g_const, hydroStaticHeight
    use particle_data ,   only: x, mass, rho, p, vx, &
        & ntotal,etotal, itype, etype, gamma_cont, gamma_discrt,del_gamma, delC, &
        & gamma_density_cont, ve_total, x_ve, vx_ve, mid_pt_for_edge
    
    implicit none
!----------------------------------------------------------------------           
!                   i:      A variable used as a counter in loop
!                   d:      A variable used to iterate dimensions
!                   xname:  Path of the file which needs to be created
!------------------------------------------------------------------------
    integer(4):: i,d,s, itimestep,f_dim
    real(8):: dt
    character (40) :: xname, x2name 
    logical :: text_print

    if (itimestep .eq. 1) call createOutputFolder(dataOutputPath)
    
    !call gamma_density_continuous_leroy
    
!     output result in tecplot format
        if(itimestep.eq.0) then
            write(xname, '(A,A, I1, A)') dataOutputPath,'/t=000000', itimestep,'step_xv.dat'
      else if(itimestep>0.and.itimestep<10) then
            write(xname, '(A,A, I1, A)') dataOutputPath,'/t=000000', itimestep,'step_xv.dat'
      else if(itimestep>9.and.itimestep<100) then
            write(xname, '(A,A, I2, A)') dataOutputPath,'/t=00000', itimestep,'step_xv.dat'
      else if(itimestep>99.and.itimestep<1000) then
            write(xname, '(A,A, I3, A)') dataOutputPath,'/t=0000', itimestep,'step_xv.dat'
      else if(itimestep>999.and.itimestep<10000) then                                   
            write(xname, '(A,A, I4, A)') dataOutputPath,'/t=000', itimestep,'step_xv.dat'
      else if(itimestep>9999.and.itimestep<100000) then
            write(xname, '(A,A, I5, A)') dataOutputPath,'/t=00', itimestep,'step_xv.dat'
      else if(itimestep>99999.and.itimestep<1000000) then
            write(xname, '(A,A, I6, A)') dataOutputPath,'/t=0', itimestep,'step_xv.dat'
      else if(itimestep>999999.and.itimestep<10000000) then
            write(xname, '(A,A, I7, A)') dataOutputPath,'/t=', itimestep,'step_xv.dat'
      endif
    

!   Open the file with the path xname, call this as 1
    open (1, file = xname)

!  Write the title on the open file
    write (1, '(A)') 'title="Starting Config"'

!   Specify variables that need to be read by tecplot to correspond to the data fed
! for 3D output 
    if(SPH_dim.eq.3)    write(1,*) 'variables = x, y, z, mass, rho, vx, vy, vz, p, gamma, alpha, dgammax, dgammay, dgammaz, delcX, delCy, delCz '
! for 2D output
    if(SPH_dim.eq.2)    write(1,*) 'variables = x, y, mass, rho, vx, vy, p, gamma, alpha, dgammax, dgammay, delcX, delCy'
! for 1D output
    if(SPH_dim.eq.1)    write(1,*) 'variables = x, mass, rho, temp, vx, p, gamma, alpha, dgammax, delcX'
  
text_print=.true.
! First export data of all real particles ( not on boundary)
    
    do i=1,ntotal
       if((itype(i) .le. itype_real_max) .and. (itype(i) .gt. itype_real_min))  then
           if(text_print) write(1,'(A)')'ZONE T="Real Particles",F=Point,C=Blue'
           write(1,1001) (x(d, i), d=1,SPH_dim), mass(i),rho(i), (vx(d, i), d = 1, SPH_dim), p(i), gamma_cont(i), gamma_discrt(i), (del_gamma(d, i), d = 1, SPH_dim), (delC(d, i), d = 1, SPH_dim)
           text_print=.false. 
       endif
       
       if(itype(i) .eq. itype_real_min) then
           rho(i) =0.D0
           vx(:,i)=0.D0
           p(i)=0.D0
       endif
       
    enddo


    text_print=.true.
! Export data of all particles representative of boundary, which are used in interpolation,
! called edge particles
    
    do s=1,etotal 
        if(text_print) write(1,'(A)')'ZONE T="Edge Mid-Points",F=Point,C=Red'
        write(1,1001) (mid_pt_for_edge(d,s), d=1,SPH_dim), 0.D0 ,0.D0, (0.D0, d = 1, SPH_dim), 0.D0, 0.D0, 0.D0, (0.D0, d = 1, SPH_dim), (0.D0, d = 1, SPH_dim) 
        text_print=.false.                           
    enddo

    text_print=.true.       
    
! Export all vertices of edges.
    
    do i=1, ve_total
        if(text_print)  write(1,'(A)')'ZONE T="Edge vertices",F=Point,C=Black'
        write(1,1001) (x_ve(d, i), d=1,SPH_dim),0,0, (vx_ve(d, i), d = 1, SPH_dim), 0, 0, 0, (0.D0, d = 1, SPH_dim), (0, d = 1, SPH_dim)    
        text_print=.false.
    enddo
   
    text_print=.true.
! Export data of all copied real particles due to periodic bc (not on boundary)
   
    do i=1,ntotal
       if( ((itype(i)-itype_periodic) .le. itype_real_max) .and. ((itype(i)-itype_periodic) .ge. itype_real_min)  )  then
           if(text_print)  write(1,'(A)')'ZONE T="Periodic Real Particles",F=Point,C=Purple'
           write(1,1001) (x(d, i), d=1,SPH_dim), mass(i),rho(i), (vx(d, i), d = 1, SPH_dim), p(i), gamma_cont(i), gamma_discrt(i), (del_gamma(d, i), d = 1, SPH_dim), (0, d = 1, SPH_dim) 
           text_print=.false.
       endif
    enddo


!Here is a readable format in which the data is stored
1001 format(14(e22.10,2x))    ! i am using Maximum, it was 9 before


1002  FORMAT(8(2X, e22.10))    ! i am using Maximum, it was 9 before
      
      close(1)

end
      