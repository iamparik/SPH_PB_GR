!****************************************************************************
!
!  SUBROUTINE: output_time_thermal 
!
!  PURPOSE: This subroutine will output the input geomtrical data,
!           and can be used to test if the input file is being
!           read correctly for only thermal diffusion
!
!   CREATED:        09/13/2017       by  PARIKSHIT BOREGOWDA
!   Last Modified:  08/25/2022        by  PARIKSHIT BOREGOWDA 
!****************************************************************************      
subroutine  output_time_thermal(itimestep,dt)
    use config_parameter, only:dataOutputPath, SPH_dim, itype_real_max, itype_real_min, &
        & etype_real_max, etype_real_min, etype_virtual, itype_virtual, itype_periodic
    use particle_data ,   only: x, mass, rho, temp, nedge_rel_edge, &
        & ntotal,etotal, itype, etype
    
    implicit none
!----------------------------------------------------------------------           
!                   i:      A variable used as a counter in loop
!                   d:      A variable used to iterate dimensions
!                   xname:  Path of the file which needs to be created
!------------------------------------------------------------------------
    integer(4):: i,d,s, itimestep,f_dim
    real(8):: dt
    character (40) :: xname
    logical :: text_print

!     output result in tecplot format

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
    if(SPH_dim.eq.3)    write(1,*) 'variables = x, y, z, mass, rho, temp '
! for 2D output
    if(SPH_dim.eq.2)    write(1,*) 'variables = x, y, mass, rho, temp'
! for 1D output
    if(SPH_dim.eq.1)    write(1,*) 'variables = x, mass, rho, temp'
  
text_print=.true.
! First export data of all real particles ( not on boundary)
    
    do i=1,ntotal
       if((itype(i) .le. itype_real_max) .and. (itype(i) .gt. itype_real_min))  then
           if(text_print) write(1,'(A)')'ZONE T="Real Particles",F=Point,C=Blue'
           write(1,1001) (x(d, i), d=1,SPH_dim), mass(i),rho(i), temp(i) 
           text_print=.false. 
       endif
       
       if(itype(i) .eq. itype_real_min) then
           temp(i)=0.D0
       endif
    enddo
    
text_print=.true.
! Export data of all marker points used to capture data
    
    do i=1,ntotal
       if(itype(i) .eq. itype_real_min) then
            if(text_print) then
                !change f_dim for different size arrays
                if(itimestep.ne.0) call marker_point_value(temp)
            endif
            
            if(text_print) write(1,'(A)')'ZONE T="Marker Points",F=Point,C=Cyan'
            write(1,1001) (x(d, i), d=1,SPH_dim), mass(i),rho(i), temp(i) 
            text_print=.false. 
       endif
    enddo
    
    

    text_print=.true.
! Export data of all particles representative of boundary, which are used in interpolation,
! called edge particles
    
    do s=1,etotal 
        if((etype(s) .le. etype_real_max) .and. (etype(s) .ge. etype_real_min)) then
            if(text_print) write(1,'(A)')'ZONE T="Edge Particles",F=Point,C=Red'
	        i=nedge_rel_edge(s)  ! this needs to be changed for edges with mroe than one poitn representation. This will then be itereated 
            write(1,1001) (x(d, i), d=1,SPH_dim), mass(i),rho(i), temp(i)
            text_print=.false.
        endif                                       
    enddo

    text_print=.true.       
    
! Export all Ghost particles/ponts.
    
    do i=1, ntotal
        if(itype(i) .eq. 0) then
            if(text_print)  write(1,'(A)')'ZONE T="Edge vertices",F=Point,C=Black'
            write(1,1001) (x(d, i), d=1,SPH_dim),0,0,0     
            text_print=.false.
        endif
    enddo
   
    text_print=.true.
! Export data of all copied real particles due to periodic bc (not on boundary)
   
    do i=1,ntotal
       if( ((itype(i)-itype_periodic) .le. itype_real_max) .and. ((itype(i)-itype_periodic) .ge. itype_real_min)  )  then
           if(text_print)  write(1,'(A)')'ZONE T="Periodic Real Particles",F=Point,C=Purple'
           write(1,1001) (x(d, i), d=1,SPH_dim), mass(i),rho(i), temp(i) 
           text_print=.false.
       endif
    enddo

    text_print=.true.
! Export data of all particles representative of boundary, which are used in interpolation,
! called edge particles
    
    do s=1,etotal 
        if( .not.((etype(s) .le. etype_real_max) .and. (etype(s) .ge. etype_real_min))) then
            if(text_print)   write(1,'(A)')'ZONE T="Virtual Edge Particles",F=Point,C=Red'
	        i=nedge_rel_edge(s)  ! this needs to be changed for edges with more than one point representation. This will then be iterated 
            write(1,1001) (x(d, i), d=1,SPH_dim), mass(i),rho(i), temp(i)
            text_print=.false.
        endif
    enddo

    text_print=.true.
! Export all Ghost particles/ponts.
    
    do i=1, ntotal
        if((itype(i) .gt. itype_virtual) .and. (mod(itype(i),itype_virtual) .eq. 0) ) then
            if(text_print)   write(1,'(A)')'ZONE T="Virtual reference points",F=Point,C=Black'
            write(1,1001) (x(d, i), d=1,SPH_dim), 0,0,0   
            text_print=.false.
        endif
    enddo
       

!Here is a readable format in which the data is stored
1001  format(9(e22.10,2x))    ! i am using Maximum, it was 9 before

      
      close(1)

end
      