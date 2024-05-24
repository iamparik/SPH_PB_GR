!****************************************************************************
!
!  SUBROUTINE: Output_initial_config 
!
!  PURPOSE: This subroutine will output the input geomtrical data,
!           and can be used to test if the input file is being
!           read correctly
!
!   CREATED:        9/13/2017       by  PARIKSHIT BOREGOWDA
!   Last Modified:  04/18/2022        by  PARIKSHIT BOREGOWDA 
!****************************************************************************      
subroutine  output_initial_config
    use config_parameter, only:dataOutputPath, DataConfigPath, SPH_dim, itype_real_max, itype_real_min, &
        & etype_real_max, etype_real_min, etype_virtual, etype_periodic, itype_virtual, itype_periodic
    use particle_data ,   only: x,vx, mass, rho, mu, temp, surf_norm, nedge_rel_edge, &
        & ntotal,etotal, itype, etype, max_vel, KE,PE, TE, delC, delCAvg, delCMax,dynamicProblem, delCL2, &
        & ga_Max,ga_Avg,ga_L2, g_a_min, g_a_max
    
    implicit none
!----------------------------------------------------------------------           
!                   i:      A variable used as a counter in loop
!                   d:      A variable used to iterate dimensions
!                   xname:  Path of the file which needs to be created
!------------------------------------------------------------------------
    integer(4) i,d,s, a,b
    real(8) dgammaMAX
    character (40) :: xname, x2name
    logical :: text_print
    
    if(.not. allocated(temp)) then
        allocate(temp(ntotal))
        temp=0.D0
    endif
    
    call createOutputFolder(dataOutputPath)
    !
    !if( dynamicProblem) then 
    !     call ConcGradient
    !     call evolvingVariables(1,1)
    !
    ! !     output result in tecplot format
    !        write(x2name, '(A,A)') dataOutputPath,'/energy.dat'
    !        open (1, file = x2name)
    !        write(1,*) 'variables = time, KE, PE, TE, u_max , delcAvg, delCMax, delCL2'        
    !        write(1,1002) 0.D0, KE(1), PE(1), TE(1), max_vel(1), delcAvg(1), delCMax(1), delCL2(1)
    !        close(1)   
    !        
    !    deallocate(delC,KE,PE,TE,max_vel, delcAvg, delCMax, delCL2, ga_Max,ga_Avg,ga_L2, g_a_min, g_a_max)
    !
    !endif
    !
    !
    

!     output result in tecplot format

!   Write the path of the file (with file name) to be edited
    write(xname, '(A,A, I7, A)') DataConfigPath,'/init_geom.dat'

!   Open the file with the path xname, call this as 1
    open (1, file = xname)

!  Write the title on the open file
    write (1, '(A)') 'title="Starting Config"'

!   Specify variables that need to be read by tecplot to correspond to the data fed
! for 3D output 
    if(SPH_dim.eq.3)    write(1,*) 'variables = x, y, z, vx, vy, vz, mass, rho, mu, Temp, surf_normx, surf_normy, surf_normz '
! for 2D output
    if(SPH_dim.eq.2)    write(1,*) 'variables = x, y, vx, vy, mass, rho, mu, Temp, surf_normx, surf_normy'
! for 1D output
    if(SPH_dim.eq.1)    write(1,*) 'variables = x, vx, mass, rho, mu, Temp, surf_normx'
  
text_print=.true.
! First export data of all real particles ( not on boundary)
    
    do i=1,ntotal
       if((itype(i) .le. itype_real_max) .and. (itype(i) .gt. itype_real_min))  then
           if(text_print) write(1,'(A)')'ZONE T="Real Particles",F=Point,C=Blue'
           write(1,1001) (x(d, i), d=1,SPH_dim), (vx(d, i), d = 1, SPH_dim),mass(i),rho(i), mu(i),temp(i),(0, d=1,SPH_dim) 
           text_print=.false. 
       endif
    enddo
    
text_print=.true.
    ! Export data of all marker points used to capture data
    do i=1,ntotal
       if(itype(i) .eq. itype_real_min)  then
           if(text_print) write(1,'(A)')'ZONE T="Marker points",F=Point,C=Cyan'
           write(1,1001) (x(d, i), d=1,SPH_dim), (vx(d, i), d = 1, SPH_dim),mass(i),rho(i), mu(i),temp(i),(0, d=1,SPH_dim) 
           text_print=.false. 
       endif
    enddo

    text_print=.true.
! Export data of all particles representative of boundary, which are used in interpolation,
! called edge particles
    
    do s=1,etotal 
        if(((etype(s) .le. etype_real_max) .and. (etype(s) .ge. etype_real_min)) .or. (etype(s) .eq. etype_periodic)) then
            if(text_print) write(1,'(A)')'ZONE T="Edge Particles",F=Point,C=Red'
	        i=nedge_rel_edge(s)  ! this needs to be changed for edges with mroe than one poitn representation. This will then be itereated 
            write(1,1001) (x(d, i), d=1,SPH_dim), (vx(d, i), d = 1, SPH_dim),mass(i),rho(i), mu(i), temp(i), (surf_norm(d,s), d=1,SPH_dim)
            text_print=.false.
        endif                                       
    enddo
    
    

    text_print=.true.
! Export all Ghost particles/ponts.
    
    do i=1, ntotal
        if(itype(i) .eq. 0) then
            if(text_print)  write(1,'(A)')'ZONE T="Edge vertices",F=Point,C=Black'
            write(1,1001) (x(d, i), d=1,SPH_dim), (0, d = 1, SPH_dim),0,0,0,0,(0, d=1,SPH_dim)     
            text_print=.false.
        endif
    enddo
   
    text_print=.true.
! Export data of all copied real particles due to periodic bc (not on boundary)
   
    do i=1,ntotal
       if( ((itype(i)-itype_periodic) .le. itype_real_max) .and. ((itype(i)-itype_periodic) .ge. itype_real_min)  )  then
           if(text_print)  write(1,'(A)')'ZONE T="Periodic Real Particles",F=Point,C=Purple'
           write(1,1001) (x(d, i), d=1,SPH_dim), (vx(d, i), d = 1, SPH_dim),mass(i),rho(i), mu(i),temp(i),(0, d=1,SPH_dim) 
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
            write(1,1001) (x(d, i), d=1,SPH_dim), (vx(d, i), d = 1, SPH_dim),mass(i),rho(i), mu(i), temp(i), (surf_norm(d,s), d=1,SPH_dim)
            text_print=.false.
        endif
    enddo

    text_print=.true.
! Export all Ghost particles/ponts.
    
    do i=1, ntotal
        if((itype(i) .gt. itype_virtual) .and. (mod(itype(i),itype_virtual) .eq. 0) ) then
            if(text_print)   write(1,'(A)')'ZONE T="Virtual reference points",F=Point,C=Black'
            write(1,1001) (x(d, i), d=1,SPH_dim), (0, d = 1, SPH_dim),0,0,0,0,(0, d=1,SPH_dim)   
            text_print=.false.
        endif
    enddo
       

!Here is a readable format in which the data is stored
1001  format(9(e22.10,2x))    ! i am using Maximum, it was 9 before
1002  FORMAT(7(2X, e22.10))    ! i am using Maximum, it was 9 before
      
      close(1)

end
      