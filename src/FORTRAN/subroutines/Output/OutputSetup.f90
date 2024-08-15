!****************************************************************************
!
!  SUBROUTINE: output_flow_simplified 
!
!  PURPOSE: This subroutine will output simple flow variables
!
!   CREATED:        09/13/2017       by  PARIKSHIT BOREGOWDA
!   Last Modified:  03/25/2024        by  PARIKSHIT BOREGOWDA 
!****************************************************************************      
subroutine  OutputSetup
    use config_parameter, only:dataOutputPath, SPH_dim, itype_real_max, itype_real_min, &
        & etype_real_max, etype_real_min, etype_virtual, itype_virtual, itype_periodic, &
        & save_step
    use particle_data ,   only: x, surf_norm, &
        & ntotal,etotal, itype, etype, ve_total, x_ve, vx_ve, mid_pt_for_edge
    
    implicit none
!----------------------------------------------------------------------           
!                   i:      A variable used as a counter in loop
!                   d:      A variable used to iterate dimensions
!                   xname:  Path of the file which needs to be created
!------------------------------------------------------------------------
    integer(4):: i,d,s,f_dim
    real(8):: dt
    character (40) :: xname, x2name 
    logical :: text_print

    call createOutputFolder(dataOutputPath)
    
    !call gamma_density_continuous_leroy
    
!     output result in tecplot format
      write(xname, '(A,A, I1, A)') dataOutputPath,'/t=000000', 0,'step_xv.dat'
    

!   Open the file with the path xname, call this as 1
    open (1, file = xname)

!  Write the title on the open file
    write (1, '(A)') 'title="Starting Config"'

!   Specify variables that need to be read by tecplot to correspond to the data fed
! for 3D output 
    if(SPH_dim.eq.3)    write(1,*) 'variables = x, y, z, surf_normX,surf_normY,surf_normZ  '
! for 2D output
    if(SPH_dim.eq.2)    write(1,*) 'variables = x, y, surf_normX,surf_normY'
! for 1D output
    if(SPH_dim.eq.1)    write(1,*) 'variables = x, surf_normX'
  
text_print=.true.
! First export data of all real particles ( not on boundary)
    
    do i=1,ntotal
       if((itype(i) .le. itype_real_max) .and. (itype(i) .gt. itype_real_min))  then
           if(text_print) write(1,'(A)')'ZONE T="Real Particles",F=Point,C=Blue'
           write(1,1001) (x(d, i), d=1,SPH_dim),(0.D0, d = 1, SPH_dim)
           text_print=.false. 
       endif       
    enddo


    text_print=.true.
! Export data of all particles representative of boundary, which are used in interpolation,
! called edge particles
    
    do s=1,etotal 
        if(text_print) write(1,'(A)')'ZONE T="Edge Mid-Points",F=Point,C=Red'
        write(1,1001) (mid_pt_for_edge(d,s), d=1,SPH_dim), (surf_norm(d,s), d = 1, SPH_dim) 
        text_print=.false.                           
    enddo

    text_print=.true.       
    
! Export all vertices of edges.
    
    do i=1, ve_total
        if(text_print)  write(1,'(A)')'ZONE T="Edge vertices",F=Point,C=Black'
        write(1,1001) (x_ve(d, i), d=1,SPH_dim),(0.D0, d = 1, SPH_dim)    
        text_print=.false.
    enddo
   
    text_print=.true.
! Export data of all copied real particles due to periodic bc (not on boundary)
   

!Here is a readable format in which the data is stored
1001 format(14(e22.10,2x))    ! i am using Maximum, it was 9 before


1002  FORMAT(8(2X, e22.10))    ! i am using Maximum, it was 9 before
      
      close(1)

end
      