!****************************************************************************
!
!  SUBROUTINE: outputPacking 
!
!  PURPOSE: This subroutine will output the input geomtrical data,
!           and can be used to test if the input file is being
!           read correctly for only thermal diffusion
!
!   CREATED:        09/13/2017       by  PARIKSHIT BOREGOWDA
!   Last Modified:  10/17/2023        by  PARIKSHIT BOREGOWDA 
!****************************************************************************      
subroutine  outputPacking(iterStep, saveStep, TPD)
    use config_parameter, only:dataPackingPath, SPH_dim, itype_real_max, itype_real_min, &
            & etype_real_max, etype_real_min, etype_virtual, itype_virtual, itype_periodic
    use particle_data ,   only: x, itype, etype, nedge_rel_edge, ntotal, etotal
    
    implicit none
!----------------------------------------------------------------------           
!                   i:      A variable used as a counter in loop
!                   d:      A variable used to iterate dimensions
!                   xname:  Path of the file which needs to be created
!------------------------------------------------------------------------
    integer(4):: i,d,s, iterStep,f_dim,saveStep
    character (40) :: xname, x2name
    real(8) :: TPD
    logical :: text_print

    
!     output result as a continous plot against time   
    if (iterStep .eq. 1) then
        write(x2name, '(A,A)') dataPackingPath,'/energy.dat'
        open (1, file = x2name)
        write(1,'(A)') 'variables = Iterations, TPD '
        close(1)
    elseif (mod(iterStep,100).eq. 0) then
        write(x2name, '(A,A)') dataPackingPath,'/energy.dat'
        open (1, file = x2name,  position='append')
        write(1,*) iterStep, TPD
        close(1)
    endif
    
! output resultsof each time step as particle position
    if ((iterStep .eq. 1) .or. (mod(iterStep,saveStep).eq.0)) then
        !     output result in tecplot format
        if(iterStep.eq.0) then
            write(xname, '(A,A, I1, A)') dataPackingPath,'/t=000000', iterStep,'step_xv.dat'
        else if(iterStep>0.and.iterStep<10) then
            write(xname, '(A,A, I1, A)') dataPackingPath,'/t=000000', iterStep,'step_xv.dat'
        else if(iterStep>9.and.iterStep<100) then
            write(xname, '(A,A, I2, A)') dataPackingPath,'/t=00000', iterStep,'step_xv.dat'
        else if(iterStep>99.and.iterStep<1000) then
            write(xname, '(A,A, I3, A)') dataPackingPath,'/t=0000', iterStep,'step_xv.dat'
        else if(iterStep>999.and.iterStep<10000) then                                   
            write(xname, '(A,A, I4, A)') dataPackingPath,'/t=000', iterStep,'step_xv.dat'
        else if(iterStep>9999.and.iterStep<100000) then
            write(xname, '(A,A, I5, A)') dataPackingPath,'/t=00', iterStep,'step_xv.dat'
        else if(iterStep>99999.and.iterStep<1000000) then
            write(xname, '(A,A, I6, A)') dataPackingPath,'/t=0', iterStep,'step_xv.dat'
        else if(iterStep>999999.and.iterStep<10000000) then
            write(xname, '(A,A, I7, A)') dataPackingPath,'/t=', iterStep,'step_xv.dat'
        endif

 
    !   Open the file with the path xname, call this as 1
        open (1, file = xname)

    !  Write the title on the open file
        write (1, '(A)') 'title="Starting Config"'

    !   Specify variables that need to be read by tecplot to correspond to the data fed
        write(1,*) 'variables = x, y'
        
        text_print=.true.
    
        ! First export data of all real particles ( not on boundary)
        do a = 1,ntotal
            if((itype(i) .le. itype_real_max) .and. (itype(i) .gt. itype_real_min))  then
                if(text_print) write(1,'(A)')'ZONE T="Real Particles",F=Point,C=Blue'
               write(1,*) (x(d, a), d=1,SPH_dim)
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
                write(1,*) (x(d, a), d=1,SPH_dim)
                text_print=.false.
            endif                                       
        enddo

        text_print=.true.       
    
    ! Export all Ghost particles/ponts.    
        do i=1, ntotal
            if(itype(i) .eq. 0) then
                if(text_print)  write(1,'(A)')'ZONE T="Edge vertices",F=Point,C=Black'
                write(1,*) (x(d, a), d=1,SPH_dim)
                text_print=.false.
            endif
        enddo
    endif
      close(1)

end
      