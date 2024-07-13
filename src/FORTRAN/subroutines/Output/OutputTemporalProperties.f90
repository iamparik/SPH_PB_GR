!****************************************************************************
!
!  SUBROUTINE: outputTemporalProperties 
!
!****************************************************************************      
subroutine  outputTemporalProperties(iterStep)
    use config_parameter, only:dataOutputPath, SPH_dim, itype_real_max, itype_real_min, &
            & etype_real_max, etype_real_min, etype_virtual, itype_virtual, itype_periodic
    use particle_data ,   only: x, vx,mass, itype, etype, ntotal, etotal, mid_pt_for_edge, xstart
    
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


    if (iterStep .eq. 1) then
        call createOutputFolder(dataOutputPath)
        allocate (xstart(SPH_dim,nreal))
        xstart=0.D0
        xstart= x(:,1:nreal)
    endif
    
    TPD=0.D0
    delC_avg=0.D0
    do a=1,nreal
        TPD = TPD + norm2(xStart(:,a)-x(:,a))/nreal
        delC_avg= delC_avg + norm2(delC(:,a))/nreal
        KE= 0.5D0*mass(a)*norm2(vx(:,a))**2
    enddo
    
    
!     output result as a continous plot against time   
    if (iterStep .eq. 1) then
        write(x2name, '(A,A)') dataOutputPath,'/energy.dat'
        open (1, file = x2name)
        write(1,'(A)') 'variables = Iterations , TPD , delC_avg , KE'
        close(1)
    elseif (mod(iterStep,1).eq. 0) then
        write(x2name, '(A,A)') dataPackingPath,'/energy.dat'
        open (1, file = x2name,  position='append')
        write(1,*) iterStep, TPD, delC_avg, KE
        close(1)
    endif
    


end
      