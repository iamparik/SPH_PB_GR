!****************************************************************************
!
!  SUBROUTINE: caseBasedOutput 
!
!****************************************************************************      
subroutine  caseBasedOutput(iterStep,dt)
    use config_parameter, only:dataOutputPath, SPH_dim
    use particle_data ,   only: x, vx,mass, nreal, xstart, delC, prsr_bdry_val, &
            & edge,etotal, x_ve, etype
    
    implicit none
!----------------------------------------------------------------------           
!                   i:      A variable used as a counter in loop
!                   d:      A variable used to iterate dimensions
!                   xname:  Path of the file which needs to be created
!------------------------------------------------------------------------
    integer(4), intent(in) :: iterStep
    real(8), intent(in) :: dt    
    integer(4):: a,d, s
    character (40) ::x2name
    real(8) :: v_max, PE, p_lft_wedge,KE

    if (iterStep .eq. 1) call createOutputFolder(dataOutputPath)

    
    PE=0.D0
    p_lft_wedge=0.D0
    KE=0.D0
    v_max=0.D0
    
    do a=1,nreal
        PE= PE + mass(a)*9.81*x(2,a)
        v_max = max(norm2(vx(:,a)),v_max)
        KE= KE+ 0.5D0*mass(a)*norm2(vx(:,a))**2
    enddo
    
    do s=1,etotal
        if(etype(s) .eq. 2) p_lft_wedge=p_lft_wedge + prsr_bdry_val(s)*norm2(x_ve(:,edge(1,s))-x_ve(:,edge(2,s))) 
    enddo
    
    
!     output result as a continous plot against time   
    if (iterStep .eq. 1) then
        write(x2name, '(A,A)') dataOutputPath,'/energy.dat'
        open (1, file = x2name)
        write(1,'(A)') 'variables = time , PE , p_lft_wedge , KE, v_max'
        write(1,*) iterStep*dt, PE, p_lft_wedge, KE, v_max
        close(1)
    else
        write(x2name, '(A,A)') dataOutputPath,'/energy.dat'
        open (1, file = x2name,  position='append')
        write(1,*) iterStep*dt, PE, p_lft_wedge, KE, v_max
        close(1)
    endif
    


end
      