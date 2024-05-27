!****************************************************************************
!
!  SUBROUTINE: Output_factors 
!
!  PURPOSE: This subroutine will output the correction factors defined,
!
!   CREATED:        05/03/2022       by  PARIKSHIT BOREGOWDA
!   Last Modified:  05/03/2022       by  PARIKSHIT BOREGOWDA 
!****************************************************************************      
subroutine  outputFncnApprx(n,x, fncn0, fncn, dfncn0, dfncn, lap_fncn0, lap_fncn, nedge_rel_edge, edge, etotal, SPH_dim, ntotal,itype,etype,dataOutputPath)
     use config_parameter, only: itype_real_max, itype_real_min, etype_real_max, etype_real_min, &
        & etype_virtual, itype_virtual, itype_periodic
    implicit none
!----------------------------------------------------------------------           
!                   i:      A variable used as a counter in loop
!                   d:      A variable used to iterate dimensions
!                   xname:  Path of the file which needs to be created
!------------------------------------------------------------------------
    integer(4), intent(in):: SPH_dim, ntotal, n, etotal
    integer(2), intent(in), dimension(ntotal):: itype
    integer(4), intent(in), dimension(etotal):: etype
    real(8), intent(in), dimension(ntotal):: fncn, fncn0, lap_fncn, lap_fncn0
    real(8), intent(in), dimension(SPH_dim, ntotal):: dfncn0, dfncn, x
    integer(4), intent(in), dimension(SPH_dim, etotal):: edge
    integer(4), intent(in), dimension(etotal):: nedge_rel_edge
    integer(4):: a,d,s  
    character(LEN=*) :: dataOutputPath
    character (40) :: xname
    logical :: text_print

    
    
!     output result in tecplot format

!   Write the path of the file (with file name) to be edited
    write(xname, '(A,A, I2, A)') dataOutputPath,'/n=',n,'step.dat'

!   Open the file with the path xname, call this as 1
    open (1, file = xname)

!  Write the title on the open file
    write (1, '(A,I2,A)') 'title="n=',n,'time-step"'
    
 

!   Specify variables that need to be read by tecplot to correspond to the data fed
! for 3D output 
    if(SPH_dim.eq.3)    write(1,'(A)') 'variables ="x", "y", "z", "f0", "f", "E(f)", "df0(x)", "df0(y)", "df0(z)", "df(x)", "df(y)", "df(z)", &
                                &  "Edf(x)", "Edf(y)", "Edf(z)","d2f0","d2f", "Ed2f" '
! for 2D output
    if(SPH_dim.eq.2)    write(1,'(A)') 'variables="x", "y", "f0", "f", "E(f)", "df0(x)", "df0(y)","df(x)", "df(y)", "Edf(x)", "Edf(y)",  "d2f0","d2f", "Ed2f"  '
! for 1D output
    if(SPH_dim.eq.1)    write(1,'(A)') 'variables = "x", "f0", "f","E(f)", "df0(x)", "df(x)","Edf(x)", "d2f0","d2f", "Ed2f"'
  
   text_print=.true.
! First export data of all real particles ( not on boundary)
    
    do a=1,ntotal
        if((itype(a) .le. itype_real_max) .and. (itype(a) .gt. itype_real_min)) then
            if(text_print) write(1,'(A)')'ZONE T="Real Particles",F=Point,C=Blue'
            write(1,1001) (x(d, a), d=1,SPH_dim), fncn0(a), fncn(a),abs(fncn0(a)-fncn(a)), (dfncn0(d, a), d=1,SPH_dim), &
                        & (dfncn(d, a), d=1,SPH_dim), (abs(dfncn0(d,a)-dfncn(d,a)), d=1,SPH_dim) ,lap_fncn0(a), lap_fncn(a),&
                        & abs(lap_fncn0(a)-lap_fncn(a)) 
            
            text_print=.false. 
        endif
        
    enddo
    
   text_print=.true.
! Export data of all marker points used to capture data    
    do a=1,ntotal
        if(itype(a) .eq. itype_real_min) then
            if(text_print) write(1,'(A)')'ZONE T="Marker points",F=Point,C=Cyan'
            write(1,1001) (x(d, a), d=1,SPH_dim), fncn0(a), fncn(a),abs(fncn0(a)-fncn(a)), (dfncn0(d, a), d=1,SPH_dim), &
                        & (dfncn(d, a), d=1,SPH_dim), (abs(dfncn0(d,a)-dfncn(d,a)), d=1,SPH_dim) ,lap_fncn0(a), lap_fncn(a),&
                        & abs(lap_fncn0(a)-lap_fncn(a)) 
            
            text_print=.false. 
        endif
        
    enddo

    text_print=.true.
! Export data of all particles representative of boundary, which are used in interpolation,
! called edge particles
    
    do s=1,etotal 
        if((etype(s) .le. etype_real_max) .and. (etype(s) .ge. etype_real_min)) then
            if(text_print)  write(1,'(A)')'ZONE T="Edge Particles",F=Point,C=Red'
	        a=nedge_rel_edge(s)  ! this needs to be changed for edges with mroe than one poitn representation. This will then be itereated
            write(1,1001) (x(d, a), d=1,SPH_dim), fncn0(a), fncn(a),abs(fncn0(a)-fncn(a)), (dfncn0(d, a), d=1,SPH_dim), &
                        & (dfncn(d, a), d=1,SPH_dim), (abs(dfncn0(d,a)-dfncn(d,a)), d=1,SPH_dim) ,lap_fncn0(a), lap_fncn(a),&
                        & abs(lap_fncn0(a)-lap_fncn(a)) 
            text_print=.false. 
        endif
    enddo

    text_print=.true.
! Export all edge vertices
    
    do a=1,ntotal
        if( itype(a) .eq. 0) then
            if(text_print) write(1,'(A)')'ZONE T="Edge vertices",F=Point,C=Black'
            write(1,1001) (x(d, a), d=1,SPH_dim), fncn0(a), fncn(a),abs(fncn0(a)-fncn(a)), (dfncn0(d, a), d=1,SPH_dim), &
                        & (dfncn(d, a), d=1,SPH_dim), (abs(dfncn0(d,a)-dfncn(d,a)), d=1,SPH_dim) ,lap_fncn0(a), lap_fncn(a),&
                        & abs(lap_fncn0(a)-lap_fncn(a))    
            text_print=.false. 
        endif
    enddo
    
    text_print=.true.
! Export data of any copied real particles (not on boundary) - which serve as virtual domain
    
    do a=1,ntotal
        if(.not.((itype(a) .le. itype_real_max) .and. (itype(a) .ge. itype_real_min)) .and. (mod(itype(a),itype_virtual) .ne. 0) ) then
            if(text_print) write(1,'(A)')'ZONE T="Virtual Domain Particles",F=Point,C=Blue'
            write(1,1001) (x(d, a), d=1,SPH_dim), fncn0(a), fncn(a),abs(fncn0(a)-fncn(a)), (dfncn0(d, a), d=1,SPH_dim), &
                        & (dfncn(d, a), d=1,SPH_dim), (abs(dfncn0(d,a)-dfncn(d,a)), d=1,SPH_dim) ,lap_fncn0(a), lap_fncn(a),&
                        & abs(lap_fncn0(a)-lap_fncn(a))  
            text_print=.false. 
        endif

    enddo

    text_print=.true.
! Edge particles/reference points that are represent the boudnary of the virtual domain           
    do s=1,etotal 
        if(.not.((etype(s) .le. etype_real_max) .and. (etype(s) .ge. etype_real_min))) then
            if(text_print)   write(1,'(A)')'ZONE T="Virtual Edge Particles",F=Point,C=Red'
	        a=nedge_rel_edge(s)  ! this needs to be changed for edges with mroe than one poitn representation. This will then be itereated
            write(1,1001) (x(d, a), d=1,SPH_dim), fncn0(a), fncn(a),abs(fncn0(a)-fncn(a)), (dfncn0(d, a), d=1,SPH_dim), &
                        & (dfncn(d, a), d=1,SPH_dim), (abs(dfncn0(d,a)-dfncn(d,a)), d=1,SPH_dim) ,lap_fncn0(a), lap_fncn(a),&
                        & abs(lap_fncn0(a)-lap_fncn(a)) 
            text_print=.false. 
        endif
   
    enddo

   text_print=.true.
! Export all Ghost particles/ponts.      
    do a=1,ntotal
        if( (itype(a) .gt. itype_virtual) .and. (mod(itype(a),itype_virtual) .eq. 0)) then
            if(text_print)  write(1,'(A)')'ZONE T="Virtual reference points",F=Point,C=Black'
            write(1,1001) (x(d, a), d=1,SPH_dim), fncn0(a), fncn(a),abs(fncn0(a)-fncn(a)), (dfncn0(d, a), d=1,SPH_dim), &
                        & (dfncn(d, a), d=1,SPH_dim), (abs(dfncn0(d,a)-dfncn(d,a)), d=1,SPH_dim) ,lap_fncn0(a), lap_fncn(a),&
                        & abs(lap_fncn0(a)-lap_fncn(a)) 
            text_print=.false. 
        endif
        
    enddo
      
!Here is a readable format in which the data is stored
1001  format(12(e22.10,2x))    ! i am using Maximum, it was 9 before
! 2X means skip 2 spaces
      
      close(1)

end
      