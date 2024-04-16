!****************************************************************************
!
!  SUBROUTINE: Output_factors 
!
!  PURPOSE: This subroutine will output the correction factors defined,
!
!   CREATED:        04/18/2022       by  PARIKSHIT BOREGOWDA
!   Last Modified:  04/22/2022        by  PARIKSHIT BOREGOWDA 
!****************************************************************************      
subroutine  output_factors
    use config_parameter, only:dataOutputErrors, SPH_dim, itype_real_max, itype_real_min, &
        & etype_real_max, etype_real_min, etype_virtual, itype_virtual, itype_periodic
    use particle_data ,   only: x, gamma_discrt,nreal, nedge, nghost, ntotal,nedge_rel_edge, &
                    & edge,epair_a,epair_s, eniac, del_gamma, gamma_cont, xi1_mat,xi1_mat_inv,&
                    & beta_mat, gamma_mat, gamma_mat_inv, xi_cont_mat, xi_cont_mat_inv, etotal,&
                    & itype, etype
    
    
    implicit none
!----------------------------------------------------------------------           
!                   i:      A variable used as a counter in loop
!                   d:      A variable used to iterate dimensions
!                   xname:  Path of the file which needs to be created
!------------------------------------------------------------------------
    integer(4) i,j,d,k,s,a
    integer(4) pairs(ntotal)
    character (40) :: xname
    logical :: text_print
    

    pairs=0
    do k=1,eniac
        a=epair_a(k)
        s=epair_s(k)
        
        pairs(a)=pairs(a)+1
        
    enddo
    
    
!     output result in tecplot format

!   Write the path of the file (with file name) to be edited
    write(xname, '(A,A)') dataOutputErrors,'/factors.dat'

!   Open the file with the path xname, call this as 1
    open (1, file = xname)

!  Write the title on the open file
    write (1, '(A)') 'title="SPH Correction factors"'
    
 

!   Specify variables that need to be read by tecplot to correspond to the data fed
! for 3D output 
    if(SPH_dim.eq.3)    write(1,'(A)') 'variables ="x", "y", "z", "alpha" ,"truncted_num", "dgammax", "dgammy", "dgammz", "gamma",&
                                    & "xi(1,1)", "xi(1,2)", "xi(1,3)", "xi(2,1)", "xi(2,2)", "xi(2,3)", "xi(3,1)", "xi(3,2)", "xi(3,3)",&
                                    & "beta(1,1)", "beta(1,2)","beta(1,3)","beta(2,1)", "beta(2,2)","beta(2,3)","beta(3,1)", "beta(3,2)","beta(3,3)",&
                                    & "GAMMA(1,1)", "GAMMA(1,2)", "GAMMA(1,3)", "GAMMA(2,1)", "GAMMA(2,2)", "GAMMA(2,3)", "GAMMA(3,1)", "GAMMA(3,2)", "GAMMA(3,3)",&
                                    & "xi_cont(1,1)", "xi_cont(1,2)", "xi_cont(1,3)", "xi_cont(2,1)", "xi_cont(2,2)", "xi_cont(2,3)", "xi_cont(3,1)", "xi_cont(3,2)", "xi_cont(3,3)",&
                                    & "xi_inv(1,1)", "xi_inv(1,2)", "xi_inv(1,3)", "xi_inv(2,1)", "xi_inv(2,2)", "xi_inv(2,3)", "xi_inv(3,1)", "xi_inv(3,2)", "xi_inv(3,3)",&
                                    & "GAMMA_inv(1,1)", "GAMMA_inv(1,2)", "GAMMA_inv(1,3)", "GAMMA_inv(2,1)", "GAMMA_inv(2,2)", "GAMMA_inv(2,3)", "GAMMA_inv(3,1)", "GAMMA_inv(3,2)", "GAMMA_inv(3,3)",&
                                    & "xi_cont_inv(1,1)", "xi_cont_inv(1,2)", "xi_cont_inv(1,3)", "xi_cont_inv(2,1)", "xi_cont_inv(2,2)", "xi_cont_inv(2,3)", "xi_cont_inv(3,1)", "xi_cont_inv(3,2)", "xi_cont_inv(3,3)"'
! for 2D output
    if(SPH_dim.eq.2)    write(1,'(A)') 'variables="x", "y", "alpha", "truncted_num", "dgammax", "dgammy", "gamma",&
                                    & "xi(1,1)", "xi(1,2)", "xi(2,1)", "xi(2,2)"&
                                    & "beta(1,1)", "beta(1,2)","beta(2,1)", "beta(2,2)",&
                                    & "GAMMA(1,1)", "GAMMA(1,2)", "GAMMA(2,1)", "GAMMA(2,2)",&
                                    & "xi_cont(1,1)", "xi_cont(1,2)", "xi_cont(2,1)", "xi_cont(2,2)"&
                                    & "xi_inv(1,1)", "xi_inv(1,2)", "xi_inv(2,1)", "xi_inv(2,2)",&
                                    & "GAMMA_inv(1,1)", "GAMMA_inv(1,2)","GAMMA_inv(2,1)", "GAMMA_inv(2,2)",&
                                    & "xi_cont_inv(1,1)", "xi_cont_inv(1,2)", "xi_cont_inv(2,1)", "xi_cont_inv(2,2)"'
! for 1D output
    if(SPH_dim.eq.1)    write(1,'(A)') 'variables = "x", "alpha", "truncted_num", "dgammax", "gamma", "xi(1,1)", "beta(1,1)", "GAMMA(1,1)","xi_cont(1,1)","xi_inv(1,1)","GAMMA_inv(1,1)","xi_cont_inv(1,1)"'
  
   text_print=.true.
! First export data of all real particles ( not on boundary)
    do a=1,ntotal
        if ((itype(a) .le. itype_real_max) .and. (itype(a) .gt. itype_real_min)) then
            if(text_print) write(1,'(A)')'ZONE T="Real Particles",F=Point,C=Blue'
            write(1,1001) (x(d, a), d=1,SPH_dim), gamma_discrt(a),pairs(a), (del_gamma(d, a), d=1,SPH_dim), gamma_cont(a), &
                            & ((xi1_mat(i,j,a), j=1,SPH_dim), i=1,SPH_dim),((beta_mat(i,j,a), j=1,SPH_dim), i=1,SPH_dim),&
                            & ((gamma_mat(i,j,a), j=1,SPH_dim), i=1,SPH_dim), ((xi_cont_mat(i,j,a), j=1,SPH_dim), i=1,SPH_dim), & 
                            & ((xi1_mat_inv(i,j,a), j=1,SPH_dim), i=1,SPH_dim), ((gamma_mat_inv(i,j,a), j=1,SPH_dim), i=1,SPH_dim), &
                            & ((xi_cont_mat_inv(i,j,a), j=1,SPH_dim), i=1,SPH_dim)
            text_print=.false.
        endif
    enddo
    
    text_print=.true.
! Export data of all marker points used to capture data
    do a=1,ntotal
        if (itype(a) .eq. itype_real_min)then
            if(text_print) write(1,'(A)')'ZONE T="Marker points",F=Point,C=Cyan'
            write(1,1001) (x(d, a), d=1,SPH_dim), gamma_discrt(a),pairs(a), (del_gamma(d, a), d=1,SPH_dim), gamma_cont(a), &
                            & ((xi1_mat(i,j,a), j=1,SPH_dim), i=1,SPH_dim),((beta_mat(i,j,a), j=1,SPH_dim), i=1,SPH_dim),&
                            & ((gamma_mat(i,j,a), j=1,SPH_dim), i=1,SPH_dim), ((xi_cont_mat(i,j,a), j=1,SPH_dim), i=1,SPH_dim), & 
                            & ((xi1_mat_inv(i,j,a), j=1,SPH_dim), i=1,SPH_dim), ((gamma_mat_inv(i,j,a), j=1,SPH_dim), i=1,SPH_dim), &
                            & ((xi_cont_mat_inv(i,j,a), j=1,SPH_dim), i=1,SPH_dim)
            text_print=.false.
        endif
    enddo

    text_print=.true.
! Export data of all particles representative of boundary, which are used in interpolation,
! called edge particles
    do s=1,etotal 
        if((etype(s) .le. etype_real_max) .and. (etype(s) .ge. etype_real_min)) then
            if(text_print) write(1,'(A)')'ZONE T="Edge Particles",F=Point,C=Red'
	        a=nedge_rel_edge(s)  ! this needs to be changed for edges with mroe than one poitn representation. This will then be itereated
             write(1,1001) (x(d, a), d=1,SPH_dim), gamma_discrt(a),pairs(a), (del_gamma(d, a), d=1,SPH_dim), gamma_cont(a), &
                            & ((xi1_mat(i,j,a), j=1,SPH_dim), i=1,SPH_dim),((beta_mat(i,j,a), j=1,SPH_dim), i=1,SPH_dim),&
                            & ((gamma_mat(i,j,a), j=1,SPH_dim), i=1,SPH_dim), ((xi_cont_mat(i,j,a), j=1,SPH_dim), i=1,SPH_dim), & 
                            & ((xi1_mat_inv(i,j,a), j=1,SPH_dim), i=1,SPH_dim), ((gamma_mat_inv(i,j,a), j=1,SPH_dim), i=1,SPH_dim), &
                            & ((xi_cont_mat_inv(i,j,a), j=1,SPH_dim), i=1,SPH_dim) 
             text_print=.false.
        endif
    enddo


    text_print=.true.
! Export all edge vertices.
    do a=1,ntotal
        if (itype(a) .eq. 0) then
            if(text_print)  write(1,'(A)')'ZONE T="Edge vertices",F=Point,C=Black'
            write(1,1001) (x(d, a), d=1,SPH_dim), gamma_discrt(a),pairs(a), (del_gamma(d, a), d=1,SPH_dim), gamma_cont(a), &
                            & ((xi1_mat(i,j,a), j=1,SPH_dim), i=1,SPH_dim),((beta_mat(i,j,a), j=1,SPH_dim), i=1,SPH_dim),&
                            & ((gamma_mat(i,j,a), j=1,SPH_dim), i=1,SPH_dim), ((xi_cont_mat(i,j,a), j=1,SPH_dim), i=1,SPH_dim), & 
                            & ((xi1_mat_inv(i,j,a), j=1,SPH_dim), i=1,SPH_dim), ((gamma_mat_inv(i,j,a), j=1,SPH_dim), i=1,SPH_dim), &
                            & ((xi_cont_mat_inv(i,j,a), j=1,SPH_dim), i=1,SPH_dim)
            text_print=.false.
        endif
    enddo
    
    text_print=.true.
 
! Export data of any copied real particles (not on boundary) - which serve as virtual domain
    do a=1,ntotal
        if (.not.((itype(a) .le. itype_real_max) .and. (itype(a) .ge. itype_real_min)).and. (mod(itype(a),itype_virtual) .ne. 0)) then
            if(text_print) write(1,'(A)')'ZONE T="Virtual Domain Particles",F=Point,C=Blue'
            write(1,1001) (x(d, a), d=1,SPH_dim), gamma_discrt(a),pairs(a), (del_gamma(d, a), d=1,SPH_dim), gamma_cont(a), &
                            & ((xi1_mat(i,j,a), j=1,SPH_dim), i=1,SPH_dim),((beta_mat(i,j,a), j=1,SPH_dim), i=1,SPH_dim),&
                            & ((gamma_mat(i,j,a), j=1,SPH_dim), i=1,SPH_dim), ((xi_cont_mat(i,j,a), j=1,SPH_dim), i=1,SPH_dim), & 
                            & ((xi1_mat_inv(i,j,a), j=1,SPH_dim), i=1,SPH_dim), ((gamma_mat_inv(i,j,a), j=1,SPH_dim), i=1,SPH_dim), &
                            & ((xi_cont_mat_inv(i,j,a), j=1,SPH_dim), i=1,SPH_dim)
            text_print=.false.
        endif
    enddo

    text_print=.true.
! Edge particles/reference points that are represent the boudnary of the virtual domain       
    do s=1,etotal 
         if(.not.((etype(s) .le. etype_real_max) .and. (etype(s) .ge. etype_real_min))) then
            if(text_print)  write(1,'(A)')'ZONE T="Virtual Edge Particles",F=Point,C=Red'
	        a=nedge_rel_edge(s)  ! this needs to be changed for edges with mroe than one poitn representation. This will then be itereated
            write(1,1001) (x(d, a), d=1,SPH_dim), gamma_discrt(a),pairs(a), (del_gamma(d, a), d=1,SPH_dim), gamma_cont(a), &
                            & ((xi1_mat(i,j,a), j=1,SPH_dim), i=1,SPH_dim),((beta_mat(i,j,a), j=1,SPH_dim), i=1,SPH_dim),&
                            & ((gamma_mat(i,j,a), j=1,SPH_dim), i=1,SPH_dim), ((xi_cont_mat(i,j,a), j=1,SPH_dim), i=1,SPH_dim), & 
                            & ((xi1_mat_inv(i,j,a), j=1,SPH_dim), i=1,SPH_dim), ((gamma_mat_inv(i,j,a), j=1,SPH_dim), i=1,SPH_dim), &
                            & ((xi_cont_mat_inv(i,j,a), j=1,SPH_dim), i=1,SPH_dim) 
            text_print=.false.
         endif
    enddo

    text_print=.true.
! Export all virtual reference points.       
    do a=1,ntotal
        if ((itype(a) .gt. itype_virtual) .and. (mod(itype(a),itype_virtual) .eq. 0) ) then
            if(text_print)  write(1,'(A)')'ZONE T="Virtual reference points",F=Point,C=Black'
            write(1,1001) (x(d, a), d=1,SPH_dim), gamma_discrt(a),pairs(a), (del_gamma(d, a), d=1,SPH_dim), gamma_cont(a), &
                            & ((xi1_mat(i,j,a), j=1,SPH_dim), i=1,SPH_dim),((beta_mat(i,j,a), j=1,SPH_dim), i=1,SPH_dim),&
                            & ((gamma_mat(i,j,a), j=1,SPH_dim), i=1,SPH_dim), ((xi_cont_mat(i,j,a), j=1,SPH_dim), i=1,SPH_dim), & 
                            & ((xi1_mat_inv(i,j,a), j=1,SPH_dim), i=1,SPH_dim), ((gamma_mat_inv(i,j,a), j=1,SPH_dim), i=1,SPH_dim), &
                            & ((xi_cont_mat_inv(i,j,a), j=1,SPH_dim), i=1,SPH_dim)
            text_print=.false.
        endif
    enddo

      
!Here is a readable format in which the data is stored
1001  format(3(e22.10,2x),(I5,2X),23(e22.10,2x))    ! i am using Maximum, it was 9 before
! 2X means skip 2 spaces
      
      close(1)

end
      