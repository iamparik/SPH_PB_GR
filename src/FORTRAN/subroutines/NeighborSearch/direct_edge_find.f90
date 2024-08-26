!****************************************************************************
!
!  SUBROUTINE: direct_edge_find 
!
!  PURPOSE:  For every edge, all particles are checked to see if  they
!           are within interactive distance such that the kernel of particle
!           is truncated at the boundary.
!
!   CREATED:        04/21/2022       by  PARIKSHIT BOREGOWDA 
!   Last Modified:  04/21/2022       by  PARIKSHIT BOREGOWDA 
!****************************************************************************

subroutine direct_edge_find
    use config_parameter, only: SPH_dim, etype_virtual, itype_virtual, etype_periodic
    use particle_data ,   only: eniac,edge,ntotal,itype, &  
        & hsml, surf_norm, epair_a,epair_s,max_e_interaction, x , etotal, &
        & etype, pBC_epair_a, pBC_epair_s, pBC_eniac, x_ve
        
    implicit none

!----------------------------------------------------------------------
!                   i,j: A variable used as a counter in loop
!                   d: A variable used to iterate dimensions
!                   interact_edge_particle: A variable when true implies, an
!                               edgeand particle will interact.
!
!----------------------------------------------------------------------

    integer(4) i,s,d, n_edges   
    real(8) scale_k,khsml, driac,dxiac(SPH_dim), tdwdx(SPH_dim) ,err_tol

!  The maximum relative distance is determined, after which the kernel will be zero,
! scale_k is not to be confused with kappa which modifies smoothening length
! scale_k is only the factor used in the kernel to approximate various sections
! as 0<R<1, 1<R<2 will all have different approximations.
    call sml_mult_factor(scale_k)

! We initailize the variable storing the number of interations between edges and particles 
    eniac=0
    
    pbc_eniac=0
    
! We retrieve the total number of edges in the computation as n_edges
    n_edges=etotal
    
! Now we check every possible pair of particles in the computation domain 
! and determine particle pairs who lie within each others compact kernel   
! The below is only for cosntant smoothing length
    do i=1,ntotal
        ! Only particles or points that are used in SPH integration/summation is to be inlcuded
        ! all itypes 0 are classified as reference points that are not used in simulation or for boundary condition
        if (mod(itype(i),itype_virtual) .ne. 0) then
            do s = 1, n_edges
                ! Here smoothing length of ith particle needs to be considered for compact distance 
                ! as its truncated kernel gives raise to boundary integral formulation
                khsml = scale_k*hsml(i)
            
                call particleEdgePair(i,s,scale_k,khsml)
            
            enddo
        endif
        
   enddo
      

end