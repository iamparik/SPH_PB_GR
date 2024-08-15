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

subroutine direct_edge_find_reduced(sel_part_list, n_reduced)
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

    integer(4) i,s,d, n_edges, ii   
    real(8) scale_k,khsml, driac,dxiac(SPH_dim), tdwdx(SPH_dim) ,err_tol
    logical :: interact_edge_particle
    integer(4), intent(in):: sel_part_list(ntotal), n_reduced    
    integer(4) max_pBC_e_interaction

    max_pBC_e_interaction=int(40*max_e_interaction/etotal)
    
    if (.NOT. Allocated(epair_a)) allocate(epair_a(max_e_interaction)) 
    if (.NOT. Allocated(epair_s)) allocate(epair_s(max_e_interaction)) 

    if (.NOT. Allocated(pBC_epair_a)) allocate(pbc_epair_a(max_pBC_e_interaction)) 
    if (.NOT. Allocated(pBC_epair_s)) allocate(pbc_epair_s(max_pBC_e_interaction)) 

    epair_a=0
    epair_s=0

    pBC_epair_a=0
    pBC_epair_s=0
    
! initialize the interaction logical coutner as false to start
    interact_edge_particle= .false.
    
    
    
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
    do i=1,n_reduced
        ii= sel_part_list(i)
        ! Only particles or points that are used in SPH integration/summation is to be inlcuded
        ! all itypes 0 are classified as reference points that are not used in simulation or for boundary condition
        if (mod(itype(ii),itype_virtual) .eq. 0) go to 2
        do s = 1, n_edges
            ! Here smoothing length of ith particle needs to be considered for compact distance 
            ! as its truncated kernel gives raise to boundary integral formulation
            khsml = scale_k*hsml(ii)
            
            err_tol= dble( khsml* 1.D-6)
            
            ! We use an algorithm to determine if the edge and particle interact,
            ! in other words,
            !the itneraction occurs only if particles kernel is truncated (atleast for symmetric kernels)
            if (SPH_dim .eq. 2) then
                call edge_particle_pair_2D(khsml,x(:,ii),x_ve(:,edge(1,s)),x_ve(:,edge(2,s)),surf_norm(:,s),interact_edge_particle,err_tol)
            else
                write(*,*) ' >>> ERROR <<< : Edge_particle pairing algorithm not defined for dimension', SPH_dim
                pause    
            endif
            
            if (interact_edge_particle) then
                if (mod(etype(s),etype_virtual) .ne. 0) then 
                    if (eniac .lt. max_e_interaction) then
                        !     Neighboring pair list, and totalinteraction number and
                        !     the interaction number for each particle 
                        eniac = eniac + 1
                        epair_a(eniac) = ii
                        epair_s(eniac) = s
                   
                    else
                        write(*,*) ' >>> ERROR <<< : Too many interactions'
                        write(*,*) ' max_e_interaction needs to be more than ', max_e_interaction
                        pause
                    endif
                endif
            
                if(etype(s) .eq. etype_periodic) then
                    pBC_eniac = pBC_eniac + 1 
                    pBC_epair_a(pBC_eniac) = ii
                    pBC_epair_s(pBC_eniac) = s
                endif
                
            
            endif
            
1       enddo
2   enddo
      

end