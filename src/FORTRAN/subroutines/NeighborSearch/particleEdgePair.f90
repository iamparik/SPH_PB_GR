subroutine particleEdgePair(a,s,scale_k,khsml)
    use config_parameter, only:SPH_dim, etype_virtual, etype_periodic
    use particle_data ,   only: eniac,max_e_interaction, &     
        & epair_a,epair_s, x, x_ve, surf_norm, &
        & pBC_eniac, pBC_epair_a, pBC_epair_s, etype, edge
     
    implicit none
    integer(4), intent(in):: a,s
    real(8), intent(in):: scale_k, khsml
    logical :: interact_edge_particle
    integer(4) d
    real(8) :: err_tol

    interact_edge_particle = .false.
    err_tol= dble( khsml* 1.D-6)
            
    ! We use an algorithm to determine if the edge and particle interact,
    ! in other words,
    !the itneraction occurs only if particles kernel is truncated (atleast for symmetric kernels)
    if (SPH_dim .eq. 2) then
        call edge_particle_pair_2D(khsml,x(:,a),x_ve(:,edge(1,s)),x_ve(:,edge(2,s)),surf_norm(:,s),interact_edge_particle,err_tol)
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
                epair_a(eniac) = a
                epair_s(eniac) = s
                   
            else
                write(*,*) ' >>> ERROR <<< : Too many interactions'
                write(*,*) ' max_e_interaction needs to be more than ', max_e_interaction
                pause
            endif
        endif
            
        if(etype(s) .eq. etype_periodic) then
            pBC_eniac = pBC_eniac + 1 
            pBC_epair_a(pBC_eniac) = a
            pBC_epair_s(pBC_eniac) = s
        endif
                
            
    endif
    
end subroutine