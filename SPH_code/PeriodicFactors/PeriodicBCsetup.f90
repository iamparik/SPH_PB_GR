!****************************************************************************
!
!  SUBROUTINE: PeriodicBCsetup
!
!  PURPOSE: For simulations having periodic boundary condition
!           the periodicity is set up in this subroutine 
!
!   CREATED:        07/25/2022       by  PARIKSHIT BOREGOWDA
!   Last Modified:  08/02/2022       by  PARIKSHIT BOREGOWDA 
!****************************************************************************

!this subroutine needs to further generalized to be applied in other scenarios
subroutine PeriodicBCsetup2D 
use config_parameter, only:SPH_dim, itype_real_max, itype_real_min, &
    & etype_periodic,etype_virtual, itype_periodic, itype_virtual, hsml_const
use particle_data ,   only: ntotal, etotal, ntotal_prev,etotal_prev, surf_norm, nedge_rel_edge, &
            & pBC_edges, pBC_epair_a, pBC_epair_s, pBC_eniac, tangent_pBC, edge_pBC_pairs, &
            & x,rho,mass, itype, hsml, edge, etype, pBC_duplicate_pair, nperiodic, simGridSize
    
    implicit none

!----------------------------------------------------------------------
!                   i,j: A variable used as a counter in loop
!                   d: A variable used to iterate dimensions
!
!----------------------------------------------------------------------

    integer(4) i,k,kn,ke,kpBC,a,s, s_pBC1, s_pBC2  
    real(8) ns_as,ts_as, scale_k, khsml, v1, v2
    real(8) ts_pBC1(SPH_dim), ts_pBC2(SPH_dim),dxr_as(SPH_dim), dxr_v1(SPH_dim), dxr_v2(SPH_dim)
    
   
    ! Find the smoothing kernel reach
    call sml_mult_factor(scale_k)
    khsml=scale_k*hsml_const
    
    ALLOCATE(pBC_duplicate_pair(2, ntotal))
    
    pBC_duplicate_pair=0
    
    kpBC=0
    nperiodic=0
    
    ! Store the total particle & edge number, in a variable that is updated for periodic particles
    kn=ntotal
    ke=etotal
    
    ! total particle & edge number before periodic bc is stored
    ntotal_prev=ntotal
    etotal_prev=etotal 
   
   ! Find all real particles that are in the reach of one of the periodic edge, and 
   ! copy them on the other side of the second periodic edge 
   do k=1,pBC_eniac

       a = pBC_epair_a(k)
       s_pBC1 = pBC_epair_s(k)
       if((itype(a) .le. itype_real_max) .and. (itype(a) .gt. itype_real_min)) then                                         
           ! Find which edge is being called as the particle pair       
           do i= 1, size(pBC_edges,2)
               ! For the identified edge, identify its periodic edge pair,
               ! and the sense of the vector
                if (s_pBC1 .eq. pBC_edges(1,i)) then
                    s_pBC2= pBC_edges(2,i)
                    ts_pBC1(:)=tangent_pBC(:,1,i)
                    ts_pBC2(:)=tangent_pBC(:,2,i)
                endif
            
                if (s_pBC1 .eq. pBC_edges(2,i)) then 
                    s_pBC2= pBC_edges(1,i)
                    ts_pBC1(:)=tangent_pBC(:,2,i)
                    ts_pBC2(:)=tangent_pBC(:,1,i)
                endif
            
           enddo
       
           !Find the distance vector from periodic edge mid point to particle 
           dxr_as(:)=x(:,a)-x(:,nedge_rel_edge(s_pBC1))
       
           !Split the distance vectors along the normal and tengential direction
           ! of the periodic edge
           ns_as= dot_product(surf_norm(:,s_pBC1),dxr_as(:))
           ts_as= dot_product(ts_pBC1(:),dxr_as(:))
       
       
           kn=kn+1
           !Determine the position of the periodic particle, using the distance vector
           x(:,kn)= -ns_as*surf_norm(:,s_pBC2)+ ts_as*ts_pBC2(:) + x(:,nedge_rel_edge(s_pBC2))
       
           !itype of copied particle is offset by a parameter for identification
           itype(kn)= mod(itype(a),itype_virtual) + itype_periodic 
           
           !copy hsml value to this particle
           hsml(kn)=hsml(a)
           mass(kn)=mass(a)
           rho(kn)=rho(a)
           
           ! The copied particle/point and the copy are paired for future reference
           kpBC=kpBC+1
           pBC_duplicate_pair(1,kpBC)=a
           pBC_duplicate_pair(2,kpBC)=kn
           
       endif
       
   enddo
   
  
! Find all edges that are in the reach of one of the periodic edge, and 
! copy them on the other side of the second periodic edge
! For this, first go over all the edges to find edges within the reach of either
! periodic edge
  do s=1,etotal
    
    !We remove the periodic edges themselves to avoid  copying them over their 
    ! periodic pair
    if( etype(s) .ne. etype_periodic) then
        
        ! identify the vertices associated to the edge
        v1=edge(1,s)
        v2=edge(2,s)
    
        ! Find if either of the vertices are in the reach of eithe periodic edge
        ! if they are then copy the associated edge of those vertices to the other
        ! side of the periodic edge pair
        do i= 1, size(pBC_edges,2)
            
            s_pBC1 = pBC_edges(1,i)
            s_pBC2 = pBC_edges(2,i)
            
            !Find the distance vector from first periodic edge mid point to edge points 
            dxr_v1(:)= x(:,v1) - x(:,nedge_rel_edge(s_pBC1))  
            dxr_v2(:)= x(:,v2) - x(:,nedge_rel_edge(s_pBC1))
            dxr_as(:)= x(:,nedge_rel_edge(s))  - x(:,nedge_rel_edge(s_pBC1))
        
            ! Check if the vertices are in reach for the first periodic edge pair 
            if( (dot_product(surf_norm(:,s_pBC1),dxr_v1(:)) .le. khsml)  &
                & .or. (dot_product(surf_norm(:,s_pBC1),dxr_v2(:)) .le. khsml) ) then 
            
                ke=ke+1
                ts_pBC1(:)=tangent_pBC(:,1,i)
                ts_pBC2(:)=tangent_pBC(:,2,i)
                
                ns_as= dot_product(surf_norm(:,s_pBC1),dxr_v1(:))
                ts_as= dot_product(ts_pBC1(:),dxr_v1(:))
                kn=kn+1
                !determine location of the copied vertice
                x(:,kn)=-ns_as*surf_norm(:,s_pBC2)+ ts_as*ts_pBC2(:) &
                        & +x(:,nedge_rel_edge(s_pBC2))
                edge(1,ke)=kn
                itype(kn)=0 !itype(edge(1,s))!+200
                hsml(kn)=hsml(v1)
                mass(kn)=mass(v1)
                rho(kn)=rho(v1)
                
                ! The copied particle/point and the copy are paired for future reference
                kpBC=kpBC+1
                pBC_duplicate_pair(1,kpBC)=v1
                pBC_duplicate_pair(2,kpBC)=kn

                ns_as= dot_product(surf_norm(:,s_pBC1),dxr_v2(:))
                ts_as= dot_product(ts_pBC1(:),dxr_v2(:))
                kn=kn+1
                !determine location of the copied vertice
                x(:,kn)=-ns_as*surf_norm(:,s_pBC2)+ ts_as*ts_pBC2(:) &
                        & +x(:,nedge_rel_edge(s_pBC2))              
                edge(2,ke)=kn
                itype(kn)=0 !itype(edge(2,s))+200
                hsml(kn)=hsml(v2)
                mass(kn)=mass(v2)
                rho(kn)=rho(v2)
                
                ! The copied particle/point and the copy are paired for future reference
                kpBC=kpBC+1
                pBC_duplicate_pair(1,kpBC)=v2
                pBC_duplicate_pair(2,kpBC)=kn
                
                ns_as= dot_product(surf_norm(:,s_pBC1),dxr_as(:))
                ts_as= dot_product(ts_pBC1(:),dxr_as(:))
                kn=kn+1
                !determine location of the copied edge reference point
                x(:,kn)=-ns_as*surf_norm(:,s_pBC2)+ ts_as*ts_pBC2(:) &
                        & +x(:,nedge_rel_edge(s_pBC2))
                itype(kn)= mod(itype(nedge_rel_edge(s)),itype_virtual)+ itype_periodic 
                hsml(kn)=hsml(nedge_rel_edge(s))
                mass(kn)=mass(nedge_rel_edge(s))
                rho(kn)=rho(nedge_rel_edge(s))
                
                kpBC=kpBC+1
                pBC_duplicate_pair(1,kpBC)=nedge_rel_edge(s)
                pBC_duplicate_pair(2,kpBC)=kn
                
                nedge_rel_edge(ke)=kn
                
                surf_norm(:,ke)=surf_norm(:,s)
                etype(ke)= mod(etype(s),etype_virtual) + etype_periodic
                
            endif
            
            !Find the distance vector from the second periodic edge mid point to edge points 
            dxr_v1(:)= x(:,v1) - x(:,nedge_rel_edge(s_pBC2))  
            dxr_v2(:)= x(:,v2) - x(:,nedge_rel_edge(s_pBC2))
            dxr_as(:)= x(:,nedge_rel_edge(s))  - x(:,nedge_rel_edge(s_pBC2))
            
            ! Check if the vertices are in reach for the second periodic edge pair          
            if( (dot_product(surf_norm(:,s_pBC2),dxr_v1(:)) .le. khsml)  &
                & .or. (dot_product(surf_norm(:,s_pBC2),dxr_v2(:)) .le. khsml) ) then 
            
                ke=ke+1
                ts_pBC2(:)=tangent_pBC(:,2,i)
                ts_pBC1(:)=tangent_pBC(:,1,i)
                
                ns_as= dot_product(surf_norm(:,s_pBC2),dxr_v1(:))
                ts_as= dot_product(ts_pBC2(:),dxr_v1(:))
                kn=kn+1
                !determine location of the copied vertice
                x(:,kn)=-ns_as*surf_norm(:,s_pBC1)+ ts_as*ts_pBC1(:) &
                        & +x(:,nedge_rel_edge(s_pBC1))              
                edge(1,ke)=kn
                itype(kn)=0 !itype(edge(1,s))+200
                hsml(kn)=hsml(v1)
                mass(kn)=mass(v1)
                rho(kn)=rho(v1)
                
                ! The copied particle/point and the copy are paired for future reference
                kpBC=kpBC+1
                pBC_duplicate_pair(1,kpBC)=v1
                pBC_duplicate_pair(2,kpBC)=kn
                
                ns_as= dot_product(surf_norm(:,s_pBC2),dxr_v2(:))
                ts_as= dot_product(ts_pBC2(:),dxr_v2(:))
                kn=kn+1
                !determine location of the copied vertice
                x(:,kn)=-ns_as*surf_norm(:,s_pBC1)+ ts_as*ts_pBC1(:) &
                        & +x(:,nedge_rel_edge(s_pBC1))                
                edge(2,ke)=kn
                itype(kn)=0! itype(edge(1,s))+200
                hsml(kn)=hsml(v2)
                mass(kn)=mass(v2)
                rho(kn)=rho(v2)
                
                ! The copied particle/point and the copy are paired for future reference
                kpBC=kpBC+1
                pBC_duplicate_pair(1,kpBC)=v2
                pBC_duplicate_pair(2,kpBC)=kn
                
                ns_as= dot_product(surf_norm(:,s_pBC2),dxr_as(:))
                ts_as= dot_product(ts_pBC2(:),dxr_as(:))
                kn=kn+1
                !determine location of the edge reference point
                x(:,kn)=-ns_as*surf_norm(:,s_pBC1)+ ts_as*ts_pBC1(:) &
                        & +x(:,nedge_rel_edge(s_pBC1))   
                itype(kn)=mod(itype(nedge_rel_edge(s)),itype_virtual)+ itype_periodic
                hsml(kn)=hsml(nedge_rel_edge(s))
                mass(kn)=mass(nedge_rel_edge(s))
                rho(kn)=rho(nedge_rel_edge(s))
                
                ! The copied particle/point and the copy are paired for future reference
                kpBC=kpBC+1
                pBC_duplicate_pair(1,kpBC)=nedge_rel_edge(s)
                pBC_duplicate_pair(2,kpBC)=kn
                
                nedge_rel_edge(ke)=kn
                
                surf_norm(:,ke)=surf_norm(:,s)
                etype(ke)=mod(etype(s),etype_virtual) + etype_periodic
            endif
            
            
        enddo
      
      
    endif
    
  enddo
  
  !update grid size:
  simGridSize(:,1)= MINVAL(x,2) 
  simGridSize(:,2)= MAXVAL(x,2) 
  
  !update the total particle and edge numbers
  ntotal=kn
  etotal=ke
  nperiodic=kpBC
  
  call  nnps_algorithm
  
  call  nnes_algorithm
  
 end
    