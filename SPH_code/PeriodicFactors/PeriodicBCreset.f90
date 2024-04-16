!***************************************************************************
!
!  SUBROUTINE: PeriodicBCreset
!
!  PURPOSE: For simulations having periodic boundary condition
!           the periodic particles are removed from numbering
!           and variables assosciated with pBC reset 
!
!   CREATED:        08/12/2022       by  PARIKSHIT BOREGOWDA
!   Last Modified:  08/12/2022       by  PARIKSHIT BOREGOWDA 
!****************************************************************************
subroutine PeriodicBCreset
    use config_parameter, only:SPH_dim, itype_real_max, itype_real_min
    use particle_data, only:ntotal, etotal, ntotal_prev,etotal_prev, &
        & pBC_duplicate_pair,surf_norm, nedge_rel_edge, &
        & pBC_edges, pBC_epair_a, pBC_epair_s, pBC_eniac, tangent_pBC, &
        & x, itype
    
    implicit none 
    
    integer(4) i,k,a,s, s_pBC1, s_pBC2  
    real(8) ts_pBC1(SPH_dim), ts_pBC2(SPH_dim),dxr_as(SPH_dim),ns_as,ts_as
    
    
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
           
           if (ns_as .lt. 0.D0) then                               
               x(:, a)= -ns_as*surf_norm(:,s_pBC2)+ ts_as*ts_pBC2(:) + x(:,nedge_rel_edge(s_pBC2))
           endif                                     

       endif
       
   enddo               
    
    
    deallocate(pBC_duplicate_pair)
    
    ntotal=ntotal_prev
    etotal=etotal_prev
    
    end
    
    
