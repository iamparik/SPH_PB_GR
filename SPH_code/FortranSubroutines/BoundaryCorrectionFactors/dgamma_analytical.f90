!****************************************************************************
!
!  SUBROUTINE: dgamma_analytical
!
!  PURPOSE: This defines the gradiant of gamma, and is calculated analytically
!           using boundary itnegration of the Kernel. 
!
!   CREATED:        04/13/2022       by  PARIKSHIT BOREGOWDA
!   Last Modified:  04/18/2022        by  PARIKSHIT BOREGOWDA 
!****************************************************************************
subroutine dgamma_analytical
use config_parameter, only:  pi
use particle_data, only: SPH_dim,del_gamma_as, del_gamma, edge, &
            &  surf_norm, eniac, epair_a, epair_s, x, hsml, ntotal

implicit none

integer(4) a,s,k,v 
real(8) xs_v(SPH_dim,size(edge,1)), xa(SPH_dim), scale_k

if( .NOT. allocated(del_gamma_as)) allocate(del_gamma_as(SPH_dim,eniac)) 
if( .NOT. allocated(del_gamma)) allocate(del_gamma(SPH_dim,ntotal)) 

del_gamma_as=0.D0
del_gamma=0.D0

call sml_mult_factor(scale_k)

do k= 1, eniac
    a=epair_a(k)
    s=epair_s(k)
    
    xa(:)=x(:,a)

    
    if (SPH_dim .eq. 2) then
        
        ! The below loop extracts postion of all vertices
        ! assosciated to an edge. This will need to be updated, if different 
        ! types of boundary elements are present (like triangles and rectangles)
        do v=1,size(edge,1)
            xs_v(:,v)=x(:,edge(v,s))    
        enddo
        
        call dgamma_analytical_Ferrand_2D(xa,xs_v,surf_norm(:,s), del_gamma_as(:,k),hsml(a), pi, scale_k) 
        !note scale_k currently commented in the above subroutine
        
    else
        write(*,*) ' >>> ERROR <<< : analytical dgamma calcuation not defined for dimension', SPH_dim
        pause    
    endif
    
    del_gamma(:,a)= del_gamma(:,a)+ del_gamma_as(:,k)

enddo

end


