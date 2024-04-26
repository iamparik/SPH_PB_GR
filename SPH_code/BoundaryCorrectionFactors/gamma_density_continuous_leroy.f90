!****************************************************************************
!
!  SUBROUTINE: gamma_continuous_Ferrand
!
!  PURPOSE: This defines the continuous wall renormalization factor gamma 
!           in Leroy et. al paper 2014 JCP
!
!   CREATED:        11/13/2022       by  PARIKSHIT BOREGOWDA
!   Last Modified:  04/01/2024        by  PARIKSHIT BOREGOWDA 
!****************************************************************************
subroutine gamma_density_continuous_leroy
    
use config_parameter, only: SPH_dim,pi, dx_r
use particle_data, only: x,gamma_density_cont, eniac, epair_a, epair_s, &
                        & edge,hsml, surf_norm, maxn 

implicit none

integer(4) k,a,s,v
real(8)  xs_v(SPH_dim,size(edge,1)), dist_as
real(8), dimension(:), allocatable:: gamma_as


! gamma_as is initialized to 0.D0
if( .NOT. allocated(gamma_as)) allocate(gamma_as(eniac))
gamma_as=0.D0

! gamma_cont is initialized to 1
if( .NOT. allocated(gamma_density_cont)) allocate(gamma_density_cont(maxn)) 
gamma_density_cont=1.D0


! This loop calculates   𝛾_𝑎=1−∑1_𝑠▒𝛾_𝑎𝑠 
do k=1,eniac
    a=epair_a(k)
    s=epair_s(k)
    
    if( SPH_dim .eq. 2) then
    
        do v=1,size(edge,1)
            xs_v(:,v) = x(:,edge(v,s))
        enddo
        
        dist_as= min(abs(dot_product(surf_norm(:,s), (xs_v(:,1)-x(:,a)))),norm2(xs_v(:,1)-x(:,a))) 
        dist_as= min(dist_as,norm2(xs_v(:,2)-x(:,a))) 
        if( dist_as .lt. dx_r/2.D0) then
            call gamma_analytical_leroy_2D(x(:,a),xs_v,surf_norm(:,s),gamma_as(k),dx_r/2.D0,pi)
        else
        gamma_as(k)=0.D0
        endif
    
        
    else
        
    endif
    
    gamma_density_cont(a)= gamma_density_cont(a) - gamma_as(k)
    
enddo

deallocate(gamma_as)
    
end    