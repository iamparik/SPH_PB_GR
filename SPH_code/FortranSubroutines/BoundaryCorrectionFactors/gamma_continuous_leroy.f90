!****************************************************************************
!
!  SUBROUTINE: gamma_continuous_Ferrand
!
!  PURPOSE: This defines the continuous wall renormalization factor gamma 
!           in Leroy et. al paper 2014 JCP
!
!   CREATED:        11/13/2022       by  PARIKSHIT BOREGOWDA
!   Last Modified:  11/13/2022        by  PARIKSHIT BOREGOWDA 
!****************************************************************************
subroutine gamma_continuous_leroy
    
use config_parameter, only: pi
use particle_data, only: SPH_dim,x,gamma_cont, eniac, epair_a, epair_s, &
                        & edge,hsml, surf_norm, maxn 

implicit none

integer(4) k,a,s,v
real(8)  xs_v(SPH_dim,size(edge,1)) 
real(8), dimension(:), allocatable:: gamma_as


! gamma_as is initialized to 0.D0
if( .NOT. allocated(gamma_as)) allocate(gamma_as(eniac))
gamma_as=0.D0

! gamma_cont is initialized to 1
if( .NOT. allocated(gamma_cont)) allocate(gamma_cont(maxn)) 
gamma_cont=1.D0


! This loop calculates   𝛾_𝑎=1−∑1_𝑠▒𝛾_𝑎𝑠 
do k=1,eniac
    a=epair_a(k)
    s=epair_s(k)
    
    if( SPH_dim .eq. 2) then
    
        do v=1,size(edge,1)
            xs_v(:,v) = x(:,edge(v,s))
        enddo
        
        call gamma_analytical_leroy_2D(x(:,a),xs_v,surf_norm(:,s),gamma_as(k),hsml(a),pi)
        
        
    else
        
    endif
    
    gamma_cont(a)= gamma_cont(a) - gamma_as(k)
    
enddo

    
    
end    