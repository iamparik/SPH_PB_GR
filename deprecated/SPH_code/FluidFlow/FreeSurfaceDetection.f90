subroutine FreeSurfaceDetection
    
use config_parameter, only:SPH_dim, itype_real_max, itype_real_min
use particle_data, only: niac,pair_i, pair_j,eniac,epair_a, epair_s, &
    & dwdx, delC, surf_norm, del_gamma_as, nedge_rel_edge, ntotal, etotal, &
    & mass, rho, x, vx, itype, gamma_cont, gamma_discrt,FreeSurfaceVar

implicit none

integer(4) :: s,a ,d
real(8),DIMENSION(:,:),ALLOCATABLE:: x_s
real(8),DIMENSION(:,:,:),ALLOCATABLE:: gamma_mat_approx_inv

ALLOCATE(x_s(SPH_dim,etotal), FreeSurfaceVar(ntotal))
x_s=0.D0
FreeSurfaceVar=0.D0

! This evaluates function value at the boundary. This can be further generalized
do s= 1, etotal
    x_s(:,s)=x(:,nedge_rel_edge(s))
enddo

!call DivergenceKCBI( FreeSurfaceVar, x, x_s, gamma_cont, SPH_dim, niac, pair_i, pair_j, dwdx, eniac, epair_a, epair_s, del_gamma_as, mass, rho, itype, ntotal, etotal)


 ALLOCATE( gamma_mat_approx_inv(SPH_dim,SPH_dim,ntotal)) 
    gamma_mat_approx_inv=0.D0
    do a=1,ntotal
        do d=1,SPH_dim
                gamma_mat_approx_inv(d,d,a)= 1.D0/gamma_cont(a)
        enddo
    enddo
        
    ! Determine function's first order derivative approximation using WCBI1 formulation
call DivergencePCBI( FreeSurfaceVar, x, x_s, gamma_mat_approx_inv, SPH_dim, niac, pair_i, pair_j, dwdx, eniac, epair_a, epair_s, del_gamma_as, mass, rho, itype, ntotal, etotal)

DEALLOCATE( gamma_mat_approx_inv )

deallocate(x_s)




end
    

