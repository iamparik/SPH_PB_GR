!****************************************************************************
!
!  SUBROUTINE: ParticleShiftingTechnique
!
!  PURPOSE:  Subroutine to shift particles
!
!   CREATED:        03/25/2024       by  PARIKSHIT BOREGOWDA
!   Last Modified:  03/25/2024       by  PARIKSHIT BOREGOWDA 
!****************************************************************************

subroutine ParticleShiftingTechnique

use config_parameter, only:SPH_dim, itype_real_max, itype_real_min, &
        & PSTCoeff,PSTtype, hsml_const, dx_r, FScutoff
use particle_data, only: niac,pair_i, pair_j,eniac,epair_a, epair_s, &
    & w, dwdx, delC, surf_norm, del_gamma_as, nedge_rel_edge, ntotal, etotal, &
    & mass, rho, x, vx, itype, gamma_cont, gamma_discrt, gamma_mat_inv, FreeSurfaceVar


implicit none
real(8) maxShift, grad_b_term, dstress(SPH_dim), r_c, Pm1, Pm2
real(8) PSTshift, delr(SPH_dim), drDW(SPH_dim), dCF, w_dxr, dx_as, w_dxas
integer(4) k,d,a,b,s
real(8),DIMENSION(:,:),ALLOCATABLE:: xprev, dv1,dv2 , x_s
real(8),DIMENSION(:),ALLOCATABLE:: f0, fs0

grad_b_term= PSTCoeff*hsml_const**2.D0

r_c=0.D0
maxShift=0.2D0*dx_r
delr=0.D0
drDW=0.D0

Allocate(xprev(SPH_dim,ntotal))
xprev=x

Allocate(delC(SPH_dim,ntotal))    
delC=0.D0

call kernel(dx_r,delr,hsml_const,w_dxr,drDW)

do k = 1,niac
    a=pair_i(k)
    b=pair_j(k)
    
   
    ! By default DCF is 1
    dCF=1.D0 
    
    ! For different PSTtype we change 
    if (mod(PSTtype,10) .eq. 2) then
        dCF=1.D0 + 0.2D0*(w(k)/w_dxr)**4.D0
    endif
    
    
    do d=1,SPH_dim 
        call fncnApproxOperator(delC(d,a),dCF,mass(b),rho(b),dwdx(d,k)/gamma_cont(a))
        call fncnApproxOperator(delC(d,b),dCF,mass(a),rho(a),-dwdx(d,k)/gamma_cont(b))
    enddo

enddo


do k= 1, eniac
    a=epair_a(k)
    s=epair_s(k)
    b = nedge_rel_edge(s)
    
    dx_as=norm2(x(:,a)-x(:,b))
    call kernel(dx_as,delr,hsml_const,w_dxas,drDW)
    
    ! initalize r_c to 0 every loop
    r_c=0.D0
    
    ! Add repulisve field at solid boundary by using below
    if (PSTType .gt. 10) then
        Pm1=2.D0
        Pm2=2.D0
        if (dx_as .le. dx_r/2.D0) then
            !r_c=0.5D0-min(0.5D0,abs(dot_product(x(:,a)-x(:,b), surf_norm(:,s)))/dx_r)
            r_c=min(0.5D0,abs(dot_product(x(:,a)-x(:,b), surf_norm(:,s)))/dx_r)
            r_c= 0.5D0*(Pm2 - r_c*(Pm2-Pm1)/0.5D0 -1.D0)
        else
            r_c=0.D0
        endif
    endif    

    
    
    
    ! For different PSTtype we change r_c as follows 
    if (mod(PSTtype,10) .eq. 2) then
        r_c= -(1.D0 + 0.2D0*(w_dxas/w_dxr)**4.D0)-r_c
    else
    ! By default do below
        r_c= -1.D0-r_c
    endif
    
    
    
    do d=1,SPH_dim
        call fncnApproxOperator(delC(d,a),r_c,1.D0,1.D0,del_gamma_as(d,k)/gamma_cont(a))
    enddo
enddo  


do a=1,ntotal    
    dstress(:)=0.D0
    if((itype(a) .le. itype_real_max) .and. (itype(a) .gt. itype_real_min)) then 
        if( FreeSurfaceVar(a) .gt. FScutoff) then 
            dstress(:) = -grad_b_term*delC(:,a)
            PSTShift = min(norm2(dstress(:)), maxShift)
            if(PSTShift .gt. 1D-10*grad_b_term) then
                x(:,a) = x(:,a)+ PSTshift*(dstress(:)/norm2(dstress(:)))
            endif
        endif        
    endif  
    delc(:,a) =dstress(:)
enddo


Allocate(f0(ntotal),fs0(etotal))
f0(:)= vx(1,:)
do s= 1, etotal
    fs0(s)=vx(1,nedge_rel_edge(s))
enddo
Allocate(dv1(SPH_dim,ntotal))    
dv1=0.D0
call GradientPCBI(dv1,f0, fs0, gamma_mat_inv, SPH_dim, niac, pair_i, pair_j, dwdx, eniac, epair_a, epair_s, del_gamma_as, mass, rho, itype, ntotal, etotal)

f0(:)= vx(2,:)
do s= 1, etotal
    fs0(s)=vx(2,nedge_rel_edge(s))
enddo
Allocate(dv2(SPH_dim,ntotal))    
dv2=0.D0
call GradientPCBI(dv2,f0, fs0, gamma_mat_inv, SPH_dim, niac, pair_i, pair_j, dwdx, eniac, epair_a, epair_s, del_gamma_as, mass, rho, itype, ntotal, etotal)


do a=1,ntotal        
    if((itype(a) .le. itype_real_max) .and. (itype(a) .gt. itype_real_min)) then
        delr=x(:,a)-xprev(:,a)
        if( FreeSurfaceVar(a) .gt. 1.5) then
            vx(1,a)= vx(1,a) + dot_product(dv1(:,a),delr)
            vx(2,a)= vx(2,a) + dot_product(dv2(:,a),delr)
        endif                
    endif     
enddo


deallocate(xprev,dv1,dv2, f0,fs0)

end

    

