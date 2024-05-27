
!
!  SUBROUTINE: bdryPressureGradient
!			    This calculates PCBI type approximations of function's first order derivative at all points or particles:
!               ⟨⟨∂i f(x)⟩⟩=[Γ_a^-1 ]_ij ∑_b〖(f_a-f_b ) ∂_j W_ab  m_b/ρ_b 〗-[Γ_a^-1 ]_ij ∑_s∫_(∂(Ω ∩ Ω_w )_s〖(f_a-f_s' ) W_(as' ) n_(s_j ) dS' 〗
!
!  PURPOSE:  Subroutine that calls various SPH operators 
!           to find interpolations of different functions
!
!   CREATED:        10/09/2023       by  PARIKSHIT BOREGOWDA
!   Last Modified:  10/09/2023       by  PARIKSHIT BOREGOWDA 
!****************************************************************************
    
subroutine bdryPressureGradient(d1f,f0, Gma)
use config_parameter, only: SPH_dim, g_const,c_sound, prsrBdryType, dx_r, rho_init, hsml_const
    use particle_data, only: epair_a, epair_s, eniac, etotal,mass,rho, &
        & etype, del_gamma_as, nedge_rel_edge, ntotal, x, vx, surf_norm, gamma_density_cont
    implicit none
    real(8), intent(in):: f0(ntotal), Gma(ntotal)
    real(8), intent(inout):: d1f(SPH_dim,ntotal)
    integer(4)::k,a,b,s, d
    real(8) ::  delf, temp_dgmas(SPH_dim), Cdgmas(SPH_dim), g(SPH_dim), r_c, Pm1, Pm2
        
    ! Use gravity with correct direction
    g(:)=0.D0
    g(2)= -abs(g_const)

    ! Call the boundary itnegral terms of the first order derivative  of function interpolations
    ! [Γ_a^-1 ]_ij ∑_s∫_(∂(Ω ∩ Ω_w )_s〖(f_a-f_(s^' ) ) W_(as' ) n_(s_j ) dS' 〗
    do k= 1, eniac
        a=epair_a(k)
        s=epair_s(k)
        b = nedge_rel_edge(s)
        
        delf=0.D0
        
        temp_dgmas(:)= del_gamma_as(:,k)
        Cdgmas= temp_dgmas/Gma(a)
        
        if(prsrBdryType .eq. 3) then
            
             delf= - (f0(a) &
                     & - rho(a)*dot_product(g,(x(:,a)-x(:,b))) &
                     & + rho(a)*0.5D0*(dot_product(vx(:,a),vx(:,a))-dot_product(vx(:,b),vx(:,b))) &
                     & )
             
        elseif(prsrBdryType .eq. 4) then
            
             delf= - (f0(a) &
                     & - rho(a)*dot_product(g, (x(:,a)-x(:,b))) &
                     & + rho(a)*0.5D0*(dot_product(vx(:,a)-vx(:,b), surf_norm(:,s))**2) &
                     & )
             
        elseif(prsrBdryType .eq. 5) then
            
             delf= - (f0(a) &
                     & - rho(a)*c_sound*dot_product(vx(:,a)-vx(:,b), surf_norm(:,s)) &
                     & )
             
        elseif(prsrBdryType .eq. 6) then
             delf= - (f0(a) &
                     & - rho(a)*dot_product(g, (x(:,a)-x(:,b))) &
                     & - rho(a)*c_sound*dot_product(vx(:,a)-vx(:,b), surf_norm(:,s)) &
                     & )
                 
        elseif(prsrBdryType .eq. 7) then
            
            if (norm2(x(:,a)-x(:,b)) .le. dx_r) then
                r_c=0.5D0-min(dx_r/2.D0,abs(dot_product(x(:,a)-x(:,b), surf_norm(:,s))))/dx_r
            else
                r_c=0.D0
            endif
            
            delf= - (max(f0(a),f0(a)/(1.D0-r_c)) &
                     & - rho(a)*c_sound*dot_product(vx(:,a)-vx(:,b), surf_norm(:,s)) &
                     & )
                
        elseif(prsrBdryType .eq. 8) then
            r_c=0.D0
            Pm1=2.D0
            Pm2=2.D0
            if (norm2(x(:,a)-x(:,b)) .le. dx_r/2.D0) then
                r_c=min(0.5D0,abs(dot_product(x(:,a)-x(:,b), surf_norm(:,s)))/dx_r)
                r_c= (Pm2 - r_c*(Pm2-Pm1)/0.5D0)
            else
                r_c=1.D0
            endif
            
            delf= - (f0(a)*r_c &
                         & - rho(a)*c_sound*dot_product(vx(:,a)-vx(:,b), surf_norm(:,s)) &
                         & )
       
        elseif(prsrBdryType .eq. 9) then
            
            if(dot_product(vx(:,a), surf_norm(:,s)) .lt. 0.D0) then
                 delf= - (f0(a) &
                         & - rho(a)*c_sound*dot_product(vx(:,a)-vx(:,b), surf_norm(:,s)) &
                         & )
            else
                delf= - (f0(a))
            endif
            
         elseif(prsrBdryType .eq. 10) then
            

            if(abs(dot_product(x(:,a)-x(:,b), surf_norm(:,s))) .le. dx_r/2.D0)then
                delf= - (f0(a) &
                         & + rho(a)*(c_sound**2/hsml_const)*dot_product(x(:,a)-x(:,b), surf_norm(:,s)))
            else
                delf= - (f0(a) &
                         & - rho(a)*c_sound*dot_product(vx(:,a)-vx(:,b), surf_norm(:,s)) &
                         & )
            endif
            
        elseif(prsrBdryType .eq. 11) then
            
             delf= - (f0(a)/(gamma_density_cont(a)**2.D0) &
                     & - rho(a)*c_sound*dot_product(vx(:,a)-vx(:,b), surf_norm(:,s)) &
                     & )
            
        endif
             
        
         do d=1,SPH_dim           
            call fncnApproxOperator(d1f(d,a),delf,1.D0,1.D0,Cdgmas(d))
        enddo
    
    enddo
     
    end
    
    
    
    
        