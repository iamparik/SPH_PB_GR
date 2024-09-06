subroutine summationDensityOperatorPtoP( rho_a, a, start_Setup, sum_dens_type)
    use config_parameter,   only: SPH_dim
    use particle_data, only: w, mass, rho, gamma_cont, gamma_discrt,niac, &
        & dgrho_prev
    
    implicit none
    
    real(8), intent(inout) :: rho_a
    integer(4), intent(in) :: a, sum_dens_type
    logical, intent(in) :: start_Setup
    real(8) :: rho_, Kbt, bt
        
    if(sum_dens_type .eq. 1) then
        rho_a= rho_a/gamma_cont(a)
        
    elseif(sum_dens_type .eq. 2) then
        rho_a= rho_a/gamma_discrt(a)
        
    elseif(sum_dens_type .eq. 3) then
        Kbt=30000.D0
        bt=exp(-Kbt*(min(gamma_discrt(a)/gamma_cont(a),1.D0)-1.D0)**2)  
        rho_a= rho_a/(gamma_cont(a)*bt+(1-bt)*gamma_discrt(a)) 
        
    elseif(sum_dens_type .eq. 4) then
        if(start_Setup) then
            Kbt=30000.D0
            bt=exp(-Kbt*(min(gamma_discrt(a)/gamma_cont(a),1.D0)-1.D0)**2)  
            rho_= rho_a/(gamma_cont(a)*bt+(1-bt)*gamma_discrt(a))   
            dgrho_prev(a)= rho_*gamma_cont(a)- rho_a !( [𝜌a *𝛾𝑎 ] - [∑𝑚𝑏*𝑊𝑎𝑏] )
            rho_a = rho_
        else
            rho_ = rho_a + dgrho_prev(a) !  [𝜌a *𝛾𝑎 ]_current =[∑𝑚𝑏*𝑊𝑎𝑏]_current + ( [𝜌a *𝛾𝑎 ]_prev - [∑𝑚𝑏*𝑊𝑎𝑏]_prev )
            dgrho_prev(a)= rho_ - rho_a
            rho_a= rho_/gamma_cont(a)    
        endif   
        
    elseif(sum_dens_type .eq. 5) then
        if(start_Setup) then
            rho_= rho(a) ! intial density setup   
            dgrho_prev(a)= rho_*gamma_cont(a)- rho_a !( [𝜌a *𝛾𝑎 ] - [∑𝑚𝑏*𝑊𝑎𝑏] )
            rho_a = rho_
        else
            
            rho_ = rho_a + dgrho_prev(a) !  [𝜌a *𝛾𝑎 ]_current =[∑𝑚𝑏*𝑊𝑎𝑏]_current + ( [𝜌a *𝛾𝑎 ]_prev - [∑𝑚𝑏*𝑊𝑎𝑏]_prev )
            dgrho_prev(a)= rho_ - rho_a
            rho_a= rho_/gamma_cont(a)    
        endif 
    endif

end subroutine
    
subroutine sumDens_ab(rho_a, rho_b, mass_a,mass_b,w_ab)    
    implicit none
    
    real(8), intent(inout) :: rho_a, rho_b
    real(8), intent(in) :: mass_a, mass_b, w_ab
    
    rho_a= rho_a + w_ab*mass_b
    rho_b= rho_b + w_ab*mass_a
end subroutine   
    
subroutine sumDens_aa(rho_a, mass_a, w_aa_)    
    implicit none
    
    real(8), intent(inout) :: rho_a
    real(8), intent(in) :: mass_a, w_aa_
    
    rho_a= rho_a + w_aa_*mass_a
end subroutine 
    
subroutine sumDens(rho_)
    use particle_data, only: w, w_aa, mass, pair_i, pair_j, ntotal, niac
    
    implicit none
    
    real(8), intent(inout):: rho_(ntotal)
    integer(4) :: a,b, k
    
    rho_=0.D0
    
    do k =1, niac
        a=pair_i(k)
        b=pair_j(k)        
        call sumDens_ab(rho_(a),rho_(b),mass(a), mass(b), w(k))        
    enddo
    
    do a= 1,ntotal
        call sumDens_aa(rho_(a), mass(a), w_aa(a))
    enddo
        
end subroutine
    
 