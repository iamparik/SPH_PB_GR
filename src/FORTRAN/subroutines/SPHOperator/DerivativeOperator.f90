subroutine DerivativeOperator(dF_a,dF_b,a,b,F_a,F_b,Cdwdx_a, Cdwdx_b, mass_a, mass_b, rho_a, rho_b, id)
    !Calculates the derivative term from the particle-particle itneractionterm
    use particle_data, only : gamma_cont
    implicit none
    integer(4), intent(in) :: id,a,b
    real(8), intent(in) :: F_a, F_b, Cdwdx_a, Cdwdx_b, mass_a, mass_b, rho_a, rho_b
    real(8), intent(inout) :: dF_a,dF_b    
    
    !if id =0 no operation is performed
    
    if(id .eq. 1) then
        
        dF_a= dF_a + F_b*Cdwdx_a*mass_b/rho_b
        dF_b= dF_b + F_a*Cdwdx_b*mass_a/rho_a
        
    elseif(id .eq. 2) then
        
        dF_a= dF_a + (F_b-F_a)*Cdwdx_a*mass_b/rho_b
        dF_b= dF_b + (F_a-F_b)*Cdwdx_b*mass_a/rho_a
        
    elseif(id .eq. 3) then
        
        dF_a= dF_a + (F_b+F_a)*Cdwdx_a*mass_b/rho_b
        dF_b= dF_b + (F_a+F_b)*Cdwdx_b*mass_a/rho_a
    
   elseif(id .eq. 4) then
        
        dF_a= dF_a + rho_a*rho_b*(F_b/rho_b**2+F_a/rho_a**2)*Cdwdx_a*mass_b/rho_b
        dF_b= dF_b + rho_a*rho_b*(F_a/rho_a**2+F_b/rho_b**2)*Cdwdx_b*mass_a/rho_a   
        
    elseif(id .eq. 5) then
        
        dF_a= dF_a + rho_a*rho_b*(F_b/(gamma_cont(b)*rho_b**2)+F_a/(gamma_cont(a)*rho_a**2))*Cdwdx_a*mass_b/rho_b
        dF_b= dF_b + rho_a*rho_b*(F_a/(gamma_cont(a)*rho_a**2)+F_b/(gamma_cont(b)*rho_b**2))*Cdwdx_b*mass_a/rho_a     
    
    endif

endsubroutine
    
    
    
subroutine DerivativeOperatorB(dF_a,a,s,F_a,F_s,Cdgmas, id)
    !Calculates the derivative term from the boundary integral part, due to aprticle edge itneraction
    use particle_data, only : rho, rho_s, gamma_cont
    implicit none
    integer(4), intent(in) :: id,a,s
    real(8), intent(in) :: F_a, F_s, Cdgmas
    real(8), intent(inout) :: dF_a  
    
    !if id =0 no operation is performed
    
    if(id .eq. 1) then
        
        dF_a= dF_a - F_s*Cdgmas
        
    elseif(id .eq. 2) then
        
        dF_a= dF_a - (F_s-F_a)*Cdgmas
        
    elseif(id .eq. 3) then
        
        dF_a= dF_a - (F_s+F_a)*Cdgmas
    
    elseif(id .eq. 4) then
        if(.not. allocated(rho_s)) write(*,*) "rho_s not deallocated, use appropriate pressure bdry type such that prsr_bdry_preCalc=.true."
        dF_a= dF_a - rho_s(s)*rho(a)*(F_s/(rho_s(s)**2)+F_a/(rho(a)**2))*Cdgmas
    
    elseif(id .eq. 5) then
        
        dF_a= dF_a - F_a*Cdgmas/gamma_cont(a)
    endif

endsubroutine
    
    