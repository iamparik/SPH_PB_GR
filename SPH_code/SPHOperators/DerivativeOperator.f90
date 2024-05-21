subroutine DerivativeOperator(dF_a,dF_b,F_a,F_b,Cdwdx_a, Cdwdx_b, mass_a, mass_b, rho_a, rho_b, id)
    !Calculates the derivative term from the particle-particle itneractionterm
    implicit none
    integer(4), intent(in) :: id
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
    endif

endsubroutine
    
    
    
subroutine DerivativeOperatorB(dF_a,F_a,F_s,Cdgmas, id)
    !Calculates the derivative term from the boundary integral part, due to aprticle edge itneraction
    implicit none
    integer(4), intent(in) :: id
    real(8), intent(in) :: F_a, F_s, Cdgmas
    real(8), intent(inout) :: dF_a  
    
    !if id =0 no operation is performed
    
    if(id .eq. 1) then
        
        dF_a= dF_a - F_s*Cdgmas
        
    elseif(id .eq. 2) then
        
        dF_a= dF_a - (F_s-F_a)*Cdgmas
        
    elseif(id .eq. 3) then
        
        dF_a= dF_a - (F_s+F_a)*Cdgmas
        
    endif

endsubroutine
    
    