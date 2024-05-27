subroutine densDiffOperatorPtoP(dF_a,dF_b,F_a,F_b,dwdx, mass_a, mass_b, rho_a, rho_b, &
                    & dim, hsml, c_sound, x_ab, delta_SPH, CF_ID)
    implicit none

    integer(4), intent(in) :: dim, CF_ID
    real(8), intent(in) :: x_ab(dim),F_a(dim), F_b(dim), dwdx(dim), &
        & mass_a, mass_b, rho_a, rho_b, hsml, c_sound, delta_SPH
    real(8), intent(inout) :: dF_a,dF_b
    integer(4) :: d 
    real(8) :: Cdwdx_a(dim), Cdwdx_b(dim), eta, diff_term(dim)
    
    Cdwdx_a(:)=dwdx(:)
    Cdwdx_b(:)=-dwdx(:)
    
    if(CF_ID .gt. 0) then
    !---------------------Monoghan artificial viscosity term -------------------------
        if((mass_a .gt. 0.D0) .and. (mass_b .gt. 0.D0)) then
            eta=0.1*hsml
        
            diff_term=2.D0*(rho_a-rho_b)*(x_ab)/(norm2(x_ab)**2)- (F_a+F_b)
            dF_a= dF_a + delta_SPH*hsml*c_sound*dot_product(diff_term,Cdwdx_a)*mass_b/rho_b
            dF_b= dF_b + delta_SPH*hsml*c_sound*dot_product(diff_term,Cdwdx_b)*mass_a/rho_a
            
            if(isNAN(dF_a) .or. isNAN(dF_b)) pause
        endif
    
    endif
    
endsubroutine