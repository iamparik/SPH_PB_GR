subroutine artViscOperatorPtoP(dF_a,dF_b,F_a,F_b,dwdx, mass_a, mass_b, rho_a, rho_b, &
                    & dim, hsml, c_sound, x_ab,CF_ID)
    implicit none

    integer(4), intent(in) :: dim, CF_ID
    real(8), intent(in) :: x_ab(dim),F_a(dim), F_b(dim), dwdx(dim), &
        & mass_a, mass_b, rho_a, rho_b, hsml, c_sound
    real(8), intent(inout) :: dF_a(dim),dF_b(dim)
    integer(4) :: d, Scalar0Matrix1 
    real(8) :: Cdwdx_a(dim), Cdwdx_b(dim), eta, alpha, beta, MU_ab, PI_ab
    
     
    Cdwdx_a(:)=dwdx(:)
    Cdwdx_b(:)=-dwdx(:)
    
    if(CF_ID .eq. 1) then
    !---------------------Monoghan artificial viscosity term -------------------------
        if((mass_a .gt. 0.D0) .and. (mass_b .gt. 0.D0)) then
            eta=0.1*hsml
        
            alpha = 1.D0
            beta = 2.D0
        
            MU_ab= hsml*dot_product(F_a-F_b,x_ab)/(norm2(x_ab)**2+eta**2)
            
            if(MU_ab .lt. 0.D0) then
        
                PI_ab=(-alpha*c_sound*MU_ab+ beta* MU_ab**2)/((rho_a+rho_b)/2.D0)
        
                do d=1,dim
                    dF_a(d)= dF_a(d) + PI_ab*Cdwdx_a(d)*mass_b*rho_a
                    dF_b(d)= dF_b(d) + PI_ab*Cdwdx_b(d)*mass_a*rho_b
                enddo
            endif
        endif
    
    endif
    
endsubroutine
    