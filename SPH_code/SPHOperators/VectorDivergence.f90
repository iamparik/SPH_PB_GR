subroutine VectorDivergencePtoP(divF_a,divF_b,F_a,F_b,Cdwdx_a, Cdwdx_b, mass_a, mass_b, rho_a, rho_b, dim )
    implicit None
    integer(4), intent(in) :: dim
    real(8), intent(in) :: F_a(dim), F_b(dim), Cdwdx_a(dim), Cdwdx_b(dim), mass_a, mass_b, rho_a, rho_b
    real(8), intent(out) :: divF_a,divF_b
    integer(4) :: d
    real(8) :: delf 
    do d= 1,dim
        delf= F_b(d)-F_a(d)
        call fncnApproxOperator(divF_a,delf,mass_b,rho_b,Cdwdx_a(d))
        delf= F_a(d)-F_b(d)
        call fncnApproxOperator(divF_b,delf,mass_a,rho_a,Cdwdx_b(d))
    enddo

end 
    




subroutine VectorDivergencePtoB(divF_a,F_a,F_b,Cdgmas,dim)
    implicit None
    integer(4), intent(in) :: dim
    real(8), intent(in) :: F_a(dim), F_b(dim), Cdgmas(dim)
    real(8), intent(out) :: divF_a
    integer(4) :: d
    real(8) :: delf 

    do d= 1,dim
        delf= F_a(d)-F_b(d)
        call fncnApproxOperator(divF_a,delf,1.D0,1.D0,Cdgmas(d))
    enddo

end 
