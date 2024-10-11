subroutine VectorDivergencePtoP(divF_a,divF_b,a,b,F_a,F_b,Cdwdx_a, Cdwdx_b, mass_a, mass_b, rho_a, rho_b, dim , divType)
    ! Divergence term due to particle particle interactions
    implicit None
    integer(4), intent(in) :: dim, divType,a,b
    real(8), intent(in) :: F_a(dim), F_b(dim), Cdwdx_a(dim), Cdwdx_b(dim), mass_a, mass_b, rho_a, rho_b
    real(8), intent(inout) :: divF_a,divF_b
    integer(4) :: d
    real(8) :: delf 
    
    do d= 1,dim
        call DerivativeOperator(divF_a,divF_b,a,b,F_a(d),F_b(d),Cdwdx_a(d), Cdwdx_b(d), mass_a, mass_b, rho_a, rho_b, divType)
    enddo

end 

subroutine VectorDivergencePtoB(divF_a,a,s,F_a,F_s,Cdgmas,dim, divType)
    ! Divergence term due to particle edge interactions
    implicit None
    integer(4), intent(in) :: dim, divType,a,s
    real(8), intent(in) :: F_a(dim), F_s(dim), Cdgmas(dim)
    real(8), intent(inout) :: divF_a
    integer(4) :: d

    do d= 1,dim
        call DerivativeOperatorB(divF_a,a,s,F_a(d),F_s(d),Cdgmas(d), divType)
    enddo

end 
