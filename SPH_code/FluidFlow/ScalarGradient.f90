subroutine ScalarGradientPtoP(delF_a,delF_b,F_a,F_b,Cdwdx_a, Cdwdx_b, mass_a, mass_b, rho_a, rho_b, dim , divType)
    ! Gradient term due to particle particle interactions
    implicit None
    integer(4), intent(in) :: dim, divType
    real(8), intent(in) :: F_a, F_b, Cdwdx_a(dim), Cdwdx_b(dim), mass_a, mass_b, rho_a, rho_b
    real(8), intent(inout) :: delF_a(dim),delF_b(dim)
    integer(4) :: d
    real(8) :: delf 
    
    do d= 1,dim
        call DerivativeOperator(delF_a(d),delF_b(d),F_a,F_b,Cdwdx_a(d), Cdwdx_b(d), mass_a, mass_b, rho_a, rho_b, divType)
    enddo

end 

subroutine ScalarGradientPtoB(delF_a,F_a,F_s,Cdgmas,dim, divType)
    ! Gradient term due to particle edge interactions
    implicit None
    integer(4), intent(in) :: dim, divType
    real(8), intent(in) :: F_a, F_s, Cdgmas(dim)
    real(8), intent(inout) :: delF_a(dim)
    integer(4) :: d

    do d= 1,dim
        call DerivativeOperatorB(delF_a(d),F_a,F_s,Cdgmas(d), divType)
    enddo

end 