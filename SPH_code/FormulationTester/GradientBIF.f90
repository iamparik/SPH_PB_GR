!
!  SUBROUTINE: ConcGradient
!			    This calculates KCBI approximations of function's first order derivative at all points or particles:
!               ⟨⟨∂i f(x)⟩⟩≅1/γ_a ∑_s∫_(∂(Ω ∩ Ω_w )_s)〖f(x_s' )W(‖x-x_s' ‖,h)n_(s_i )dS' 〗-1/γ_a ∑_b〖f(x_b )  ∂_(b_i ) W(‖x-x_b ‖,h)  m_b/ρ_b 〗
!
!  PURPOSE:  Subroutine that calls various SPH operators 
!           to find interpolations of different functions
!
!   CREATED:        06/08/2022       by  PARIKSHIT BOREGOWDA
!   Last Modified:  11/20/2023       by  PARIKSHIT BOREGOWDA 
!****************************************************************************
    
subroutine ConcGradient
    use config_parameter, only: SPH_dim
    use config_geometry, only: PSTGradtype,dx_r
    use particle_data, only: niac, pair_i, pair_j, dwdx,x, &
        & eniac, epair_a,epair_s, del_gamma_as, mass, rho, itype, ntotal, etotal, &
        & nedge_rel_edge, surf_norm, gamma_cont,delC
        
    implicit none    
    real(8), intent(inout):: d1f(SPH_dim,ntotal)
    integer(4)::k,a,b,s,d
    real(8) ::  delf, r_c
    
    if( .not. Allocated(delC)) Allocate(delC(SPH_dim,ntotal))
    
    delC=0.D0

    do k = 1,niac
        a=pair_i(k)
        b=pair_j(k)
        
        do d=1,SPH_dim 
            call fncnApproxOperator(delC(d,a),1.D0,mass(b),rho(b),dwdx(d,k)/gamma_cont(a))
            call fncnApproxOperator(delC(d,b),1.D0,mass(a),rho(a),-dwdx(d,k)/gamma_cont(a))
        enddo

    enddo
     
    do k= 1, eniac
        a=epair_a(k)
        s=epair_s(k)
        do d=1,SPH_dim
            call fncnApproxOperator(delC(d,a),-1.D0,1.D0,1.D0,del_gamma_as(d,k)/gamma_cont(a))
        enddo
    enddo   
end