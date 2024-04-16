!
!  SUBROUTINE: GradientKCBI
!			    This calculates KCBI approximations of function's first order derivative at all points or particles:
!               ⟨⟨∂i f(x)⟩⟩≅1/γ_a ∑_s∫_(∂(Ω ∩ Ω_w )_s)〖f(x_s' )W(‖x-x_s' ‖,h)n_(s_i )dS' 〗-1/γ_a ∑_b〖f(x_b )  ∂_(b_i ) W(‖x-x_b ‖,h)  m_b/ρ_b 〗
!
!  PURPOSE:  Subroutine that calls various SPH operators 
!           to find interpolations of different functions
!
!   CREATED:        06/08/2022       by  PARIKSHIT BOREGOWDA
!   Last Modified:  12/04/2022       by  PARIKSHIT BOREGOWDA 
!****************************************************************************
    
subroutine GradientKCBI(d1f, f0, fs0, gma, SPH_dim, niac, pair_i, pair_j, dwdx, eniac, epair_a, epair_s, dgmas, mass, rho, itype, ntotal, num_edges)

    implicit none    
    integer(4), intent(in):: ntotal, SPH_dim, niac, eniac, num_edges
    integer(2), intent(in):: itype(ntotal)
    real(8), intent(in):: f0(ntotal), gma(ntotal), mass(ntotal), rho(ntotal), fs0(num_edges) 
    real(8), intent(in):: dwdx(SPH_dim,niac), dgmas(SPH_dim,eniac) 
    integer(4), intent(in):: pair_i(niac), pair_j(niac)
    integer(4), intent(in):: epair_a(eniac), epair_s(eniac)
    real(8), intent(inout):: d1f(SPH_dim,ntotal)
    integer(4)::k,a,b,s,d
    real(8) ::  delf
    
          
    ! Call the interpolation operation for first order derivative of a function
     ! 1/γ ∑_b〖f(x_b )  ∂_(b_i ) W(‖x-x_b ‖,h)  m_b/ρ_b 〗
     do k = 1,niac
        a=pair_i(k)
        b=pair_j(k)
        
         do d=1,SPH_dim 
            call fncnApproxOperator(d1f(d,a),f0(b),mass(b),rho(b),dwdx(d,k)/gma(a))
            call fncnApproxOperator(d1f(d,b),f0(a),mass(a),rho(a),-dwdx(d,k)/gma(b))
        enddo

     enddo
     
   
    ! Call the boundary itnegral terms of the first order derivative  of function interpolations
    ! 1/γ ∑_s∫_(∂(Ω ∩ Ω_w )_s)〖f(x_s' )W(‖x-x_s'‖,h)n_(s_i)dS' 〗
     do k= 1, eniac
        a=epair_a(k)
        s=epair_s(k)

        do d=1,SPH_dim
            call fncnApproxOperator(d1f(d,a),-fs0(s),1.D0,1.D0,dgmas(d,k)/gma(a))
        enddo
   
    enddo
     
end
    
    