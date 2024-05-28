!
!  SUBROUTINE: zerothOrderFunction
!			    This calculated boundary integral approximations of function at all points or particles:
!               ⟨〈f(x)〉⟩= 1/γ_a ∑_b〖f(x_b ) W_ab  m_b/ρ_b 〗
!
!  PURPOSE:  Subroutine that calls various SPH operators 
!           to find interpolations of different functions
!
!   CREATED:        06/08/2022       by  PARIKSHIT BOREGOWDA
!   Last Modified:  12/06/2022       by  PARIKSHIT BOREGOWDA 
!****************************************************************************
    
subroutine zerothOrderFunction(f0, f, gma, niac, pair_i, pair_j, w, w_aa, mass, rho, itype, ntotal)
    use config_parameter, only: itype_virtual, itype_real_max, itype_real_min
    implicit none
    
    integer(4), intent(in):: ntotal, niac
    integer(2), intent(in):: itype(ntotal)
    real(8), intent(in):: f0(ntotal), gma(ntotal), mass(ntotal), rho(ntotal),w(niac) , w_aa(ntotal)
    integer(4), intent(in):: pair_i(niac), pair_j(niac) 
    real(8), intent(inout):: f(ntotal)
    integer(4)::k,a,b
    
    
     do k = 1,niac
        a=pair_i(k)
        b=pair_j(k)
        
        ! Call the BI SPH operation for zeroth order derivative  of a function
        ! For ⟨〈f(x)〉⟩= 1/γ ∑_b〖f(x_b )W(‖x-x_b ‖,h)  m_b/ρ_b 〗
        call fncnApproxOperator(f(a),f0(b),mass(b),rho(b),w(k)/gma(a))
        call fncnApproxOperator(f(b),f0(a),mass(a),rho(a),w(k)/gma(b))
        
     enddo
     
     do a=1,ntotal
        if ((mod(itype(a),itype_virtual) .le. itype_real_max) .and. (mod(itype(a),itype_virtual) .gt. itype_real_min)) &
            & call fncnApproxOperator(f(a),f0(a),mass(a),rho(a),w_aa(a)/gma(a))
    enddo
     
end
    
    