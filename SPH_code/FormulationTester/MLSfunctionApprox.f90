!
!  SUBROUTINE: MLSfunctionApprox
!			    This calculated boundary integral approximations of function at all points or particles:
!               ⟨〈f(x)〉⟩= ∑_b〖f(x_b ) MM_ab W_ab  m_b/ρ_b 〗
!
!  PURPOSE:  Subroutine to find MLS approximation
!
!   CREATED:        07/20/2023       by  PARIKSHIT BOREGOWDA
!   Last Modified:  07/20/2023       by  PARIKSHIT BOREGOWDA 
!****************************************************************************
    
subroutine MLSfunctionApprox( f, f0, b_MLS,x, niac, pair_i, pair_j, w, w_aa, mass, rho, itype, gamma_cont, ntotal)
    use config_parameter, only: itype_virtual, itype_real_max, itype_real_min, SPH_dim
    implicit none
    
    integer(4), intent(in):: ntotal, niac
    integer(2), intent(in):: itype(ntotal)
    real(8), intent(in):: f0(ntotal), b_MLS(SPH_dim+1,ntotal), mass(ntotal), rho(ntotal),w(niac) , w_aa(ntotal), x(SPH_dim,ntotal), gamma_cont(ntotal)
    integer(4), intent(in):: pair_i(niac), pair_j(niac) 
    real(8), intent(inout):: f(ntotal)
    integer(4)::k,a,b, d
    real(8)::MM, Px(SPH_dim+1)
    
    Px=1.D0
    f=0.D0
     do k = 1,niac
        a=pair_i(k)
        b=pair_j(k)
        
        ! Call the BI SPH operation for zeroth order derivative  of a function
        ! For ⟨〈f(x)〉⟩= 1/γ ∑_b〖f(x_b )W(‖x-x_b ‖,h)  m_b/ρ_b 〗
        Px(2:SPH_dim+1)= x(:,b)-x(:,a)
        !Px(2:SPH_dim+1)= x(:,a)-x(:,b)
        MM=0.D0
        do d=1,SPH_dim+1
            MM=MM+PX(d)*b_MLS(d,a)
        enddo        
        
        call fncnApproxOperator(f(a),f0(b),mass(b),rho(b),w(k)*MM)
        
        Px(2:SPH_dim+1)= x(:,a)-x(:,b)
        !Px(2:SPH_dim+1)= x(:,b)-x(:,a)
        MM=0.D0
        do d=1,SPH_dim+1
            MM=MM+PX(d)*b_MLS(d,b)
        enddo
        call fncnApproxOperator(f(b),f0(a),mass(a),rho(a),w(k)*MM)
        
     enddo
     
     do a=1,ntotal
        if ((mod(itype(a),itype_virtual) .le. itype_real_max) .and. (mod(itype(a),itype_virtual) .gt. itype_real_min)) then
            if(gamma_cont(a) .lt. 1) then
                MM=PX(1)*b_MLS(1,a)
                call fncnApproxOperator(f(a),f0(a),mass(a),rho(a),w_aa(a)*MM)
            else
                f(a)=f0(a)
            endif
            
        endif
    enddo
     
end
    
    