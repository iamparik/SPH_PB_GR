!
!  SUBROUTINE: GradientFerrand
!			    This calculates PCBI type approximations of function's first order derivative at all points or particles:
!               ⟨⟨∂i f(x)⟩⟩=[Γ_a^-1 ]_ij ∑_b〖(f_a-f_b ) ∂_j W_ab  m_b/ρ_b 〗-[Γ_a^-1 ]_ij ∑_s∫_(∂(Ω ∩ Ω_w )_s〖(f_a-f_s' ) W_(as' ) n_(s_j ) dS' 〗
!
!  PURPOSE:  Subroutine that calls various SPH operators 
!           to find interpolations of different functions
!
!   CREATED:        06/08/2022       by  PARIKSHIT BOREGOWDA
!   Last Modified:  01/25/2023       by  PARIKSHIT BOREGOWDA 
!****************************************************************************
    
subroutine GradientFerrand(d1f,f0, fs0, Gma, SPH_dim, niac, pair_i, pair_j, dwdx, eniac, epair_a, epair_s, dgmas, mass, rho, rhos, itype, ntotal, num_edges)
    implicit none
    integer(4), intent(in):: ntotal, SPH_dim, niac, eniac, num_edges
    integer(2), intent(in):: itype(ntotal)
    real(8), intent(in):: f0(ntotal), Gma(ntotal), mass(ntotal), rho(ntotal), rhos(num_edges), fs0(num_edges)
    real(8), intent(in):: dwdx(SPH_dim,niac), dgmas(SPH_dim,eniac) 
    integer(4), intent(in):: pair_i(niac), pair_j(niac) 
    integer(4), intent(in):: epair_a(eniac), epair_s(eniac)
    real(8), intent(inout):: d1f(SPH_dim,ntotal)
    integer(4)::k,a,b,s, d
    real(8) ::  delf,temp_dwdx(SPH_dim),CdwdxA(SPH_dim),CdwdxB(SPH_dim)
    real(8) ::  temp_dgmas(SPH_dim),Cdgmas(SPH_dim)
    
          
    ! Call the interpolation operation for first order derivative of a function
    ! (ρ_a/𝛾_𝑎)∑_b〖(f_a/ρ_b**2 + f_b/ρ_b**2 ) ∂_j W_ab m_b 〗
    
     do k = 1,niac
        a=pair_i(k)
        b=pair_j(k)

        temp_dwdx(:)= dwdx(:,k)
        CdwdxA= temp_dwdx/Gma(a)
        
        temp_dwdx(:)= -dwdx(:,k)
        CdwdxB= temp_dwdx/Gma(b)
        
        do d= 1,SPH_dim
            delf= rho(a)*(f0(b)/(rho(b)**2)+f0(a)/(rho(a)**2))
            call fncnApproxOperator(d1f(d,a),delf,mass(b),1.D0,CdwdxA(d))
            delf= rho(b)*(f0(b)/(rho(b)**2)+f0(a)/(rho(a)**2))
            call fncnApproxOperator(d1f(d,b),delf,mass(a),1.D0,CdwdxB(d))
        enddo
        
        
     enddo
     
   
    ! Call the boundary itnegral terms of the first order derivative  of function interpolations
    ! [Γ_a^-1 ]_ij ∑_s∫_(∂(Ω ∩ Ω_w )_s〖(f_a-f_(s^' ) ) W_(as' ) n_(s_j ) dS' 〗
    do k= 1, eniac
        a=epair_a(k)
        s=epair_s(k)
        
        temp_dgmas(:)= dgmas(:,k)
        Cdgmas= temp_dgmas/Gma(a)
        
         do d=1,SPH_dim
             delf= rhos(s)*rho(a)*(-f0(a)/(rho(a)**2)-fs0(s)/(rhos(s)**2))
             !delf= rho(a)*rho(a)*(-f0(a)/(rho(a)**2)-f0(a)/(rho(a)**2))
            call fncnApproxOperator(d1f(d,a),delf,1.D0,1.D0,Cdgmas(d))
        enddo
    
    enddo
     
    end
    
    
    
    
    