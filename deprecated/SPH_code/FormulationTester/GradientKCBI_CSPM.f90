!
!  SUBROUTINE: GradientKCBI_CSPM
!			    This calculates PCBI type approximations of function's first order derivative at all points or particles:
!               ⟨⟨∂i f(x)⟩⟩=[ξ_a^-1 ]_ij ∑_b〖-f_b ∂_j W_ab  m_b/ρ_b 〗+[ξ_a^-1 ]_ij ∑_s∫_(∂(Ω ∩ Ω_w )_s〖f_a W_(as' ) n_(s_j ) dS' 〗
!
!  PURPOSE:  Subroutine that calls various SPH operators 
!           to find interpolations of different functions
!
!   CREATED:        06/08/2022       by  PARIKSHIT BOREGOWDA
!   Last Modified:  12/06/2022       by  PARIKSHIT BOREGOWDA 
!****************************************************************************
    
subroutine GradientKCBI_CSPM( d1f, f0, Xi_inv, SPH_dim, niac, pair_i, pair_j, dwdx, eniac, epair_a, dgmas, mass, rho, itype, ntotal, num_edges)
    implicit none
    integer(4), intent(in):: ntotal, SPH_dim, niac, eniac, num_edges
    integer(2), intent(in):: itype(ntotal)
    real(8), intent(in):: f0(ntotal), Xi_inv(SPH_dim,SPH_dim,ntotal), mass(ntotal), rho(ntotal)
    real(8), intent(in):: dwdx(SPH_dim,niac), dgmas(SPH_dim,eniac)
    integer(4), intent(in)::pair_i(niac), pair_j(niac) 
    integer(4), intent(in):: epair_a(eniac)
    real(8), intent(inout):: d1f(SPH_dim,ntotal)
    integer(4)::k,a,b,d
    real(8) ::  XiInvtemp(SPH_dim,SPH_dim),temp_dwdx(SPH_dim),CdwdxA(SPH_dim),CdwdxB(SPH_dim)
    real(8) ::  temp_dgmas(SPH_dim),Cdgmas(SPH_dim)
          
    ! Call the interpolation operation for first order derivative of a function
     ! [ξ_a^-1 ]_ij ∑_b f_b ∂_j W_ab m_b/ρ_b 〗
     do k = 1,niac
        a=pair_i(k)
        b=pair_j(k)
        
        XiInvtemp(:,:) = Xi_inv(:,:,a)
        temp_dwdx(:)= dwdx(:,k)
        CdwdxA= MATMUL(XiInvtemp,temp_dwdx)
        
        XiInvtemp(:,:) = Xi_inv(:,:,b)
        temp_dwdx(:)= -dwdx(:,k)
        CdwdxB= MATMUL(XiInvtemp,temp_dwdx)

        do d=1,SPH_dim 

            call fncnApproxOperator(d1f(d,a),f0(b),mass(b),rho(b),CdwdxA(d))
            call fncnApproxOperator(d1f(d,b),f0(a),mass(a),rho(a),CdwdxB(d))

        enddo
        
     enddo
     
   
    ! Call the boundary itnegral terms of the first order derivative  of function interpolations
    ! [ξ_a^-1 ]_ij ∑_s∫f_a W_(as' ) n_(s_j ) dS'
     do k= 1, eniac
        a=epair_a(k)
        
        XiInvtemp(:,:) = Xi_inv(:,:,a)
        temp_dgmas(:)= dgmas(:,k)
        Cdgmas= MATMUL(XiInvtemp,temp_dgmas)

        do d=1,SPH_dim 
            
            call fncnApproxOperator(d1f(d,a),-f0(a),1.D0,1.D0,Cdgmas(d))
            
        enddo
   
    enddo
     
end
    
    