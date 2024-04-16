!
!  SUBROUTINE: CSPM_gradient
!			    This calculates CSPM type approximations of function's first order derivative at all points or particles:
!               ⟨⟨∂i f(x)⟩⟩=[ξ_a^-1 ]_ij ∑_b〖(f_a-f_b ) ∂_j W_ab  m_b/ρ_b 〗
!
!  PURPOSE:  Subroutine that calls various SPH operators 
!           to find interpolations of different functions
!
!   CREATED:        06/09/2022       by  PARIKSHIT BOREGOWDA
!   Last Modified:  12/06/2022       by  PARIKSHIT BOREGOWDA 
!****************************************************************************
    
subroutine GradientCSPM(d1f,f0, Xi_inv, SPH_dim, niac, pair_i, pair_j, dwdx, mass, rho, itype, ntotal)
    use config_parameter, only: itype_real_min, itype_real_max, itype_virtual
    implicit none
    integer(4), intent(in):: ntotal, SPH_dim, niac 
    integer(2), intent(in):: itype(ntotal)
    real(8), intent(in):: f0(ntotal), Xi_inv(SPH_dim,SPH_dim,ntotal), mass(ntotal), rho(ntotal)
    real(8), intent(in):: dwdx(SPH_dim,niac)
    integer(4), intent(in):: pair_i(niac), pair_j(niac) 
    real(8), intent(inout):: d1f(SPH_dim,ntotal)
    integer(4)::k,a,b, d
    real(8) ::  delf,XiInvtemp(SPH_dim,SPH_dim),temp_dwdx(SPH_dim),CdwdxA(SPH_dim),CdwdxB(SPH_dim)
    
          
    ! Call the interpolation operation for first order derivative of a function
    ! [ξ_a^-1 ]_ij ∑_b〖(f_a-f_b ) ∂_j W_ab m_b/ρ_b 〗
    
     do k = 1,niac
        a=pair_i(k)
        b=pair_j(k)
        
        XiInvtemp(:,:) = Xi_inv(:,:,a)
        temp_dwdx(:)= dwdx(:,k)
        CdwdxA= MATMUL(XiInvtemp,temp_dwdx)
        
        XiInvtemp(:,:) = Xi_inv(:,:,b)
        temp_dwdx(:)= -dwdx(:,k)
        CdwdxB= MATMUL(XiInvtemp,temp_dwdx)

        do d= 1,SPH_dim
            delf= f0(b)-f0(a)
                call fncnApproxOperator(d1f(d,a),delf,mass(b),rho(b),CdwdxA(d))
             delf= f0(a)-f0(b)
                call fncnApproxOperator(d1f(d,b),delf,mass(a),rho(a),CdwdxB(d))
        enddo
        
     enddo
     
end
    
    