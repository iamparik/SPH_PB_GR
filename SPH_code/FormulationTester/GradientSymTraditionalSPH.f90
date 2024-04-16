!
!  SUBROUTINE: GradientSymTraditionalSPH
!
!  PURPOSE:  Subroutine that calls various SPH operators 
!           to find interpolations of different functions
!
!   CREATED:        06/08/2022       by  PARIKSHIT BOREGOWDA
!   Last Modified:  01/25/2023       by  PARIKSHIT BOREGOWDA 
!****************************************************************************
    
subroutine GradientSymTraditionalSPH(d1f,f0, SPH_dim, niac, pair_i, pair_j, dwdx, mass, rho, ntotal)
    implicit none
    integer(4), intent(in):: ntotal, SPH_dim, niac
    real(8), intent(in):: f0(ntotal), mass(ntotal), rho(ntotal)
    real(8), intent(in):: dwdx(SPH_dim,niac)
    integer(4), intent(in):: pair_i(niac), pair_j(niac) 
    real(8), intent(inout):: d1f(SPH_dim,ntotal)
    integer(4)::k,a,b,d
    real(8) ::  delf,temp_dwdx(SPH_dim),CdwdxA(SPH_dim),CdwdxB(SPH_dim)

          
    ! Call the interpolation operation for first order derivative of a function
    ! (ρ_a)∑_b〖(f_a/ρ_a**2 + f_b/ρ_b**2 ) ∂_j W_ab m_b 〗
    
     do k = 1,niac
        a=pair_i(k)
        b=pair_j(k)

        temp_dwdx(:)= dwdx(:,k)
        CdwdxA= temp_dwdx!/Gma(a)
        
        temp_dwdx(:)= -dwdx(:,k)
        CdwdxB= temp_dwdx!/Gma(b)
        
        do d= 1,SPH_dim
            delf= rho(a)*(f0(b)/(rho(b)**2)+f0(a)/(rho(a)**2))
            call fncnApproxOperator(d1f(d,a),delf,mass(b),1.D0,CdwdxA(d))
            delf= rho(b)*(f0(b)/(rho(b)**2)+f0(a)/(rho(a)**2))
            call fncnApproxOperator(d1f(d,b),delf,mass(a),1.D0,CdwdxB(d))
        enddo
        
        
     enddo
     
   
     
    end
    
    
    
    
    