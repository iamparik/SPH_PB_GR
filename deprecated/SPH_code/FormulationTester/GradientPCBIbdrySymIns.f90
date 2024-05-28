!
!  SUBROUTINE: GradientPCBIbdrySymIns
!			    This calculates PCBI type approximations of function's first order derivative at bdry
!               and symmetric/conservative form inside
!
!  PURPOSE:  Subroutine that calls various SPH operators 
!           to find interpolations of different functions
!
!   CREATED:        06/08/2022       by  PARIKSHIT BOREGOWDA
!   Last Modified:  06/18/2023       by  PARIKSHIT BOREGOWDA 
!****************************************************************************
    
subroutine GradientPCBIbdrySymIns(d1f,f0, fs0, Gma, Gma_inv, SPH_dim, niac, pair_i, pair_j, dwdx, eniac, epair_a, epair_s, dgmas, mass, rho, rhos, itype, ntotal, num_edges)
    implicit none
    integer(4), intent(in):: ntotal, SPH_dim, niac, eniac, num_edges
    integer(2), intent(in):: itype(ntotal)
    real(8), intent(in):: f0(ntotal), Gma(ntotal), Gma_inv(SPH_dim,SPH_dim,ntotal), mass(ntotal), rho(ntotal), rhos(num_edges), fs0(num_edges)
    real(8), intent(in):: dwdx(SPH_dim,niac), dgmas(SPH_dim,eniac) 
    integer(4), intent(in):: pair_i(niac), pair_j(niac) 
    integer(4), intent(in):: epair_a(eniac), epair_s(eniac)
    real(8), intent(inout):: d1f(SPH_dim,ntotal)
    integer(4)::k,a,b,s, d
    real(8) ::  delf,GmaInvtemp(SPH_dim,SPH_dim),temp_dwdx(SPH_dim),CdwdxA(SPH_dim),CdwdxB(SPH_dim)
    real(8) ::  temp_dgmas(SPH_dim),Cdgmas(SPH_dim)
    
          
    ! Call the interpolation operation for first order derivative of a function
    ! (ρ_a/𝛾_𝑎)∑_b〖(f_a/ρ_b**2 + f_b/ρ_b**2 ) ∂_j W_ab m_b 〗
    
     do k = 1,niac
        a=pair_i(k)
        b=pair_j(k)
        
        if(Gma(a) .eq. 1.D0) then
            temp_dwdx(:)= dwdx(:,k)
            CdwdxA= temp_dwdx/Gma(a)
            do d= 1,SPH_dim
                delf= rho(a)*(f0(b)/(rho(b)**2)+f0(a)/(rho(a)**2))
                call fncnApproxOperator(d1f(d,a),delf,mass(b),1.D0,CdwdxA(d))
            enddo
        else
            GmaInvtemp(:,:) = Gma_inv(:,:,a)
            temp_dwdx(:)= dwdx(:,k)
            CdwdxA= MATMUL(GmaInvtemp,temp_dwdx)
            do d= 1,SPH_dim
                delf= f0(b)-f0(a)
                call fncnApproxOperator(d1f(d,a),delf,mass(b),rho(b),CdwdxA(d))
            enddo
        endif
        
        
        if(Gma(b) .eq. 1.D0) then
            temp_dwdx(:)= -dwdx(:,k)
            CdwdxB= temp_dwdx/Gma(b)
            do d= 1,SPH_dim
                delf= rho(b)*(f0(b)/(rho(b)**2)+f0(a)/(rho(a)**2))
                call fncnApproxOperator(d1f(d,b),delf,mass(a),1.D0,CdwdxB(d))
            enddo
        else
            GmaInvtemp(:,:) = Gma_inv(:,:,b)
            temp_dwdx(:)= -dwdx(:,k)
            CdwdxB= MATMUL(GmaInvtemp,temp_dwdx)
            do d= 1,SPH_dim
                delf= f0(a)-f0(b)
                call fncnApproxOperator(d1f(d,b),delf,mass(a),rho(a),CdwdxB(d))
            enddo
            
        endif
        
     enddo
     
   
    ! Call the boundary itnegral terms of the first order derivative  of function interpolations
    ! [Γ_a^-1 ]_ij ∑_s∫_(∂(Ω ∩ Ω_w )_s〖(f_a-f_(s^' ) ) W_(as' ) n_(s_j ) dS' 〗
    do k= 1, eniac
        a=epair_a(k)
        s=epair_s(k)
        
        GmaInvtemp(:,:) = Gma_inv(:,:,a)
        temp_dgmas(:)= dgmas(:,k)
        Cdgmas=MATMUL(GmaInvtemp,temp_dgmas)
        
         do d=1,SPH_dim
             delf= f0(a)-fs0(s)
            call fncnApproxOperator(d1f(d,a),delf,1.D0,1.D0,Cdgmas(d))
        enddo
   
    enddo
     
    end
    
    
    
    
    