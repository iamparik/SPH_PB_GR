!
!  SUBROUTINE: DivergencePCBIBdry
!			    This calculates PCBI approximations of divergence of vector field near the boudnary  
!               but WCBI1 approx inside
!  PURPOSE:  Subroutine that calls various SPH operators 
!           to find interpolations of different functions
!
!   CREATED:        06/08/2022       by  PARIKSHIT BOREGOWDA
!   Last Modified:  12/04/2022       by  PARIKSHIT BOREGOWDA 
!****************************************************************************
    
subroutine DivergencePCBIBdry(d1f,f0, fs0, gamma_cont, Gma_inv, SPH_dim, niac, pair_i, pair_j, dwdx, eniac, epair_a, epair_s, dgmas, mass, rho, itype, ntotal, num_edges)
    implicit none
    integer(4), intent(in):: ntotal, SPH_dim, niac, eniac, num_edges
    integer(2), intent(in):: itype(ntotal)
    real(8), intent(in):: f0(SPH_dim,ntotal), gamma_cont(ntotal), Gma_inv(SPH_dim,SPH_dim,ntotal), mass(ntotal), rho(ntotal), fs0(SPH_dim,num_edges)
    real(8), intent(in):: dwdx(SPH_dim,niac), dgmas(SPH_dim,eniac) 
    integer(4), intent(in):: pair_i(niac), pair_j(niac) 
    integer(4), intent(in):: epair_a(eniac), epair_s(eniac)
    real(8), intent(inout):: d1f(ntotal)
    integer(4)::k,a,b,s, d
    real(8) ::  delf,GmaInvtemp(SPH_dim,SPH_dim),temp_dwdx(SPH_dim),CdwdxA(SPH_dim),CdwdxB(SPH_dim)
    real(8) ::  temp_dgmas(SPH_dim),Cdgmas(SPH_dim)
    
          
    ! Call the interpolation operation for first order derivative of a function
    ! [Γ_a^-1 ]_ij ∑_b〖(f_a-f_b ) ∂_j W_ab m_b/ρ_b 〗
    
     do k = 1,niac
        a=pair_i(k)
        b=pair_j(k)
        
        if(gamma_cont(a) .eq. 1.D0) then
            temp_dwdx(:)= dwdx(:,k)
            CdwdxA= temp_dwdx/gamma_cont(a) 
        else
            GmaInvtemp(:,:) = Gma_inv(:,:,a)
            temp_dwdx(:) = dwdx(:,k)
            CdwdxA= MATMUL(GmaInvtemp,temp_dwdx)
        endif
        
        if(gamma_cont(b) .eq. 1.D0) then
            temp_dwdx(:)= -dwdx(:,k)
            CdwdxB= temp_dwdx/gamma_cont(b)
            
        else
            GmaInvtemp(:,:) = Gma_inv(:,:,b)
            temp_dwdx(:) = -dwdx(:,k)
            CdwdxB= MATMUL(GmaInvtemp,temp_dwdx)          
        endif

        do d= 1,SPH_dim
            delf= f0(d,b)-f0(d,a)
            call fncnApproxOperator(d1f(a),delf,mass(b),rho(b),CdwdxA(d))
            delf= f0(d,a)-f0(d,b)
            call fncnApproxOperator(d1f(b),delf,mass(a),rho(a),CdwdxB(d))
        enddo
        
        
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
             delf= f0(d,a)-fs0(d,s)
            call fncnApproxOperator(d1f(a),delf,1.D0,1.D0,Cdgmas(d))
        enddo
   
    enddo
     
    end
    
    
    
    
    