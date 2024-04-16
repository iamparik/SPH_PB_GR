!
!  SUBROUTINE: LaplacianSPHTraditional
!               ⟨⟨𝜵∙ 𝜵𝑓_𝑎 ⟩⟩= −2/𝛾_𝑎  ∑_𝑏〖(𝑓_𝑎 − 𝑓_𝑏)((𝒙_𝑎 − 𝒙_b ))/(𝑟_𝑎𝑏^2 )∙𝛁𝑊 𝑚_𝑏/𝜌_𝑏 〗
!                                                   + 2/𝛾_𝑎 ∑_𝑠〖(𝛁𝑓_𝑎)∙𝛁𝛾_𝑎𝑠〗 
!
!  PURPOSE:  Subroutine that calls various SPH operators 
!           to find interpolations of different functions
!
!   CREATED:        09/02/2022       by  PARIKSHIT BOREGOWDA
!   Last Modified:  01/09/2023       by  PARIKSHIT BOREGOWDA 
!****************************************************************************
    
subroutine LaplacianSPHTraditionalCorrected(lapf, f0, G, x, Gma_inv, SPH_dim, niac, pair_i, pair_j, dwdx, mass, rho, itype, ntotal)
    implicit none
    integer(4), intent(in):: ntotal, SPH_dim, niac
    integer(2), intent(in):: itype(ntotal)
    real(8), intent(in):: f0(ntotal), Gma_inv(SPH_dim,SPH_dim,ntotal), mass(ntotal), rho(ntotal)
    real(8), intent(in):: x(SPH_dim,ntotal), G(ntotal)
    real(8), intent(in):: dwdx(SPH_dim,niac)
    integer(4), intent(in):: pair_i(niac), pair_j(niac)
    real(8), intent(inout):: lapf(ntotal)
    integer(4)::k,a,b,s,d
    real(8) ::  GmaInvtemp(SPH_dim,SPH_dim),temp_dwdx(SPH_dim),CdwdxA(SPH_dim),CdwdxB(SPH_dim)
    real(8) :: r_ab, lap_coeff
    
          
    ! Call the interpolation operation for finding the Laplacian of a function
    ! −1/𝛾_𝑎  ∑_𝑏〖(𝑘_𝑎+𝑘_𝑏)(𝑓_𝑎 − 𝑓_𝑏)((𝒙_𝑎 − 𝒙_b ))/(𝑟_𝑎𝑏^2 )∙𝛁𝑊 𝑚_𝑏/𝜌_𝑏 〗
     do k = 1,niac
        a=pair_i(k)
        b=pair_j(k)
        
        GmaInvtemp(:,:) = Gma_inv(:,:,a)
        temp_dwdx(:)= dwdx(:,k)
        CdwdxA= MATMUL(GmaInvtemp,temp_dwdx)
        
        GmaInvtemp(:,:) = Gma_inv(:,:,b)
        temp_dwdx(:)= -dwdx(:,k)
        CdwdxB= MATMUL(GmaInvtemp,temp_dwdx)
        
        r_ab=norm2(x(:,a)-x(:,b))
        if ( .not. (isNAN(1.D0/r_ab))) then
            do d=1,SPH_dim
            
                lap_coeff=(G(a)+G(b))*(f0(a)-f0(b))*(x(d,a)-x(d,b))/(r_ab**2)
                call fncnApproxOperator(lapf(a),lap_coeff,mass(b),rho(b),CdwdxA(d))
                
                lap_coeff=(G(a)+G(b))*(f0(b)-f0(a))*(x(d,b)-x(d,a))/(r_ab**2)
                call fncnApproxOperator(lapf(b),lap_coeff,mass(a),rho(a),CdwdxB(d))
            enddo
        endif

     enddo
     
     
end



